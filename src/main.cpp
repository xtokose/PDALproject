#include <pdal/StageFactory.hpp>
#include <pdal/Options.hpp>
#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>
#include <pdal/Dimension.hpp>
#include <pdal/Stage.hpp>
#include <pdal/io/BufferReader.hpp>
#include <pdal/PointRef.hpp>
#include <pdal/filters/OutlierFilter.hpp>
#include <pdal/filters/RangeFilter.hpp>
#include <pdal/Filter.hpp>
#include <pdal/filters/NormalFilter.hpp>
#include <pdal/filters/DBSCANFilter.hpp>
#include <pdal/PipelineManager.hpp>
#include <pdal/filters/DelaunayFilter.hpp>
#include <pdal/Mesh.hpp> 
#include <pdal/io/PlyWriter.hpp>

//#include <pdal/PointViewId.hpp>
#include <pdal/DimUtil.hpp>

#include <Eigen/Dense>
#include <stdexcept>
#include <cstddef>




#include <cstdlib>
#include <iomanip>

#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <unordered_set>
#include <thread>
#include <map>
#include <limits>

const double PI = 3.14159265359;
using namespace std;
using namespace pdal;
using namespace Eigen;

struct BestFitPlane
{
    
    Eigen::Vector3d normal;   // unit-length (A,B,C)
    double d;                 // D in Ax + By + Cz + D = 0
    Eigen::Vector3d centroid; // point on plane
};
struct Line3D {
    Vector3d p;   // point on line
    Vector3d dir; // direction (not necessarily normalized)
    //vector<int> origin;

    //vector<int> planeId;
    Vector2d origin_planes;
};
struct mEdge {
    PointId id1;
    PointId id2;

    mEdge(PointId i, PointId j)
        : id1(i), id2(j)
    {
    }
    bool operator<(const mEdge& other) const {
        if (id1 != other.id1) return id1 < other.id1;
        return id2 < other.id2;
    }
};





struct MeshReadyView
{
    std::shared_ptr<pdal::PointTable> table; // must stay alive with the view
    pdal::PointViewPtr view;
};
struct RoofData {
    PointTable table;
    PointViewSet view; // contains the PointViewPtrs made from `table`


    Vector3d centroid;
};
struct MainRoof {
    unique_ptr<RoofData> roof;


    vector<unique_ptr<RoofData>> subroof;
    vector<unique_ptr<RoofData>> msubroof;  //modified subroofs (points are normal vectors)
    vector<unique_ptr<RoofData>> nsubroof;  //subroofs clustered by normal vectors
    vector<unique_ptr<RoofData>> fsubroof;



    vector<BestFitPlane> plane;
	vector<Line3D> lines;

    vector<TriangularMesh*> initial_mesh;
    vector<TriangularMesh> new_mesh;
    vector<MeshReadyView> mesh_ready_view;
    
    vector<vector<Vector3d>> jacket;


};
class WallFilter : public pdal::Filter
{
public:
    std::string getName() const override { return "filters.walls"; }

private:
    double m_curvature = 0.0;
    double m_angle = 0.0;
    double ratio;
    vector<double> excludedZ;

    void addArgs(pdal::ProgramArgs& args) override
    {
        args.add("curvature", "Curvature ratio threshold", m_curvature);
        args.add("angle", "Angle ratio threshold", m_angle);
    }

    pdal::PointViewSet run(pdal::PointViewPtr in) override
    {
        auto start = chrono::high_resolution_clock::now();
        pdal::PointViewPtr out(new pdal::PointView(in->table()));
        ratio = sin(m_angle * PI / 180.0);
        
        for (pdal::PointId i = 0; i < in->size(); ++i)
        {
            if (keep(*in, i))                 // your criteria
                out->appendPoint(*in, i);     // copies all dimensions for point i
        }

        //debug
        auto duration = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start);
        //cout << "Points remaining: " << get_points_size(buildpoint) << endl;
        double min = *min_element(excludedZ.begin(), excludedZ.end());
        double max = *max_element(excludedZ.begin(), excludedZ.end());
        double average = accumulate(excludedZ.begin(), excludedZ.end(), 0.0) / excludedZ.size();
        /*cout << "min " << min << endl;
        cout << "max " << max << endl;
        cout << "average " << average << endl;*/




        pdal::PointViewSet s;
        s.insert(out);
        return s;
    }
    bool keep(const pdal::PointView& v, pdal::PointId i)
    {
        
        const double nz = v.getFieldAs<double>(pdal::Dimension::Id::NormalZ, i);
        const double curv = v.getFieldAs<double>(pdal::Dimension::Id::Curvature, i);


        if (nz > ratio) {
            /*if (curv < m_curvature) {
                return true;
            }
            else {
				excludedZ.push_back(abs(nz));
            }*/

            return true;
        }
        else {
            excludedZ.push_back(abs(nz));
            return false;
        }
    }
};



class PointReader {
private:
    bool print_progress = false;
    StageFactory factory;
    PipelineManager manager;
    PointTable mainTable;


    size_t number_of_roofs;
    size_t get_points_size(PointViewSet& set) {
        size_t size = 0;
        for (auto v : set) size += v->size();
        return size;
    }
    void recomputeNormals() {
        ////not using anywhere currently
        ////not using anywhere currently

        PointViewPtr inputView = *buildpoint.begin();

        BufferReader bufferReader;
        bufferReader.addView(inputView);

        PointTable table;

        // Create the normal filter
        NormalFilter normalFilter;
        Options options;
        options.add("knn", 16);          // typical K-nearest neighbors
        normalFilter.setOptions(options);
        normalFilter.setInput(bufferReader);
        normalFilter.prepare(table);
        buildpoint = normalFilter.execute(table);


    }
    unordered_set<double> getRoofsIDs(PointViewPtr view){


        unordered_set<double> set;
        for (PointId i = 0; i < view->size(); ++i)
        {
            double id = view->getFieldAs<double>(Dimension::Id::ClusterID, i);


            if (id > 0)               // ignore noise
                set.insert(id);
        }
        return set;
    }


    double distanceToLine(const Line3D& line, const Eigen::Vector3d& point)
    {
        double dirLen = line.dir.norm();
        if (dirLen == 0.0)
            throw std::runtime_error("Line direction is zero.");

        return (point - line.p).cross(line.dir).norm() / dirLen;
    }


public:
    PointViewSet buildpoint;
    vector<MainRoof> roofs;



    void printSchema(PointTable table)
    {
        auto layout = table.layout();
        cout << "Dimensions (" << layout->dims().size() << "):\n";
        for (auto id : layout->dims())
            cout << "  - " << Dimension::name(id) << "\n";
    }
    void printAllDimensions(size_t N = 10)
    {
        PointViewPtr v = *buildpoint.begin();
        auto layout = mainTable.layout();
        const PointId n = min<PointId>(N, v->size());

        for (PointId i = 0; i < n; ++i)
        {
            cout << "Point " << i << ":\n";
            for (auto id : layout->dims())
            {
                // getFieldAs<double> is convenient for numeric dims (applies scale/offset)
                double val = v->getFieldAs<double>(id, i);
                cout << "  " << left << setw(18)
                    << Dimension::name(id) << " = " << val << "\n";
            }
        }
    }


        //STAGES, all points
    Stage* loadFile(string filepath) {
        auto start = chrono::high_resolution_clock::now();

        Stage* reader = factory.createStage("readers.las");
        Options opts;                                        ////options
        opts.add("filename", filepath);
        reader->setOptions(opts);


        PointTable table;
        reader->prepare(table);
        PointViewSet all_points = reader->execute(table);

            //debug
        auto duration = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start);
        if(print_progress){
            cout << "\nLoading finished in " << duration.count()/1000. << " seconds\n";
            cout << "Points remaining: " << get_points_size(all_points) << endl;
        }
        return reader;       
    }
    Stage* filterClass6(Stage* input, string out_file) {
        auto start = chrono::high_resolution_clock::now();

        Stage* range = factory.createStage("filters.range");
        Options opts;
        opts.add("limits", "Classification[6:6]");
        range->setOptions(opts);
        range->setInput(*input);

        


            //debug
        auto duration = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start);
        cout << "\n-- Filtering class 6 finished in " << duration.count() / 1000. << " miliseconds\n";
        string path = out_file + "/points_after_filtering_class_6.txt";
        if(print_progress){ writeToTextFile(*range, path); }
    

        return range;
    }
    Stage* computeNormals(Stage* input, size_t nearest_neighbors) {
        auto start = chrono::high_resolution_clock::now();

        Stage* normal = factory.createStage("filters.normal");
        PointTable table;
        Options opts;
        opts.add("knn", nearest_neighbors);                         // number of neighbors
        normal->setOptions(opts);
        normal->setInput(*input);        /////////input
 

            //debug      
        auto duration = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start);
        cout << "\n-- Computing normals finished in " << duration.count() / 1000. << " miliseconds\n";
        return normal;
    }
    Stage* filterWalls(Stage* input, string out_file, double ratio_curvature, double ratio_angle) {
        auto start = chrono::high_resolution_clock::now();

        WallFilter* my = new WallFilter();
        Options o;
        o.add("curvature", ratio_curvature);
        o.add("angle", ratio_angle);
        my->setOptions(o);
        my->setInput(*input);


        //debug
        auto duration = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start);
        cout << "\n-- Wall filtering finished in " << duration.count() / 1000. << " miliseconds\n";
        string path = out_file + "/points_after_wall_filter.txt";
        if (print_progress){ writeToTextFile(*my, path); }

        return my;
    }
    Stage* filterOutliers(Stage* input, string out_file, string method, size_t mean_k, float multiplier, float radius, size_t min_k) {
        auto start = chrono::high_resolution_clock::now();



        Stage& outlier = manager.makeFilter("filters.outlier", *input);

        Options opts;
        if (method == "statistical") {
            opts.add("method", "statistical");
            opts.add("mean_k", mean_k);                      //mean number of neighbors  6
            opts.add("multiplier", multiplier);                //Standard deviation threshold 0.5
            outlier.setOptions(opts);
        }
        else if (method == "radius") {
            opts.add("method", "radius");
            opts.add("radius", radius);                     //1
            opts.add("min_k", min_k);                       //min number of neighbors in radius  5
        }



        Stage& rangeb = manager.makeFilter("filters.range", outlier);
        Options rangeOpts;
        rangeOpts.add("limits", "Classification![7:7]");
        rangeb.setOptions(rangeOpts);

        //debug
        auto duration = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start);
        cout << "\n-- Outlier filtering finished in " << duration.count() / 1000. << " miliseconds\n";
        string path = out_file + "/points_after_outliers_filter.txt";
        if (print_progress) { writeToTextFile(rangeb, path); }
        return &rangeb;
    }
    Stage* clusterPoints(Stage* input, float tolerance, size_t min_points) {
        auto start = chrono::high_resolution_clock::now();
        Stage* cluster = factory.createStage("filters.cluster");
        Options opts;
        opts.add("tolerance", tolerance);                     //max distance point to be added to the cluster
        opts.add("min_points", min_points);                    //minimum number of points in a cluster
        cluster->setOptions(opts);
        cluster->setInput(*input);
        



      
        //debug      
        auto duration = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start);
        cout << "\n-- Clustering finished in " << duration.count() / 1000. << " miliseconds\n";
        return cluster;
    }
    Stage* zsmooth(Stage* input, string out_file, size_t radius) {
        auto start = chrono::high_resolution_clock::now();

        Options zsmoothOptions;
        zsmoothOptions.add("radius", radius);
        zsmoothOptions.add("dim", "Zsmooth");

        Stage* zsmooth = factory.createStage("filters.zsmooth");
        zsmooth->setInput(*input);
        zsmooth->setOptions(zsmoothOptions);


        Options opts;
        opts.add("value", "Z = Zsmooth");
        Stage* assign = factory.createStage("filters.assign");
        assign->setInput(*zsmooth);
        assign->setOptions(opts);


        //debug      
        auto duration = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start);
        cout << "\n-- Smoothing finished in " << duration.count() / 1000. << " miliseconds\n";
        if (print_progress){ writeToTextFile(*assign, out_file + "/points_after_smoothing.txt"); }

        return assign;
    }

    Stage* addTextWriter(Stage& input, const string& filename)
    {
        Stage* w = factory.createStage("writers.text");

        Options o;
        o.add("filename", filename);   // or "STDOUT"
        o.add("format", "csv");       // "csv" or "geojson"
        o.add("order", "X,Y,Z");         // e.g. "X,Y,Z,Intensity,ClusterID"
        o.add("keep_unspecified", false);   // <-- only write the dims in `order`
        o.add("write_header", false);
        o.add("delimiter", ",");
        o.add("precision", 2);


        w->setOptions(o);
        w->setInput(input);
        return w;
    }
    void writeToTextFile(Stage& input, const string& filename) {
        
            //writing
        auto writer = addTextWriter(input, filename);
        PointTable table;
        writer->prepare(table);
        writer->execute(table);


            //counting the lines
        ifstream in(filename);
        string line;
        size_t lines = 0;
        while (getline(in, line)) {
            ++lines;
        }
        
		cout << "Points written to " << filename << endl;
        cout << "Points remaining " << lines << endl;
	}
    void execute(Stage* input) {
        auto start2 = chrono::high_resolution_clock::now();

        input->prepare(mainTable);
        buildpoint = input->execute(mainTable);


        unordered_set<double> clusterIDs = getRoofsIDs(*buildpoint.begin());
        number_of_roofs = clusterIDs.size();
        

            //debug
        auto duration = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start2);
        cout << "\Executing point table finished in " << duration.count() / 1000. << " seconds\n";
        cout << "Point cloud has been clustered into " << number_of_roofs << " roofs" << endl;
    }
    void filerByZValue() {
        //auto start = chrono::high_resolution_clock::now();
        //PointViewPtr view = *buildpoint.begin();
        //PointViewSet temporaryset;
        //PointViewPtr outView = view->makeNew();

        //

        //for (PointId i = 0; i < view->size(); ++i) {
        //    bool append = true;
        //    for (PointId j = 0; j < view->size(); ++j) {
        //        double xi = view->getFieldAs<double>(Dimension::Id::X, i);
        //        double yi = view->getFieldAs<double>(Dimension::Id::Y, i);
        //        double zi = view->getFieldAs<double>(Dimension::Id::Z, i);

        //        double xj = view->getFieldAs<double>(Dimension::Id::X, j);
        //        double yj = view->getFieldAs<double>(Dimension::Id::Y, j);
        //        double zj = view->getFieldAs<double>(Dimension::Id::Z, j);

        //        
        //        if (((abs(xi - xj) < ratio_Zvalue) &&
        //            (abs(yi - yj) < ratio_Zvalue) &&
        //            (zj > zi))) {
        //            j = view->size();
        //            append = false;
        //        }
        //    }
        //    if (append) {
        //        outView->appendPoint(*view, i);
        //    }
        //    
        //}

        //temporaryset.insert(outView);
        //buildpoint = temporaryset;

        ////debug      
        //if (debug) {
        //    auto duration = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start);
        //    cout << "\nFiltering by Z value finished in " << duration.count() / 1000. << " seconds\n";
        //    cout << "Points remaining: " << get_points_size(buildpoint) << endl;

        //}
        //makePointFile("points_after_Zfilter.txt");
    }
 


    void makePointViewSetRoofs() {
        PointViewPtr old_view = *buildpoint.begin();
        cout << "Making pointViewSet   roofs " << endl;


        for (int k = 0; k < number_of_roofs; k++) {
            auto rd = make_unique<RoofData>();

            PointLayoutPtr layout = rd->table.layout();
            layout->registerDim(Dimension::Id::X);
            layout->registerDim(Dimension::Id::Y);
            layout->registerDim(Dimension::Id::Z);
            layout->registerDim(Dimension::Id::NormalX);
            layout->registerDim(Dimension::Id::NormalY);
            layout->registerDim(Dimension::Id::NormalZ);
            layout->registerDim(Dimension::Id::Curvature);
            layout->registerDim(Dimension::Id::ClusterID);

            int count = 0;
            PointViewPtr view(new PointView(rd->table));
            for (int j = 0; j < old_view->size(); j++) {

                size_t clusterID = old_view->getFieldAs<double>(Dimension::Id::ClusterID, j);
                if (clusterID == (k + 1)) {


                    double x = old_view->getFieldAs<double>(Dimension::Id::X, j);
                    double y = old_view->getFieldAs<double>(Dimension::Id::Y, j);
                    double z = old_view->getFieldAs<double>(Dimension::Id::Z, j);
                    double nx = old_view->getFieldAs<double>(Dimension::Id::NormalX, j);
                    double ny = old_view->getFieldAs<double>(Dimension::Id::NormalY, j);
                    double nz = old_view->getFieldAs<double>(Dimension::Id::NormalZ, j);
                    double curv = old_view->getFieldAs<double>(Dimension::Id::Curvature, j);


                    view->setField(Dimension::Id::X, count, x);
                    view->setField(Dimension::Id::Y, count, y);
                    view->setField(Dimension::Id::Z, count, z);
                    view->setField(Dimension::Id::NormalX, count, nx);
                    view->setField(Dimension::Id::NormalY, count, ny);
                    view->setField(Dimension::Id::NormalZ, count, nz);
                    view->setField(Dimension::Id::Curvature, count, curv);
                    view->setField(Dimension::Id::ClusterID, count, 0);

                    count++;
                }

            }
            cout << "Roof " << (k + 1) << " --- contains " << view->size() << " points " << endl;

            rd->view.insert(view);
            MainRoof r;
            r.roof = move(rd);
            roofs.push_back(move(r));
        }

    }





    void printMainRoofs(const PointViewPtr& view, string out_file, string folder) {
     
        //getting number of clusters
        unordered_set<double> clusterIDs = getRoofsIDs(view);
        double numClusters = clusterIDs.size();


            //creating empty files
        string path = out_file + "/" + folder;
        filesystem::create_directories(path);
        vector<ofstream> files;
        files.reserve(numClusters);
        for (int i = 0; i < numClusters; i++) {
            files.emplace_back(path + "/roof" + to_string(i + 1) + ".txt");
            files[i] << fixed << setprecision(3);
        }
  
            //fillilng files
        for (PointId i = 0; i < view->size(); ++i) {
            double id = view->getFieldAs<double>(Dimension::Id::ClusterID, i);
            double x = view->getFieldAs<double>(Dimension::Id::X, i);
            double y = view->getFieldAs<double>(Dimension::Id::Y, i);
            double z = view->getFieldAs<double>(Dimension::Id::Z, i);
            


            if ((id > 0) && (id < 10000)) {
                files[id - 1] << x << "," << y << "," << z << endl;
            }
        }

            //closing files
        for (int i = 0; i < numClusters; i++) {
            files[i].close();
        }

            //debug
        cout << numClusters << " files have been created" << endl;
    }





    using RoofVec = std::vector<std::unique_ptr<RoofData>>;


    void createNewCloud(PointViewPtr input, RoofVec& out) {
        //PointViewPtr input = *roofs[i].roof->view.begin();


            //getting number of clusters
        unordered_set<double> clusterIDs = getRoofsIDs(input);
        double numClusters = clusterIDs.size();


        for (int k = 0; k < numClusters; k++) {

            auto rd = std::make_unique<RoofData>();
            PointLayoutPtr layout = rd->table.layout();

            layout->registerDim(Dimension::Id::X);
            layout->registerDim(Dimension::Id::Y);
            layout->registerDim(Dimension::Id::Z);
            layout->registerDim(Dimension::Id::NormalX);
            layout->registerDim(Dimension::Id::NormalY);
            layout->registerDim(Dimension::Id::NormalZ);
            layout->registerDim(Dimension::Id::Curvature);
            layout->registerDim(Dimension::Id::ClusterID);

            int count = 0;
            PointViewPtr view(new PointView(rd->table));
            for (int j = 0; j < input->size(); j++) {

                size_t clusterID = input->getFieldAs<double>(Dimension::Id::ClusterID, j);
                if (clusterID == (k + 1)) {


                    double x = input->getFieldAs<double>(Dimension::Id::X, j);
                    double y = input->getFieldAs<double>(Dimension::Id::Y, j);
                    double z = input->getFieldAs<double>(Dimension::Id::Z, j);
                    double nx = input->getFieldAs<double>(Dimension::Id::NormalX, j);
                    double ny = input->getFieldAs<double>(Dimension::Id::NormalY, j);
                    double nz = input->getFieldAs<double>(Dimension::Id::NormalZ, j);
                    double curv = input->getFieldAs<double>(Dimension::Id::Curvature, j);


                    view->setField(Dimension::Id::X, count, x);
                    view->setField(Dimension::Id::Y, count, y);
                    view->setField(Dimension::Id::Z, count, z);
                    view->setField(Dimension::Id::NormalX, count, nx);
                    view->setField(Dimension::Id::NormalY, count, ny);
                    view->setField(Dimension::Id::NormalZ, count, nz);
                    view->setField(Dimension::Id::Curvature, count, curv);
                    view->setField(Dimension::Id::ClusterID, count, 0);

                    count++;
                }
            }

            rd->view.insert(view);
            //roofs[i].subroof.push_back(move(rd));
            out.push_back(move(rd));

        }

        //return rd;

    }


    void fillSubroofs() {


        for (int i = 0; i < number_of_roofs; i++) {

            PointViewPtr input = *roofs[i].roof->view.begin();
            createNewCloud(input, roofs[i].subroof);
        }

        
    }
    void endSubroofs() {

        for (int i = 0; i < number_of_roofs; i++) {

            for (int j = 0; j < roofs[i].nsubroof.size(); j++) {

                PointViewPtr input = *roofs[i].nsubroof[j]->view.begin();
                createNewCloud(input, roofs[i].fsubroof);
            }

        }

    }

    void modifySubroofs() {

        for (int i = 0; i < number_of_roofs; i++) {


            for (int k = 0; k < roofs[i].subroof.size(); k++) {
                PointViewPtr subroof = *roofs[i].subroof[k]->view.begin();

                auto rd = std::make_unique<RoofData>();
                PointLayoutPtr layout = rd->table.layout();

                layout->registerDim(Dimension::Id::X);
                layout->registerDim(Dimension::Id::Y);
                layout->registerDim(Dimension::Id::Z);
                layout->registerDim(Dimension::Id::Curvature);
                layout->registerDim(Dimension::Id::ClusterID);


                PointViewPtr view(new PointView(rd->table));
                for (int j = 0; j < subroof->size(); j++) {

                    double nx = subroof->getFieldAs<double>(Dimension::Id::NormalX, j);
                    double ny = subroof->getFieldAs<double>(Dimension::Id::NormalY, j);
                    double nz = subroof->getFieldAs<double>(Dimension::Id::NormalZ, j);

                    view->setField(Dimension::Id::X, j, nx);
                    view->setField(Dimension::Id::Y, j, ny);
                    view->setField(Dimension::Id::Z, j, nz);
                }

                rd->view.insert(view);
                roofs[i].msubroof.push_back(move(rd));

            }

            

        }


    }
    void appointClusterID() {

        /*for (int i = 0; i < roofs.size(); i++) {

            for (int k = 0; k < roofs[i].subroof.size(); k++) {
                
                PointViewPtr view = *roofs[i].subroof[k]->view.begin();
                *roofs[i].subroof[k]->view.begin()
            }

        }*/
    }

    void finalizeSubroofs() {

        for (int i = 0; i < roofs.size(); i++) {

            for (int k = 0; k < roofs[i].msubroof.size(); k++) {
                
                PointViewPtr msubroof = *roofs[i].msubroof[k]->view.begin();
                PointViewPtr subroof = *roofs[i].subroof[k]->view.begin();

                //getting number of clusters
                unordered_set<double> clusterIDs = getRoofsIDs(msubroof);
                double numClusters = clusterIDs.size();




                for (int l = 0; l < numClusters; l++) {

                    auto rd = std::make_unique<RoofData>();
                    PointLayoutPtr layout = rd->table.layout();
                    layout->registerDim(Dimension::Id::X);
                    layout->registerDim(Dimension::Id::Y);
                    layout->registerDim(Dimension::Id::Z);
                    layout->registerDim(Dimension::Id::NormalX);
                    layout->registerDim(Dimension::Id::NormalY);
                    layout->registerDim(Dimension::Id::NormalZ);
                    layout->registerDim(Dimension::Id::Curvature);
                    layout->registerDim(Dimension::Id::ClusterID);
                    PointViewPtr view(new PointView(rd->table));


                    int count = 0;
                    for (int m = 0; m < msubroof->size(); m++) {

                        
                        size_t clusterID = msubroof->getFieldAs<double>(Dimension::Id::ClusterID, m);
                        if (clusterID == (l + 1)) {

                            double x = subroof->getFieldAs<double>(Dimension::Id::X, m);
                            double y = subroof->getFieldAs<double>(Dimension::Id::Y, m);
                            double z = subroof->getFieldAs<double>(Dimension::Id::Z, m);
                            double nx = subroof->getFieldAs<double>(Dimension::Id::NormalX, m);
                            double ny = subroof->getFieldAs<double>(Dimension::Id::NormalY, m);
                            double nz = subroof->getFieldAs<double>(Dimension::Id::NormalZ, m);
                            double curv = subroof->getFieldAs<double>(Dimension::Id::Curvature, m);



                            view->setField(Dimension::Id::X, count, x);
                            view->setField(Dimension::Id::Y, count, y);
                            view->setField(Dimension::Id::Z, count, z);
                            view->setField(Dimension::Id::NormalX, count, nx);
                            view->setField(Dimension::Id::NormalY, count, ny);
                            view->setField(Dimension::Id::NormalZ, count, nz);
                            view->setField(Dimension::Id::Curvature, count, curv);
                            view->setField(Dimension::Id::ClusterID, count, 0);
                            count++;
                        }
                    }
                    rd->view.insert(view);
                    roofs[i].nsubroof.push_back(move(rd));


                }

                


            }
        }



        for (int i = 0; i < roofs.size(); i++) {


            cout << "oroof " << i + 1 << " contains " << roofs[i].subroof.size() << " subroofs " << endl;
            cout << "froof " << i + 1 << " contains " << roofs[i].nsubroof.size() << " subroofs " << endl;
            cout << " " << endl;

        }



    }

    void clusterByPoints(string rooftype, float tolerance, size_t min_points) {



        cout << "Point clustering begins" << endl;

       
        if (rooftype == "mainroofs") {
            vector<RoofData*> tt;
            for (int i = 0; i < roofs.size(); i++) {
                tt.push_back(roofs[i].roof.get());
            }
            cluster(tt, tolerance, min_points);
        }
        else if (rooftype == "nroofs") {

            for (int i = 0; i < roofs.size(); i++) {

                vector<RoofData*> tt;
                for (int j = 0; j < roofs[i].nsubroof.size(); j++) {
                    cout << "j  " << j << endl;
                    tt.push_back(roofs[i].nsubroof[j].get());

                }
                cout << "DONE 1 " << endl;
                cluster(tt, tolerance, min_points);
                cout << "DONW 2 " << endl;
            }

        }

        


        //for (auto& r : roofs) {
        //    if (r.roof) {
        //        tt.push_back(r.roof.get());   // <-- NO move
        //    }
        //}

        
    }
    void clusterByNormals(float tolerance, size_t min_points) {

       
        
        cout << "Normal clustering begins " << endl;
        for (int i = 0; i < roofs.size(); i++) {

            cout << "clustering roof " << i << endl;
            vector<RoofData*> tt;
            for (int j = 0; j < roofs[i].msubroof.size(); j++) {
                tt.push_back(roofs[i].msubroof[j].get());

            }
            cluster(tt, tolerance, min_points);

        }
     
        

    }
    void cluster(const vector<RoofData*>& roof, float tolerance, size_t min_points) {
 

        vector<thread> threads;
        
        Options opts;
        opts.add("tolerance", tolerance);                     //max distance point to be added to the cluster
        opts.add("min_points", min_points);                    //minimum number of points in a cluster
        for (int i = 0; i < roof.size(); ++i) {
            //threads.emplace_back(worker, i, cluster, opts, roof); // starts immediately
            threads.emplace_back(&PointReader::worker, this, (int)i, opts, std::cref(roof));
        }


        for (auto& t : threads) {
            t.join(); // wait for all threads to finish
        }
    }
    void worker(int id, Options opts, const std::vector<RoofData*>& roof) {
        Stage* stage = factory.createStage("filters.cluster");
        auto& rd = *roof[id];
        BufferReader reader;
        reader.addView(*rd.view.begin());

        stage->setOptions(opts);
        stage->setInput(reader);
        stage->prepare(rd.table);
        stage->execute(rd.table);


        //std::cout << "Hello from thread " << id << "\n";
    }

    void makeClusteredRoofsFiles(string folder) {
        std::error_code ec;

        // Relative path: this is relative to the *current working directory*
        filesystem::path dir = folder;

        ec.clear();


        auto n = filesystem::remove_all(dir, ec); // deletes directory + all contents
        if (ec) {
            std::cerr << "remove_all failed: " << ec.message() << "\n";
            return;
        }
        filesystem::create_directories(folder);

        /*string folderName = "folderujem";
        if (!filesystem::exists(folderName)) {
            if (filesystem::create_directories(folderName)) {
                std::cout << "Folder created: " << folderName << '\n';
            }
            else {
                std::cout << "Failed to create folder\n";
            }
        }*/

        cout << "whaat " << endl;


        //std::error_code ec;
        //if (filesystem::exists(folder, ec))
        {
            //filesystem::remove_all(folder, ec);   // recursive
            //if (ec) throw filesystem::filesystem_error("remove_all failed", folder, ec);
        }

        
        //filesystem::create_directories(folder, ec);
        //if (ec) throw filesystem::filesystem_error("create_directories failed", folder, ec);



        /*filesystem::path path = folder;
        for (const auto& entry : filesystem::directory_iterator(path)) {
            if (!filesystem::is_regular_file(entry.status()))
                filesystem::remove(entry.path());

            if (entry.path().extension() == ".txt") {
                filesystem::remove(entry.path());
            }
        }*/



        size_t count = 1;
        for (int i = 0; i < roofs.size(); i++) {
            string folderName = folder + "/roof" + to_string(count);
            cout << "roof " << to_string(count);

            auto& rd = *roofs[i].roof;
            PointViewPtr ff = *rd.view.begin();
            cout << " cluster size  " << ff->size() << endl;

            if (ff->size() != 0) {
                //makeClusteredFiles(*rd.view.begin(), folderName);
                count++;
            }

            

        }






    }
    
    void printFinalRoofs(string folder) {




        error_code ec;
        filesystem::path dir = folder;
        ec.clear();

        filesystem::create_directories(folder);
        auto n = filesystem::remove_all(dir, ec); // deletes directory + all contents
        if (ec) {
            std::cerr << "remove_all failed: " << ec.message() << "\n";
            return;
        }
        filesystem::create_directories(folder);
        cout << "done " << endl;



        for (int i = 0; i < roofs.size(); i++) {

            string path = folder + "/roof" + to_string(i + 1);
 
            filesystem::create_directories(path);

            for (int j = 0; j < roofs[i].nsubroof.size(); j++) {

                string fileName = path + "/subroof" + to_string(j + 1) + ".txt";
                size_t n = printCloud(roofs[i].nsubroof[j]->view, fileName);
                cout << "roof " << i << "subroof " << j << " ----- " << n << endl;

                /*PointViewPtr inputView = *roofs[i].fsubroof[j]->view.begin();
                std::cerr
                    << "has X=" << inputView->hasDim(pdal::Dimension::Id::X)
                    << " Y=" << inputView->hasDim(pdal::Dimension::Id::Y)
                    << " Z=" << inputView->hasDim(pdal::Dimension::Id::Z)
                    << "\n";
                double x = inputView->getFieldAs<double>(Dimension::Id::X, 0);
                double y = inputView->getFieldAs<double>(Dimension::Id::Y, 0);
                double z = inputView->getFieldAs<double>(Dimension::Id::Z, 0);
                cout << "x " << x << " " << y << " " << z << endl;


                cout << "1 " << endl;
                BufferReader reader;
                reader.addView(inputView);
                

                string name = path + "/roof" + to_string(j + 1);
                writeToTextFile(reader, name);*/

            }

        }
    }

    void printFinalRoofs22(string folder) {




        error_code ec;
        filesystem::path dir = folder;
        ec.clear();

        filesystem::create_directories(folder);
        auto n = filesystem::remove_all(dir, ec); // deletes directory + all contents
        if (ec) {
            std::cerr << "remove_all failed: " << ec.message() << "\n";
            return;
        }
        filesystem::create_directories(folder);
        cout << "done " << endl;



        for (int i = 0; i < roofs.size(); i++) {

            string path = folder + "/roof" + to_string(i + 1);

            filesystem::create_directories(path);

            for (int j = 0; j < roofs[i].fsubroof.size(); j++) {

                string fileName = path + "/subroof" + to_string(j + 1) + ".txt";
                size_t n = printCloud(roofs[i].fsubroof[j]->view, fileName);
                cout << "roof " << i << "subroof " << j << " ----- " << n << endl;

                /*PointViewPtr inputView = *roofs[i].fsubroof[j]->view.begin();
                std::cerr
                    << "has X=" << inputView->hasDim(pdal::Dimension::Id::X)
                    << " Y=" << inputView->hasDim(pdal::Dimension::Id::Y)
                    << " Z=" << inputView->hasDim(pdal::Dimension::Id::Z)
                    << "\n";
                double x = inputView->getFieldAs<double>(Dimension::Id::X, 0);
                double y = inputView->getFieldAs<double>(Dimension::Id::Y, 0);
                double z = inputView->getFieldAs<double>(Dimension::Id::Z, 0);
                cout << "x " << x << " " << y << " " << z << endl;


                cout << "1 " << endl;
                BufferReader reader;
                reader.addView(inputView);


                string name = path + "/roof" + to_string(j + 1);
                writeToTextFile(reader, name);*/

            }

        }
    }



    size_t printCloud(const PointViewSet& pvs, const string& filePath) {

        std::ofstream out(filePath);
        //if (!out.is_open())
          //  throw std::runtime_error("Failed to open output file: " + filePath);

 
        out.setf(std::ios::fixed);
        out << std::setprecision(2);

 
        size_t written = 0;

        for (const PointViewPtr& view : pvs)
        {
            if (!view) continue;

            for (PointId i = 0; i < view->size(); ++i)
            {
                // PDAL standard coordinate dimensions:
                const double x = view->getFieldAs<double>(pdal::Dimension::Id::X, i);
                const double y = view->getFieldAs<double>(pdal::Dimension::Id::Y, i);
                const double z = view->getFieldAs<double>(pdal::Dimension::Id::Z, i);
                out << x << "," << y << "," << z << "\n";
                ++written;
            }
        }

        return written;
    }

    BestFitPlane countPlane(vector<Eigen::Vector3d> points) {
		BestFitPlane plane;


            //centroid
        double sumx = 0;
        double sumy = 0; 
        double sumz = 0;;
        for (int i = 0; i < points.size(); i++) {
            sumx += points[i].x();
            sumy += points[i].y();
            sumz += points[i].z();
        }


		plane.centroid = Eigen::Vector3d(sumx / points.size(), sumy / points.size(), sumz / points.size());

		    //centered points
		vector<Eigen::Vector3d> centeredPoints;
        for (int i = 0; i < points.size(); i++) {
			centeredPoints.push_back(points[i] - plane.centroid);
        }



		    //covariance matrix
        Eigen::Matrix3d cov = Eigen::Matrix3d::Zero();
        const int N = static_cast<int>(centeredPoints.size());
        for (const auto& p : centeredPoints) {
            cov.noalias() += p * p.transpose(); // outer product
        }


            //normal
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(cov);
        Eigen::Vector3d normal = solver.eigenvectors().col(0);
		
        plane.normal = normal.normalized();

       

            //d
		plane.d = -plane.normal.dot(plane.centroid);




        



		return plane;
    }
    void computatePlanes() {



        for (int i = 0; i < roofs.size(); i++) {
            for (int j = 0; j < roofs[i].fsubroof.size(); j++) {

                PointViewPtr view = *roofs[i].fsubroof[j]->view.begin();

                vector<Vector3d> points;
                for (PointId id = 0; id < view->size(); ++id) {
                    double x = view->getFieldAs<double>(Dimension::Id::X, id);
                    double y = view->getFieldAs<double>(Dimension::Id::Y, id);
                    double z = view->getFieldAs<double>(Dimension::Id::Z, id);
                    Vector3d point(x, y, z);
                    points.push_back(point);
                }
                BestFitPlane plane = countPlane(points);
                roofs[i].plane.push_back(plane);



            }

        }




    }
    void printPlanes(string folder) {
        error_code ec;
        filesystem::path dir = folder;
        ec.clear();

        cout << "EJ " << endl;
        filesystem::create_directories(folder);
        auto n = filesystem::remove_all(dir, ec); // deletes directory + all contents
        if (ec) {
            std::cerr << "remove_all failed: " << ec.message() << "\n";
            return;
        }
        filesystem::create_directories(folder);
        cout << "HAJ " << endl;



        for (int i = 0; i < roofs.size(); i++) {
            string path = folder + "/roof" + to_string(i + 1);
            cout << "LOOL " << endl;
            filesystem::create_directories(path);

            for (int j = 0; j < roofs[i].fsubroof.size(); j++) {

                string fileName = path + "/plane" + to_string(j + 1) + ".txt";
                
                std::ofstream out(fileName);
                //if (!out.is_open())
                  //  throw std::runtime_error("Failed to open output file: " + filePath);


                out.setf(std::ios::fixed);
                out << std::setprecision(7);

                const double x = roofs[i].plane[j].normal.x();
                const double y = roofs[i].plane[j].normal.y();
                const double z = roofs[i].plane[j].normal.z();
                const double d = roofs[i].plane[j].d;
                out << x << "," << y << "," << z << "," << d << "\n";


            }

        }

    }

    double distanceBetweenPoints(double x1, double y1, double x2, double y2) {
        double dx = x2 - x1;
        double dy = y2 - y1;
        return sqrt(dx * dx + dy * dy);
    }


    MeshReadyView buildDelaunayMesh(MeshReadyView data)
    {

        BufferReader reader;
        reader.addView(data.view);

        DelaunayFilter filter;
        filter.setInput(reader);

        filter.prepare(*data.table);
        PointViewSet outViews = filter.execute(*data.table);

        if (outViews.empty())
            throw std::runtime_error("Delaunay returned no output view.");

        data.view = *outViews.begin();
        return data;
    }
    void writeMeshToPly(const MeshReadyView& data, const string& filename)
    {

        BufferReader reader;
        reader.addView(data.view);

        PlyWriter writer;
        writer.setInput(reader);

        Options wopts;
        wopts.add("filename", filename);
        wopts.add("faces", true);
        wopts.add("storage_mode", "little endian");
        writer.setOptions(wopts);



        cout << "wr " << endl;
        writer.prepare(*data.table);
        writer.execute(*data.table);

        //PointTableRef table = view->table();

        //BufferReader reader;
        //reader.addView(view);

        //PlyWriter writer;
        //writer.setInput(reader);

        //Options wopts;
        //wopts.add("filename", filename);
        //wopts.add("faces", true);                 // write triangle faces
        //wopts.add("storage_mode", "little endian"); // binary PLY (keeps precision)
        //writer.setOptions(wopts);

        //cout << "-1 "  << endl;
        //writer.prepare(table);

        //writer.execute(table);
        //cout << "-3 " << endl;  
    
    }



    TriangularMesh cutMesh(PointViewPtr view, TriangularMesh* mesh, double alpha) {

      

        TriangularMesh outMesh;
        cout << "mesh->size() " << mesh->size() << endl;
        for (const Triangle& t : *mesh)
        {
            PointId a = t.m_a, b = t.m_b, c = t.m_c;


            double ax = view->getFieldAs<double>(Dimension::Id::X, a);
            double ay = view->getFieldAs<double>(Dimension::Id::Y, a);
            double az = view->getFieldAs<double>(Dimension::Id::Z, a);

            double bx = view->getFieldAs<double>(Dimension::Id::X, b);
            double by = view->getFieldAs<double>(Dimension::Id::Y, b);
            double bz = view->getFieldAs<double>(Dimension::Id::Z, b);

            double cx = view->getFieldAs<double>(Dimension::Id::X, c);
            double cy = view->getFieldAs<double>(Dimension::Id::Y, c);
            double cz = view->getFieldAs<double>(Dimension::Id::Z, c);


            double da = distanceBetweenPoints(ax, ay, bx, by);
            double db = distanceBetweenPoints(bx, by, cx, cy);
            double dc = distanceBetweenPoints(cx, cy, ax, ay);


            double s = (da + db + dc) / 2.0;
            double area = sqrt(s * (s - da) * (s - db) * (s - dc));
			double R = (da * db * dc) / (4.0 * area); // circumradius


            // Example policy: keep triangle only if ALL 3 vertices pass
            if (R < alpha)
            {
                outMesh.add(a, b, c);
            }

            
            
        }
        return outMesh;

    }

    void keep_only_values_with_freq_at_most_4(std::vector<int>& v) {

		int limit = 3;
        std::unordered_map<int, int> cnt;
        for (int x : v) cnt[x]++;

        v.erase(std::remove_if(v.begin(), v.end(),
            [&](int x) { return cnt[x] > limit; }),
            v.end());
    }
    void remove_repeats_sorted(std::vector<int>& v) {
        std::sort(v.begin(), v.end());
        v.erase(std::unique(v.begin(), v.end()), v.end());
    }



    vector<int> hollowMesh(PointViewPtr view, TriangularMesh mesh) {

        vector<mEdge> allEdges;
        for (const Triangle& t : mesh) {

            PointId a = t.m_a, b = t.m_b, c = t.m_c;

            mEdge e1(a < b ? a : b, a < b ? b : a);
			mEdge e2(b < c ? b : c, b < c ? c : b);
			mEdge e3(a < c ? a : c, a < c ? c : a);
			allEdges.push_back(e1);
			allEdges.push_back(e2);
			allEdges.push_back(e3);
        }
        
        map<mEdge, int> freq;

        for (const mEdge& e : allEdges) {
            freq[e]++;
        }


		vector<mEdge> borderEdges;
        for (const mEdge& e : allEdges ) {
            if (freq[e] == 1) {
				borderEdges.push_back(e);
            }
        }

        vector<int> vertices;
        for (const mEdge& e : borderEdges) {
            vertices.push_back(e.id1);
            vertices.push_back(e.id2);
        }

        vector<int> uniqueVertices;
        for (int x : vertices) {
            if (find(uniqueVertices.begin(), uniqueVertices.end(), x) == uniqueVertices.end()) {
                uniqueVertices.push_back(x);
            }
        }

 
        return uniqueVertices;
    }
    vector<Vector3d> createJacket(PointViewPtr inputView, vector<int> vertices) {

       /* pdal::PointTable table;
        table.layout()->registerDim(pdal::Dimension::Id::X);
        table.layout()->registerDim(pdal::Dimension::Id::Y);
        table.layout()->registerDim(pdal::Dimension::Id::Z);
        pdal::PointViewPtr view(new pdal::PointView(table));*/

		vector<Vector3d> points;
        for (int i = 0; i < vertices.size(); i++) {
			double x = inputView->getFieldAs<double>(pdal::Dimension::Id::X, vertices[i]);
			double y = inputView->getFieldAs<double>(pdal::Dimension::Id::Y, vertices[i]);
			double z = inputView->getFieldAs<double>(pdal::Dimension::Id::Z, vertices[i]);
			Vector3d point(x, y, z);
			points.push_back(point);

            /*view->setField(pdal::Dimension::Id::X, i, x);
            view->setField(pdal::Dimension::Id::Y, i, y);
            view->setField(pdal::Dimension::Id::Z, i, z);*/
        }
		return points;

       // pdal::PointViewSet views;
        //views.insert(view);



        /*std::ofstream out("hranicka.txt");
        out << std::fixed << std::setprecision(2);

        for (size_t i = 0; i < view->size(); ++i) {

            double ax = view->getFieldAs<double>(Dimension::Id::X, i);
            double ay = view->getFieldAs<double>(Dimension::Id::Y, i);
            double az = view->getFieldAs<double>(Dimension::Id::Z, i);
            out << ax << "," << ay << "," << az << "\n";
        }
        out << '\n';*/

    }




    vector<Vector3d> clipBorder(vector<Vector3d> points) {
        double upper_limit = 0.20 * points.size();
        double lower_limit = 0.05 * points.size();

        vector<Vector3d> points0 = points;
        vector<vector<Vector3d>> edge_points;
        vector<Vector3d> jacket;
  

        bool stop = false;
		vector<Vector3d> rest_points;

		int count = 1;
		string filename = "rest" + to_string(count) + ".txt";
        while (points0.size() > 0) {
            
            
            vector<Vector3d> localBorder = findEdgePoints(points0, 0.5, lower_limit, upper_limit, stop);
  

            if(stop){
				rest_points = localBorder;
                break;
			}
            edge_points.push_back(localBorder);

            
            points0.erase(
                std::remove_if(points0.begin(), points0.end(),
                    [&](const Vector3d& p)
                    {
                        return std::find(localBorder.begin(), localBorder.end(), p) != localBorder.end();
                    }),
                points0.end());

            printPoints(points0, filename);
            count++;
            filename = "rest" + to_string(count) + ".txt";

            //printPoints(rest_points, "rest.txt");

        }
		
 

		stop = false;
        upper_limit = 0.15 * rest_points.size();
        lower_limit = 0.05 * rest_points.size();
        while (rest_points.size() > 0) {
            vector<Vector3d> localBorder = findEdgePoints(rest_points, 0.5, lower_limit, upper_limit, stop);


            if (stop) {
                rest_points = localBorder;
                break;
            }
            edge_points.push_back(localBorder);
            
            rest_points.erase(
                std::remove_if(rest_points.begin(), rest_points.end(),
                    [&](const Vector3d& p)
                    {
                        return std::find(localBorder.begin(), localBorder.end(), p) != localBorder.end();
                    }),
                rest_points.end());

        }
  




        for (int i = 0; i < edge_points.size(); i++) {
    
			Line3D line;
            fitLine3D(edge_points[i], line);

			string filename = "lines/line" + to_string(i + 1) + ".txt";
            
            writeLineToFile(filename, line.p, line.dir);
        }

        for (int i = 0; i < edge_points.size(); i++) {
            Vector3d p1;
            Vector3d p2;
            findCornerPoints(edge_points[i], p1, p2);

            jacket.push_back(p1);
            jacket.push_back(p2);

        }
        for (int i = 0; i < rest_points.size(); i++) {
            jacket.push_back(rest_points[i]);

        }

        return jacket;


      

        

    }

    void findCornerPoints(vector<Vector3d> points, Vector3d &p1, Vector3d& p2) {

        double bestDist2 = -1.0;
        //pair<Vector3d, Vector3d> bestPair;

        for (size_t i = 0; i < points.size(); i++) {
            for (size_t j = i + 1; j < points.size(); j++) {
                double dist2 = (points[i] - points[j]).squaredNorm();

                if (dist2 > bestDist2) {
                    bestDist2 = dist2;
                    p1 = points[i];
                    p2 = points[j];   
                }
            }
        }

    }
    vector<Vector3d> uniteClosePoints(vector<Vector3d> points, double ratio) {
        vector<Vector3d> result;
        vector<bool> used(points.size(), false);

        for (size_t i = 0; i < points.size(); i++) {
            if (used[i]) continue;

            Vector3d sum = points[i];
            int count = 1;
            used[i] = true;

            for (size_t j = i + 1; j < points.size(); j++) {
                if (!used[j] && (points[i] - points[j]).norm() < ratio) {
                    sum += points[j];
                    count++;
                    used[j] = true;
                }
            }

            result.push_back(sum / count);
        }

        /*vector<Vector3d> united_points;
        for (int i = 0; i < points.size(); i++) {

            bool has_pair = false;
            for (int j = i + 1; j < points.size(); j++) {

                Vector3d p1 = points[i];
                Vector3d p2 = points[j];
                double d = (p1 - p2).norm();
                cout << "d " << d << endl;

                if (d < ratio) {
                    cout << "truuuuuuu " << endl;
                    Vector3d unity = (p1 + p2) / 2.0;
                    united_points.push_back(unity);
                    has_pair = true;
                }
            }

            if (!has_pair) {
                cout << "whyyy " << endl;
                united_points.push_back(points[i]);
            }

            
        }*/

        return result;

    }
    vector<Vector3d> orderPoints(vector<Vector3d> jacket, Vector3d centroid) {



        Vector3d anchor_point = jacket[0];
        vector<double> angles;
        angles.push_back(0.);
        double b = (centroid - anchor_point).norm();
        Vector3d AB = anchor_point - centroid;
        for (int i = 1; i < jacket.size(); i++)
        {
            double a = (centroid - jacket[i]).norm();
            double c = (jacket[i] - anchor_point).norm();

            double gama;
            Vector3d AP = jacket[i] - centroid;
            if ((AB.x() * AP.y() - AB.y() * AP.x()) > 0) {
                gama = acos((a * a + b * b - c * c) / (2 * a * b));
                angles.push_back(gama);
            }
            else {

                gama = acos((a * a + b * b - c * c) / (2 * a * b));
                gama = (2 * PI) - gama;
                angles.push_back(gama);
            }

        }

        vector<size_t> indices(jacket.size());
        iota(indices.begin(), indices.end(), 0);

        sort(indices.begin(), indices.end(),
            [&](size_t a, size_t b) {
                return angles[a] < angles[b];
            });

        vector<Vector3d> sortedPoints;
        sortedPoints.reserve(jacket.size());

        for (size_t i : indices) {
            sortedPoints.push_back(jacket[i]);
        }


        return sortedPoints;



        /*int n = jacket.size();
        vector<bool> used(n, false);
        vector<Vector3d> ordered;
        ordered.reserve(n);


        int startIndex = 0;
        int current = startIndex;
        ordered.push_back(jacket[current]);
        used[current] = true;

        for (int step = 1; step < n; step++)
        {
            int bestIndex = -1;
            double bestDist = (std::numeric_limits<double>::max)();

            for (int i = 0; i < n; i++)
            {
                if (used[i])
                    continue;

                double d = (jacket[i] - jacket[current]).squaredNorm();
                if (d < bestDist)
                {
                    bestDist = d;
                    bestIndex = i;
                }
            }

            current = bestIndex;
            used[current] = true;
            ordered.push_back(jacket[current]);
        }

        return ordered;*/


    }
    void calculateCentroid() {

        auto start = chrono::high_resolution_clock::now();
        cout << "BEGINNING centroid CALCULATING " << endl;

        for (int i = 0; i < roofs.size(); i++)
        {

            for (int j = 0; j < roofs[i].fsubroof.size(); j++)
            {
                
                double sumx = 0;
                double sumy = 0;
                double sumz = 0;
                PointViewPtr view = *roofs[i].fsubroof[j]->view.begin();
                for (int k = 0; k < view->size(); k++)
                {
                    

                    const double x = view->getFieldAs<double>(Dimension::Id::X, k);
                    const double y = view->getFieldAs<double>(Dimension::Id::Y, k);
                    const double z = view->getFieldAs<double>(Dimension::Id::Z, k);
                    sumx += x;
                    sumy += y;
                    sumz += z;
                }
                roofs[i].fsubroof[j]->centroid.x() = sumx / view->size();
                roofs[i].fsubroof[j]->centroid.y() = sumy / view->size();
                roofs[i].fsubroof[j]->centroid.z() = sumz / view->size();

                
            }


        }






        chrono::duration<double> elapsedd = chrono::high_resolution_clock::now() - start;
        cout << "\n\nFINISHING " << elapsedd.count() << " seconds.\n";
    }
    void printCentroids(string foldername) {

        




        error_code ec;
        filesystem::path dir = foldername;
        auto n = filesystem::remove_all(dir, ec); // deletes directory + all contents
        filesystem::create_directories(foldername);


        for (int i = 0; i < roofs.size(); i++)
        {
            string path = foldername + "/roof" + to_string(i + 1);

            filesystem::create_directories(path);
            for (int j = 0; j < roofs[i].fsubroof.size(); j++)
            {
                Vector3d centroid = roofs[i].fsubroof[j]->centroid;

                string fileName = path + "/centroid" + to_string(j + 1) + ".txt";
                ofstream out(fileName);
                out << fixed << std::setprecision(3);

                out << centroid.x() << "," << centroid.y() << "," << centroid.z() << endl;
            }


        }

    }


    void writeLineToFile(const std::string& filename, const Vector3d& p0, const Vector3d& dir)
    {
        std::ofstream file(filename, std::ios::app);
  
  

        file << std::fixed << std::setprecision(6) << p0.x() << ", " << p0.y() << ", " << p0.z() << "\n";
        file << std::fixed << std::setprecision(6) << dir.x() << ", " << dir.y() << ", " << dir.z() << "\n";
    }

    vector<Vector3d> findEdgePoints(vector<Vector3d> all_points, double ratio, double points_lower_limit, double points_upper_limit, bool& stop) {
        cout << "finding Edge Points  " << endl;
        cout << "all_points " << all_points.size()   << endl;

        for (int i = 0; i < all_points.size(); i++) {
            for (int j = i + 1; j < all_points.size(); j++) {



                //first round
				size_t points_count = 0;
				vector<Vector3d> line_points;
                Line3D line;
                line.p = all_points[i];
                line.dir = all_points[j] - all_points[i];
                for (int k = 0; k < all_points.size(); k++) {
                    if(k!= i && k != j){

                        Vector3d C = all_points[k];
                        double d = distanceToLine(line, C);
                        if (d < ratio) {
							line_points.push_back(C);
                            points_count++;
                        }
					}
                }
                line_points.push_back(all_points[i]);
                line_points.push_back(all_points[j]);
                points_count++;
                points_count++;


                

                ////second round - better line fit
                if (points_count > points_lower_limit) {
                                   
                    points_count = 0;
                    fitLine3D(line_points, line);
                    line_points.clear();
                    for (int k = 0; k < all_points.size(); k++) {
                        if (k != i && k != j) {

         
                            Vector3d C = all_points[k];
                            double d = distanceToLine(line, C);
                            if (d < ratio) {
                                line_points.push_back(C);
                                points_count++;

                            }
                        }
                    }
                    if (points_count > points_upper_limit) {
                        cout << "--------points_count   " << points_count << endl;
						line_points.push_back(all_points[i]);
						line_points.push_back(all_points[j]);
                        return line_points;
                    }
                }



            }
        }

		stop = true;
		return all_points;
    }
    void fitLine3D(const vector<Vector3d>& points, Line3D& line) {
        int n = static_cast<int>(points.size());
        if (n == 0)
            throw std::runtime_error("No points.");

        // 1. centroid
        line.p.setZero();
        for (const auto& p : points)
            line.p += p;
        line.p /= static_cast<double>(n);

        // 2. covariance-like matrix
        Matrix3d S = Matrix3d::Zero();
        for (const auto& p : points)
        {
            Vector3d d = p - line.p;
            S += d * d.transpose();
        }

        // 3. dominant eigenvector by power iteration
        Vector3d v(1.0, 1.0, 1.0);
        v.normalize();

        for (int i = 0; i < 100; ++i)
        {
            v = S * v;              // matrix * vector
            if (v.norm() == 0.0)    // safety check
                break;
            v.normalize();          // normalize in place
        }

        line.dir = v;
    }





    void attachTriangularMesh(PointViewPtr view, TriangularMesh mesh)
    {
        TriangularMesh* dst = view->createMesh("");
        *dst = std::move(mesh);
        cout << "reduced size " << dst->size() << endl;
        //return dst;
    }
    void getBorder(PointViewPtr view, vector<int> *input, string filename) {
 
        std::ofstream out(filename);
        out << std::fixed << std::setprecision(2);

        for (size_t i = 0; i < input->size(); ++i) {

            double ax = view->getFieldAs<double>(Dimension::Id::X, input->at(i));
            double ay = view->getFieldAs<double>(Dimension::Id::Y, input->at(i));
            double az = view->getFieldAs<double>(Dimension::Id::Z, input->at(i));

			out << ax << "," << ay << "," << az << "\n";
        }
        out << '\n';
	

    }
    void printPoints(vector<Vector3d> input, string filename) {


        std::ofstream out(filename);
        out << std::fixed << std::setprecision(2);

        for (size_t i = 0; i < input.size(); ++i) {

            double ax = input[i].x();
            double ay = input[i].y();
            double az = input[i].z();

            out << ax << "," << ay << "," << az << "\n";
        }
        out << '\n';
    }



    PointViewPtr createFilteredMeshCompact(
        pdal::PointViewPtr srcView,
		double alpha /*= 0.5*/,
        const std::string& srcMeshName /*= "delaunay2d"*/,
        const std::string& dstMeshName /*= "filtered"*/)
    {

        TriangularMesh* srcMesh = srcView->mesh(srcMeshName);
        PointViewPtr outView = srcView->makeNew();

        // New mesh attached to outView
        TriangularMesh* outMesh = outView->createMesh(dstMeshName);
  

        // Map old point IDs -> new point IDs (for remapping triangle indices)
        unordered_map<PointId, PointId> idMap;
        idMap.reserve(srcView->size());

        auto getOrCopyPoint = [&](PointId oldId) -> PointId
            {
                auto it = idMap.find(oldId);
                if (it != idMap.end())
                    return it->second;

                // Copy one point from srcView to outView (same dimensions)
                pdal::PointId newId = outView->size();
                outView->appendPoint(*srcView, oldId);   // <-- PDAL API used for point copy
                idMap.emplace(oldId, newId);
                return newId;
            };

        for (const Triangle& t : *srcMesh)
        {
            PointId a = t.m_a, b = t.m_b, c = t.m_c;

            double ax = srcView->getFieldAs<double>(Dimension::Id::X, a);
            double ay = srcView->getFieldAs<double>(Dimension::Id::Y, a);
            double az = srcView->getFieldAs<double>(Dimension::Id::Z, a);

            double bx = srcView->getFieldAs<double>(Dimension::Id::X, b);
            double by = srcView->getFieldAs<double>(Dimension::Id::Y, b);
            double bz = srcView->getFieldAs<double>(Dimension::Id::Z, b);

            double cx = srcView->getFieldAs<double>(Dimension::Id::X, c);
            double cy = srcView->getFieldAs<double>(Dimension::Id::Y, c);
            double cz = srcView->getFieldAs<double>(Dimension::Id::Z, c);


            double da = distanceBetweenPoints(ax, ay, bx, by);
            double db = distanceBetweenPoints(bx, by, cx, cy);
            double dc = distanceBetweenPoints(cx, cy, ax, ay);

			double s = (da + db + dc) / 2.0;
			double area = sqrt(s * (s - da) * (s - db) * (s - dc));
            // Example policy: keep triangle only if ALL 3 vertices pass
            if (area < alpha)
            {
                continue;
            }

            // Copy points lazily and remap indices
            PointId na = getOrCopyPoint(a);
            PointId nb = getOrCopyPoint(b);
            PointId nc = getOrCopyPoint(c);

            outMesh->add(na, nb, nc);
        }

        return outView;
    }


    bool intersectPlanes(const BestFitPlane& p1, const BestFitPlane& p2, Line3D& out, double eps = 1e-12)
    {
        Eigen::Vector3d u = p1.normal.cross(p2.normal);
        double denom = u.squaredNorm();

        if (denom < eps) {
            // Planes are parallel (or nearly). Could be distinct or the same plane.
            return false;
        }

        // Point on intersection line:
        Eigen::Vector3d p0 = ((p2.d * p1.normal - p1.d * p2.normal).cross(u)) / denom;

        out.p = p0;
        out.dir = u; // you can normalize if you want: u.normalized()
        return true;
    }
    void computateLines() {

        for (int i = 0; i < roofs.size(); i++) {
            cout << "COMPUTING LINES OF ROOF....." << i << endl;


			int linescount = 0;
            for (int j = 0; j < roofs[i].plane.size() - 1; j++) {
                for (int k = j + 1; k < roofs[i].plane.size(); k++) {

					linescount++;
                    Line3D line;
                    intersectPlanes(roofs[i].plane[j], roofs[i].plane[k], line);
                    Vector2d origin(j, k);
                    line.origin_planes = origin;
                    //cout << "Line " << linescount << " of plane " << j + 1 << " and plane " << k + 1 << endl;
                    //cout << line.p.x() << "," << line.p.y() << "," << line.p.z() << endl;
					//cout << line.dir.x() << "," << line.dir.y() << "," << line.dir.z() << endl << endl;


                    if (neighbors(roofs[i].jacket[j], roofs[i].jacket[k])) {
                        roofs[i].lines.push_back(line);
                        cout << "valuable line found "  << endl;
                    }             
                }

                
            }
        }
	}
    bool neighbors(vector<Vector3d> one, vector<Vector3d> two) {
        double ratio = 2.5;



        for (int i = 0; i < one.size(); i++)
        {
            for (int j = 0; j < two.size(); j++)
            {


                double distance = (one[i] - two[j]).norm();
                cout << "distance " << distance << endl;
                if (distance < ratio) {
                    return true;
                }
            }
        }
        return false;



    }


    void computateMeshes() {


        for (int i = 0; i < roofs.size(); i++) {

            for (int j = 0; j < roofs[i].fsubroof.size(); j++) {

                auto clean = makeMeshReadyViewDeep(*roofs[i].fsubroof[j]->view.begin());
                clean = buildDelaunayMesh(clean);
                //TriangularMesh* mesh = clean.view->mesh("delaunay2d");
                
                roofs[i].mesh_ready_view.push_back(clean);

				TriangularMesh* mesh = roofs[i].mesh_ready_view[j].view->mesh("");
                
                roofs[i].initial_mesh.push_back(mesh);
            }
        }

    }
    MeshReadyView makeMeshReadyView(const pdal::PointViewPtr& src)
    {
        using namespace pdal;

        auto table = std::make_shared<PointTable>();
        auto layout = table->layout();

        auto keepIfExists = [&](Dimension::Id id)
            {
                if (src->layout()->hasDim(id))
                    layout->registerDim(id);
            };

        // Keep only dimensions safe/useful for mesh export
        keepIfExists(Dimension::Id::X);
        keepIfExists(Dimension::Id::Y);
        keepIfExists(Dimension::Id::Z);
        keepIfExists(Dimension::Id::NormalX);
        keepIfExists(Dimension::Id::NormalY);
        keepIfExists(Dimension::Id::NormalZ);
        keepIfExists(Dimension::Id::Curvature);

        pdal::PointViewPtr dst(new pdal::PointView(*table));

        for (pdal::PointId i = 0; i < src->size(); ++i)
            dst->appendPoint(*src, i);   // copies only dimensions registered above

        return { table, dst };
    }
    MeshReadyView makeMeshReadyViewDeep(const pdal::PointViewPtr& src)
    {
        using namespace pdal;

        MeshReadyView out;
        out.table = std::make_shared<PointTable>();
        auto layout = out.table->layout();

        auto keep = [&](Dimension::Id id)
            {
                if (src->layout()->hasDim(id))
                    layout->registerDim(id);
            };

        keep(Dimension::Id::X);
        keep(Dimension::Id::Y);
        keep(Dimension::Id::Z);
        keep(Dimension::Id::NormalX);
        keep(Dimension::Id::NormalY);
        keep(Dimension::Id::NormalZ);
        keep(Dimension::Id::Curvature);

        out.view.reset(new PointView(*out.table));

        for (pdal::PointId i = 0; i < src->size(); ++i)
        {
            pdal::PointId id = out.view->size();

            if (src->layout()->hasDim(pdal::Dimension::Id::X))
                out.view->setField(pdal::Dimension::Id::X, id,
                    src->getFieldAs<double>(pdal::Dimension::Id::X, i));

            if (src->layout()->hasDim(pdal::Dimension::Id::Y))
                out.view->setField(pdal::Dimension::Id::Y, id,
                    src->getFieldAs<double>(pdal::Dimension::Id::Y, i));

            if (src->layout()->hasDim(pdal::Dimension::Id::Z))
                out.view->setField(pdal::Dimension::Id::Z, id,
                    src->getFieldAs<double>(pdal::Dimension::Id::Z, i));

            if (src->layout()->hasDim(pdal::Dimension::Id::NormalX))
                out.view->setField(pdal::Dimension::Id::NormalX, id,
                    src->getFieldAs<double>(pdal::Dimension::Id::NormalX, i));

            if (src->layout()->hasDim(pdal::Dimension::Id::NormalY))
                out.view->setField(pdal::Dimension::Id::NormalY, id,
                    src->getFieldAs<double>(pdal::Dimension::Id::NormalY, i));

            if (src->layout()->hasDim(pdal::Dimension::Id::NormalZ))
                out.view->setField(pdal::Dimension::Id::NormalZ, id,
                    src->getFieldAs<double>(pdal::Dimension::Id::NormalZ, i));

            if (src->layout()->hasDim(pdal::Dimension::Id::Curvature))
                out.view->setField(pdal::Dimension::Id::Curvature, id,
                    src->getFieldAs<double>(pdal::Dimension::Id::Curvature, i));
        }

        return out;
    }

    void printMeshes(string folder) {

        filesystem::create_directories(folder);
        for (int i = 0; i < roofs.size(); i++) {


            string path = folder + "/roof" + to_string(i + 1);
            filesystem::create_directories(path);
            for (int j = 0; j < roofs[i].fsubroof.size(); j++) {

				auto clean = roofs[i].mesh_ready_view[j];
               // auto clean = makeMeshReadyViewDeep(*roofs[i].fsubroof[j]->view.begin());
                //clean = buildDelaunayMesh(clean);

                std::cout << "points: " << clean.view->size() << '\n';
                std::cout << "mesh: " << (clean.view->mesh("delaunay2d") ? "yes" : "no") << '\n';

                cout << "j " << j << endl;
                string name = folder + "/roof" + to_string(i + 1) + "/subroof" + to_string(j + 1) + "mesh.ply";
                cout << "name     " << name << endl;

                writeMeshToPly(clean, name);
     
            }
        }
    }
    void printLines(string folder) {


        filesystem::create_directories(folder);
        for (int i = 0; i < roofs.size(); i++) {
            string path = folder + "/roof" + to_string(i + 1);
            filesystem::create_directories(path);

            for (int j = 0; j < roofs[i].lines.size(); j++) {
                string fileName = path + "/line" + to_string(j + 1) + ".txt";
                Line3D line = roofs[i].lines[j];
                
                ofstream out(fileName);
                out << std::fixed << std::setprecision(8);


                out << line.p.x() << "," << line.p.y() << "," << line.p.z() << "\n";
                out << line.dir.x() << "," << line.dir.y() << "," << line.dir.z() << "\n";
                out << '\n';
            }
        }
    }

    void cutRoofs(double alpha) {


        cout << "roofs.size() " << roofs.size() << endl;

        for (int i = 0; i < roofs.size(); i++) {

            cout << "roofs[i].fsubroof.size() " << roofs[i].fsubroof.size() << endl;
            for (int j = 0; j < roofs[i].fsubroof.size(); j++) {

                //TriangularMesh nm = cutMesh(*roofs[i].fsubroof[j]->view.begin(), roofs[i].initial_mesh[j], 0.25);
                TriangularMesh nm = cutMesh(roofs[i].mesh_ready_view[j].view, roofs[i].initial_mesh[j], alpha);
    
                //attachTriangularMesh(*roofs[i].fsubroof[j]->view.begin(), nm);
                attachTriangularMesh(roofs[i].mesh_ready_view[j].view, nm);

				roofs[i].new_mesh.push_back(nm);

            }
        }
    }
    void getJackets() {

        for (int i = 0; i < roofs.size(); i++) {
            for (int j = 0; j < roofs[i].fsubroof.size(); j++) {

                TriangularMesh m = roofs[i].new_mesh[j];
                vector<int> jacketIDs = hollowMesh(*roofs[i].fsubroof[j]->view.begin(), m);
                roofs[i].jacket.push_back(createJacket(*roofs[i].fsubroof[j]->view.begin(), jacketIDs));
                
                cout << "JACKET SIZE    " << roofs[i].jacket[j].size() << endl;
            }
        }
    }
    void printJackets(string folder) {

        filesystem::create_directories(folder);
        for (int i = 0; i < roofs.size(); i++) {
            string path = folder + "/roof" + to_string(i + 1);
            filesystem::create_directories(path);

            for (int j = 0; j < roofs[i].fsubroof.size(); j++) {
                string fileName = path + "/jacket" + to_string(j + 1) + ".txt";
                printPoints(roofs[i].jacket[j], fileName);
            }
        }
    }


    void modifyJackets() {
        
        
        for (int i = 0; i < roofs.size(); i++) {
            for (int j = 0; j < roofs[i].fsubroof.size(); j++) {
                
                
                vector<Vector3d> jacket = roofs[i].jacket[j];
                Vector3d centroid = roofs[i].fsubroof[j]->centroid;


                jacket = clipBorder(jacket);
                jacket = uniteClosePoints(jacket, 2.);
                jacket = orderPoints(jacket, centroid);

                roofs[i].jacket[j] = jacket;
            }
        }



        
    }
    void uniteAllPoints(double ratio) {


        
        for (int i = 0; i < roofs.size(); i++)
        {
            
            //collecting all main roof points
            vector<Vector3d> mainroof_points;
            for (int j = 0; j < roofs[i].jacket.size(); j++)
            {
                for (int k = 0; k < roofs[i].jacket[j].size(); k++)
                {

                    Vector3d point = roofs[i].jacket[j][k];
                    mainroof_points.push_back(point);
                }     
                
            }


            //uniting similar points
            for (int j = 0; j < mainroof_points.size(); j++)
            {

                vector<Vector3d> united;
                vector<int> ids;
                bool found = false;
                for (int k = j + 1; k < mainroof_points.size(); k++)
                {
                    
                    double distance = (mainroof_points[j] - mainroof_points[k]).norm();
                    if (distance < ratio) {
                        united.push_back(mainroof_points[k]);
                        ids.push_back(k);
                        found = true;
                    }
                }
                if (found) {
                    united.push_back(mainroof_points[j]);
                    ids.push_back(j);
                }


                //rewriting the points
                double sumx = 0;
                double sumy = 0;
                double sumz = 0;
                for (int k = 0; k < united.size(); k++)
                {
                    sumx += united[k].x();
                    sumy += united[k].y();
                    sumz += united[k].z();
                }
                Vector3d united_point(sumx / united.size(), sumy / united.size(), sumz / united.size());
                for (int k = 0; k < united.size(); k++)
                {
                    mainroof_points[ids[k]] = united_point;
                }

            }



            //regrouping
            int count = 0;
            for (int j = 0; j < roofs[i].jacket.size(); j++)
            {
                for (int k = 0; k < roofs[i].jacket[j].size(); k++)
                {
                    roofs[i].jacket[j][k] = mainroof_points[count];
                    count++;
                }

            }


        }

    }

    void shiftPoints(double ratio) {

        for (int i = 0; i < roofs.size(); i++) {

            //collecting all main roof points
            vector<Vector3d> mainroof_points;
            vector<int> subroof;
            for (int j = 0; j < roofs[i].jacket.size(); j++)
            {
                subroof.push_back(i);

                for (int k = 0; k < roofs[i].jacket[j].size(); k++)
                {

                    Vector3d point = roofs[i].jacket[j][k];
                    mainroof_points.push_back(point);
                }
            }



            //cout << "ROOF " << i << endl;
            for (int j = 0; j < mainroof_points.size(); j++) {



              // cout << "j " << j << endl;
                //cout << "roofs[i].lines.size()   " << roofs[i].lines.size() << endl;
                for (int k = 0; k < roofs[i].lines.size(); k++)
                {

                    Line3D line = roofs[i].lines[k];
                    int origin_plane1 = roofs[i].lines[k].origin_planes.x();
                    int origin_plane2 = roofs[i].lines[k].origin_planes.y();
                    if ((origin_plane1 == subroof[j]) || (origin_plane2 == subroof[j])) {

                    }


                    double denom = line.dir.dot(line.dir);
                    if (denom == 0.0) {
                        throw std::runtime_error("Line direction vector has zero length");
                    }
                    double dist = (mainroof_points[j] - line.p).cross(line.dir).norm() / denom;

                    //cout << "dist " << dist << endl;


                    if (dist < ratio) {

                        double t = line.dir.dot(mainroof_points[j] - line.p) / denom;
                        Vector3d unite = line.p + t * line.dir;
                        mainroof_points[j] = unite;

                        //cout << "ITEM FOUND "  << endl;
                    }
                    
                }
            }



            //regrouping
            int count = 0;
            for (int j = 0; j < roofs[i].jacket.size(); j++)
            {
                for (int k = 0; k < roofs[i].jacket[j].size(); k++)
                {
                    roofs[i].jacket[j][k] = mainroof_points[count];
                    count++;
                }
            }


        }


    }


};


struct LoadedView
{
    pdal::PointTable table;   // MUST outlive the view
    pdal::PointViewPtr view;
};
pdal::PointViewPtr loadXYZ_withReaderText(pdal::PointTable& table, const std::string& filename)
{
    pdal::StageFactory factory;


    Stage* reader(factory.createStage("readers.text"));
    pdal::Options opts;
    opts.add("filename", filename);
    opts.add("separator", ",");
    opts.add("header", "X,Y,Z");
    opts.add("skip", 1); // skip "x,y,z"

    reader->setOptions(opts);

    reader->prepare(table);
    pdal::PointViewSet views = reader->execute(table);

    //if (views.empty())
    //    throw std::runtime_error("No points read from file");

    return *views.begin();
}

int main() {
    _putenv_s("PROJ_LIB", "C:\\vcpkg\\installed\\x64-windows\\share\\proj");
    _putenv_s("PROJ_DATA", "C:\\vcpkg\\installed\\x64-windows\\share\\proj");
    _putenv_s("GDAL_DATA", "C:\\vcpkg\\installed\\x64-windows\\share\\gdal");  
    auto start = chrono::high_resolution_clock::now();
    PointReader p;

    
    /*pdal::PointTable table;
    pdal::PointViewPtr view = loadXYZ_withReaderText(table, "C:/C/PointReader/build/finalRoofs/roof1/subroof2.txt");
    cout << "size " << view->size() << endl;






    auto clean = p.makeMeshReadyViewDeep(view);
    clean = p.buildDelaunayMesh(clean);
    TriangularMesh* mesh = clean.view->mesh("");

    TriangularMesh nm = p.cutMesh(clean.view, mesh, 1);
    p.attachTriangularMesh(clean.view, nm);


  
    vector<int> k = p.hollowMesh(view, nm);
    p.getBorder(view, &k, "hranicka.txt");
	vector<Vector3d> jacket = p.createJacket(view, k);
	
    jacket = p.clipBorder(jacket);
    cout << "jacket size " << jacket.size() << endl;
    jacket = p.uniteClosePoints(jacket, 2.);
    jacket = p.orderPoints(jacket);

    cout << "jacket size " << jacket.size() << endl;
    p.printPoints(jacket, "finalJacket.txt");*/
    





 //   p.writeMeshToPly(view, "C:/C/PointReader/build/mesh_2.ply");
 //   TriangularMesh nm = p.cutMesh(view, m, 0.25);
	//p.attachTriangularMesh(view, nm);
 //   p.writeMeshToPly(view, "C:/C/PointReader/build/mesh2_2.ply");
	//vector<PointId> k = p.hollowMesh(view, nm);
 //   p.getBorder(view, &k, "border2.txt");




    string out_file = "output";
    error_code ec;
    filesystem::path dir = out_file;
    ec.clear();

    filesystem::remove_all(dir, ec); // deletes directory + all contents
    filesystem::create_directories(out_file);
    


    Stage* stage;
    stage = p.loadFile("LiDAR.laz");
    stage = p.filterClass6(stage, out_file);
    stage = p.computeNormals(stage, 16);
    stage = p.filterWalls(stage, out_file, 0.01, 20);
    stage = p.filterOutliers(stage, out_file, "statistical", 6, 0.5, 1, 4);
    stage = p.clusterPoints(stage, 0.25, 60);            //klasterizacia vsetkych bodov na hlavne strechy
    stage = p.zsmooth(stage, out_file, 1);
    p.execute(stage);
    //p.filerByZValue();




    p.makePointViewSetRoofs();      //   buildpoint --> roofs
    p.printMainRoofs(*p.buildpoint.begin(), out_file, "mainroofs");


    ////p.clusterByPoints("mainroofs", 0.25, 60);                //klasterizacia hlavnych striech na podstrechy
    p.clusterByPoints("mainroofs", 1, 60); 
    //p.makeClusteredRoofsFiles("subroofs");
     
    
    p.fillSubroofs();
    //p.modifySubroofs();
    //p.clusterByNormals(0.01, 40);            //klasterizacia podstriech na mensie podstrechy
    //p.finalizeSubroofs();
    //p.printFinalRoofs("finalRoofs");



    //p.clusterByPoints("nroofs", 0.5, 70);
    //p.endSubroofs();
    //p.printFinalRoofs22("finalRoofsTotal");



    //p.calculateCentroid();
    //p.printCentroids("Centroids");
    //cout << "1 " << endl;
    //p.computatePlanes();
    //cout << "2 " << endl;
    //p.printPlanes("planes");
    //cout << "3 " << endl;


    //p.computateMeshes();
    ////p.printMeshes("initialMeshes");
    //p.cutRoofs(0.25);
    ////p.printMeshes("cutMeshes25");
    //p.getJackets();
    //p.modifyJackets();

    //p.printJackets("jacketsFINAL");

    //p.computateLines();
    //p.printLines("ciary");




    //p.shiftPoints(3.2);
    //p.uniteAllPoints(1.25);
    //p.printJackets("jacketsUnitedFINAL");







    chrono::duration<double> elapsedd = chrono::high_resolution_clock::now() - start;
    cout << "\n\nProgram ran in " << elapsedd.count() << " seconds.\n";

}
