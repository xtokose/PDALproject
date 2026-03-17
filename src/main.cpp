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

const double PI = 3.14159265359;
using namespace std;
using namespace pdal;

struct BestFitPlane
{
    
    Eigen::Vector3d normal;   // unit-length (A,B,C)
    double d;                 // D in Ax + By + Cz + D = 0
    Eigen::Vector3d centroid; // point on plane
};
struct Line3D {
    Eigen::Vector3d p;   // point on line
    Eigen::Vector3d dir; // direction (not necessarily normalized)
    vector<int> origin;
};


struct MeshReadyView
{
    std::shared_ptr<pdal::PointTable> table; // must stay alive with the view
    pdal::PointViewPtr view;
};
struct RoofData {
    PointTable table;
    PointViewSet view; // contains the PointViewPtrs made from `table`
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
    vector<vector<PointId>> jacket;


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
    bool debug = false;
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
        if(debug){
            cout << "\nLoading finished in " << duration.count()/1000. << " seconds\n";
            cout << "Points remaining: " << get_points_size(all_points) << endl;
        }
        return reader;       
    }
    Stage* filterClass6(Stage* input) {
        auto start = chrono::high_resolution_clock::now();

        Stage* range = factory.createStage("filters.range");
        Options opts;
        opts.add("limits", "Classification[6:6]");
        range->setOptions(opts);
        range->setInput(*input);

        


            //debug
        auto duration = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start);
        cout << "\n-- Filtering class 6 finished in " << duration.count() / 1000. << " miliseconds\n";
        if(debug){ writeToTextFile(*range, "buildpoints/points_after_filtering_class_6.txt"); }
    

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
    Stage* filterWalls(Stage* input, double ratio_curvature, double ratio_angle) {
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
        if (debug){ writeToTextFile(*my, "buildpoints/points_after_wall_filter.txt"); }

        return my;
    }
    Stage* filterOutliers(Stage* input, string method, size_t mean_k, float multiplier, float radius, size_t min_k) {
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
        if (debug) { writeToTextFile(rangeb, "buildpoints/points_after_outliers_filter.txt"); }
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
    Stage* zsmooth(Stage* input, size_t radius) {
        auto start = chrono::high_resolution_clock::now();

        Options zsmoothOptions;
        zsmoothOptions.add("radius", radius);
        zsmoothOptions.add("dim", "Zsmooth");
        //zsmoothOptions.add("medianpercent", 0);

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
        if (debug){ writeToTextFile(*assign, "buildpoints/points_after_smoothing.txt"); }

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
        cout << "3 " << endl;

        w->setOptions(o);
        w->setInput(input);
        return w;
    }
    void writeToTextFile(Stage& input, const string& filename) {
        
        cout << "2 " << endl;
        auto writer = addTextWriter(input, filename);
        cout << "4 " << endl;

        PointTable table;
        writer->prepare(table);
        cout << "5 " << endl;


        auto* layout = table.layout();
        auto dimX = layout->findDim("X");       // returns Dimension::Id::Unknown if missing
        std::cerr << "X dim id = " << (int)dimX << "\n";

        try {
            //writer->prepare(table);
            writer->execute(table);
        }
        catch (const pdal::pdal_error& e) {
            std::cerr << "PDAL error: " << e.what() << "\n";
        }


        //writer->execute(table);
        cout << "6 " << endl;

        ifstream in(filename);
        string line;
        size_t lines = 0;
        while (std::getline(in, line)) {
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
 


    void makeClusteredFiles(const PointViewPtr& view , string folder) {
     
        //getting number of clusters
        unordered_set<double> clusterIDs = getRoofsIDs(view);
        double numClusters = clusterIDs.size();



    
            //removing existing files
        filesystem::create_directories(folder);
        filesystem::path path = folder;
        for (const auto& entry : filesystem::directory_iterator(path)) {
            if (!filesystem::is_regular_file(entry.status()))
                continue; // skip subdirs etc.

            if (entry.path().extension() == ".txt") {
                filesystem::remove(entry.path());

            }
        }




            //creating empty files
        vector<ofstream> files;
        files.reserve(numClusters);
        for (int i = 0; i < numClusters; i++) {
            files.emplace_back(folder + "/roof" + to_string(i + 1) + ".txt");
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

    void makePointViewSetRoofs() {
        PointViewPtr old_view = *buildpoint.begin();

		//int sum = 0;

        for (int k = 0; k < number_of_roofs; k++) {
            auto rd = std::make_unique<RoofData>();

            //PointTable table;
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
            cout << "Roof " << (k + 1) << " ----- " << view->size() << " points " << endl;

			//sum += view->size();
            //PointViewSet temp;
            //temp.insert(view);

            rd->view.insert(view);

            MainRoof r;
            r.roof = move(rd);
            //roofs.push_back(std::move(rd));
            roofs.push_back(move(r));
        }

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
        //StageFactory factory;


        for (int k = 0; k < roof.size(); k++) {
            //auto& rd = *roof[k];

            //
            //BufferReader reader;
            //reader.addView(*rd.view.begin());


            //Stage* cluster = factory.createStage("filters.cluster");
            //Options opts;
            //opts.add("tolerance", tolerance);                     //max distance point to be added to the cluster
            //opts.add("min_points", min_points);                    //minimum number of points in a cluster
            //cluster->setOptions(opts);
            //cluster->setInput(reader);
            //cluster->prepare(rd.table);
            //cluster->execute(rd.table);
        }
        

        std::vector<std::thread> threads;
        
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
                makeClusteredFiles(*rd.view.begin(), folderName);
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


            for (int j = 0; j < roofs[i].nsubroof.size(); j++) {

                PointViewPtr view = *roofs[i].nsubroof[j]->view.begin();

                vector<Eigen::Vector3d> points;
                for (PointId id = 0; id < view->size(); ++id) {
                    double x = view->getFieldAs<double>(Dimension::Id::X, id);
                    double y = view->getFieldAs<double>(Dimension::Id::Y, id);
                    double z = view->getFieldAs<double>(Dimension::Id::Z, id);
                    Eigen::Vector3d point(x, y, z);
                    points.push_back(point);
                }
                BestFitPlane plane = countPlane(points);
                //planes.push_back(plane);
                roofs[i].plane.push_back(plane);

                /*cout << "\n\nPLANE " << cc << endl;
                cout << "x " << plane.normal.x() << endl;
                cout << "y " << plane.normal.y() << endl;
                cout << "z " << plane.normal.z() << endl;
                cout << "d " << plane.d << endl;
                cc++;*/



            }

        }




    }
    void printPlanes(string folder) {
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




        for (int i = 0; i < roofs.size(); i++) {
            string path = folder + "/roof" + to_string(i + 1);

            filesystem::create_directories(path);

            for (int j = 0; j < roofs[i].nsubroof.size(); j++) {

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

    //TriangularMesh* buildDelaunayMesh(PointViewPtr inView)
    //{

    //    // If mesh already exists, do not run Delaunay again (would throw on duplicate mesh name).
    //    if (auto* existing = inView->mesh("delaunay2d"))
    //        return existing;

    //    PointTableRef table = inView->table();

    //    BufferReader reader;
    //    reader.addView(inView);

    //    DelaunayFilter delaunay;
    //    delaunay.setInput(reader);

    //    // Run ONLY the Delaunay stage.
    //    delaunay.prepare(table);
    //    delaunay.execute(table);

    //    // Delaunay attaches mesh to the same view (inView).
    //    inView->mesh("delaunay2d");
    //    

    //    TriangularMesh* mesh = inView->mesh("delaunay2d");

    //    return mesh;
    //}
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



    vector<PointId> hollowMesh(PointViewPtr view, TriangularMesh mesh) {

        vector<int> ids;
        for (const Triangle& t : mesh) {

            PointId a = t.m_a, b = t.m_b, c = t.m_c;
			ids.push_back(a);
			ids.push_back(b);
			ids.push_back(c);
        }
        
        cout << "SIZE " << ids.size() << endl;
        keep_only_values_with_freq_at_most_4(ids);
        cout << "SIZE NEW " << ids.size() << endl;


		remove_repeats_sorted(ids);
        cout << "SIZE NEWNEW " << ids.size() << endl;


        vector<PointId> out;
        for (int i = 0; i < ids.size(); i++) {
            PointId newId = ids[i];
			out.push_back(newId);
        }

        return out;
    }

    void attachTriangularMesh(PointViewPtr view, TriangularMesh mesh)
    {
        TriangularMesh* dst = view->createMesh("");
        *dst = std::move(mesh);
        cout << "reduced size " << dst->size() << endl;
        //return dst;
    }
    void getBorder(PointViewPtr view, vector<PointId> *input, string filename) {
 
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
            cout << "ROOF " << i << endl;


			int linescount = 0;
            for (int j = 0; j < roofs[i].plane.size() -1; j++) {
                for (int k = j + 1; k < roofs[i].plane.size(); k++) {

					linescount++;
                    Line3D line;
                    intersectPlanes(roofs[i].plane[j], roofs[i].plane[k], line);
                    line.origin.push_back(j);
                    line.origin.push_back(k);
                    cout << "Line " << linescount << " of plane " << j + 1 << " and plane " << k + 1 << endl;
                    cout << line.p.x() << "," << line.p.y() << "," << line.p.z() << endl;
					cout << line.dir.x() << "," << line.dir.y() << "," << line.dir.z() << endl << endl;
                }

                
            }
        }
	}

    void computateMeshes() {


        for (int i = 0; i < roofs.size(); i++) {

            for (int j = 0; j < roofs[i].fsubroof.size(); j++) {
               //TriangularMesh* m = buildDelaunayMesh(*roofs[i].fsubroof[j]->view.begin());
               //MeshReadyView m = buildDelaunayMesh(*roofs[i].fsubroof[j]->view.begin());
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
                vector<PointId> k = hollowMesh(*roofs[i].fsubroof[j]->view.begin(), m);
                roofs[i].jacket.push_back(k);

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
                getBorder(*roofs[i].fsubroof[j]->view.begin(), &roofs[i].jacket[j], fileName);
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

   /* pdal::PointTable table;
    pdal::PointViewPtr view = loadXYZ_withReaderText(table, "C:/C/PointReader/build/finalRoofs/roof1/subroof2.txt");
    cout << "size " << view->size() << endl;



    TriangularMesh* m = p.buildDelaunayMesh(view);
    cout << "mesh size " << m->size() << endl;
    
    p.writeMeshToPly(view, "C:/C/PointReader/build/mesh_2.ply");


    TriangularMesh nm = p.cutMesh(view, m, 0.25);
	p.attachTriangularMesh(view, nm);
    p.writeMeshToPly(view, "C:/C/PointReader/build/mesh2_2.ply");


	vector<PointId> k = p.hollowMesh(view, nm);
    p.getBorder(view, &k, "border2.txt");*/





    Stage* stage;
    stage = p.loadFile("LiDAR.laz");
    stage = p.filterClass6(stage);
    stage = p.computeNormals(stage, 16);
    stage = p.filterWalls(stage, 0.01, 20);
    stage = p.filterOutliers(stage, "statistical", 6, 0.5, 1, 4);
    stage = p.clusterPoints(stage, 1.5, 60);            //klasterizacia vsetkych bodov na hlavne strechy
    stage = p.zsmooth(stage, 1);
    p.execute(stage);



    //p.filerByZValue();

    p.makePointViewSetRoofs();
    p.makeClusteredFiles(*p.buildpoint.begin(), "roofs");


    p.clusterByPoints("mainroofs", 0.25, 60);                //klasterizacia hlavnych striech na podstrechy
	p.makeClusteredRoofsFiles("subroofs");
    p.fillSubroofs();
    p.modifySubroofs();
    p.clusterByNormals(0.01, 40);            //klasterizacia podstriech na mensie podstrechy
    p.finalizeSubroofs();
    p.printFinalRoofs("finalRoofs");



    p.clusterByPoints("nroofs", 0.25, 60);
    p.endSubroofs();
    p.printFinalRoofs22("finalRoofsTotal");

    
    p.computatePlanes();
    p.printPlanes("planes");
	//p.computateLines();



    p.computateMeshes();
    cout << "1 " << endl;
    p.printMeshes("initialMeshes");
    cout << "2 " << endl;
    p.cutRoofs(0.25);
    cout << "3 " << endl;
    p.printMeshes("cutMeshes25");
    cout << "4 " << endl;
    p.getJackets();
    p.printJackets("jackets25");    



    chrono::duration<double> elapsedd = chrono::high_resolution_clock::now() - start;
    cout << "\n\nProgram ran in " << elapsedd.count() << " seconds.\n";

}
