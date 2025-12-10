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

const double PI = 3.14159265359;
using namespace std;
using namespace pdal;

class PointReader {
private:
    bool debug = true;
    double ratio_Zvalue = 0;



    


    StageFactory factory;
    PointTable mainTable;
    PointTable rooftable;


    size_t get_points_size(PointViewSet& set) {
        size_t size = 0;
        for (auto v : set) size += v->size();
        return size;
    }
    void makePointFile(string filename) {
    
        ofstream file("buildpoints/" + filename);
        file << fixed << setprecision(2);
        PointViewPtr view = *buildpoint.begin();
        for (PointId i = 0; i < view->size(); ++i) {
            double x = view->getFieldAs<double>(Dimension::Id::X, i);
            double y = view->getFieldAs<double>(Dimension::Id::Y, i);
            double z = view->getFieldAs<double>(Dimension::Id::Z, i);
            file << x << "," << y << "," << z << endl;
        }
        file.close();
    };
       
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
    unordered_set<int> getRoofsIDs(PointViewPtr view){

        unordered_set<int> set;
        for (PointId i = 0; i < view->size(); ++i)
        {
            int id = view->getFieldAs<int>(Dimension::Id::ClusterID, i);

            if (id >= 0)               // ignore noise
                set.insert(id);
        }

        return set;
    }


public:
    PointViewSet buildpoint;
    PointViewSet roofs;
    void printSchema()
    {
        auto layout = mainTable.layout();
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
        PointTable table;
        Options opts;                                        ////options
        opts.add("filename", filepath);
        reader->setOptions(opts);


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
        if (debug) { cout << "\nFiltering class 6 finished in " << duration.count()/1000. << " miliseconds\n"; }
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
        //normal->prepare(table);
        //buildpoint = normal->execute(table);



            //debug      
        if (debug) {
            auto duration = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start);
            cout << "\nComputing normals finished in " << duration.count() / 1000. << " miliseconds\n";

        }
        return normal;
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
        if (debug) {
            auto duration = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start);
            cout << "\nClustering finished in " << duration.count() / 1000. << " miliseconds\n";
        }
        return cluster;
    }
    Stage* zsmooth(Stage* input, size_t radius) {
        
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

        return assign;
    }

    void execute(Stage* input) {
        auto start2 = chrono::high_resolution_clock::now();
        input->prepare(mainTable);
        buildpoint = input->execute(mainTable);



        if (debug) {
            auto duration = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start2);
            cout << "\Executing point table finished in " << duration.count() / 1000. << " seconds\n";
            cout << "Points remaining: " << get_points_size(buildpoint) << endl;
        }
        makePointFile("points_after_filtering_class_6.txt");
    }
    void filterWalls(double ratio_curvature, double ratio_angle) {
        auto start = chrono::high_resolution_clock::now();
        PointViewSet temporaryset;
        PointViewPtr view = *buildpoint.begin();
        PointViewPtr outView = view->makeNew();
        vector<double> excludedZ;
        double ratio = sin(ratio_angle * PI / 180.0);
        
        for (PointId i = 0; i < view->size(); ++i) {

            double x = view->getFieldAs<double>(Dimension::Id::X, i);
            double y = view->getFieldAs<double>(Dimension::Id::Y, i);
            double z = view->getFieldAs<double>(Dimension::Id::Z, i);
            double nx = view->getFieldAs<double>(Dimension::Id::NormalX, i);
            double ny = view->getFieldAs<double>(Dimension::Id::NormalY, i);
            double nz = view->getFieldAs<double>(Dimension::Id::NormalZ, i);
            double curv = view->getFieldAs<double>(Dimension::Id::Curvature, i);

            if (nz > ratio) {
                //if (curv < ratio_curvature) {
                //    outView->appendPoint(*view, i);
                //}

                outView->appendPoint(*view, i);
            }
            else {
                excludedZ.push_back(abs(nz));
            }

        }
        temporaryset.insert(outView);
        buildpoint = temporaryset;

        


        //debug       
        if (debug) {
            auto duration = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start);
            cout << "\nWall filtering finished in " << duration.count() / 1000. << " seconds\n";
            cout << "Points remaining: " << get_points_size(buildpoint) << endl;
            double min = *min_element(excludedZ.begin(), excludedZ.end());
            double max = *max_element(excludedZ.begin(), excludedZ.end());
            double average = accumulate(excludedZ.begin(), excludedZ.end(), 0.0) / excludedZ.size();
            /*cout << "min " << min << endl;
            cout << "max " << max << endl;
            cout << "average " << average << endl;*/

        }
        makePointFile("points_after_wall_filter_" + to_string((int)(ratio_angle)) + "degrees.txt");
    }
    void filerByZValue() {
        auto start = chrono::high_resolution_clock::now();
        PointViewPtr view = *buildpoint.begin();
        PointViewSet temporaryset;
        PointViewPtr outView = view->makeNew();

        

        for (PointId i = 0; i < view->size(); ++i) {
            bool append = true;
            for (PointId j = 0; j < view->size(); ++j) {
                double xi = view->getFieldAs<double>(Dimension::Id::X, i);
                double yi = view->getFieldAs<double>(Dimension::Id::Y, i);
                double zi = view->getFieldAs<double>(Dimension::Id::Z, i);

                double xj = view->getFieldAs<double>(Dimension::Id::X, j);
                double yj = view->getFieldAs<double>(Dimension::Id::Y, j);
                double zj = view->getFieldAs<double>(Dimension::Id::Z, j);

                
                if (((abs(xi - xj) < ratio_Zvalue) &&
                    (abs(yi - yj) < ratio_Zvalue) &&
                    (zj > zi))) {
                    j = view->size();
                    append = false;
                }
            }
            if (append) {
                outView->appendPoint(*view, i);
            }
            
        }

        temporaryset.insert(outView);
        buildpoint = temporaryset;

        //debug      
        if (debug) {
            auto duration = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start);
            cout << "\nFiltering by Z value finished in " << duration.count() / 1000. << " seconds\n";
            cout << "Points remaining: " << get_points_size(buildpoint) << endl;

        }
        makePointFile("points_after_Zfilter.txt");
    }
    void filterOutliers(string method, size_t mean_k, float multiplier, float radius, size_t min_k) {
        auto start = chrono::high_resolution_clock::now();
        PointViewPtr view = *buildpoint.begin();


        BufferReader breader;
        breader.addView(view);
        PipelineManager manager;



        Stage& outlier = manager.makeFilter("filters.outlier", breader);
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
        manager.execute();

        PointViewSet result = manager.views();
        buildpoint = result;

        //debug      
        if (debug) {
            auto duration = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start);
            cout << "\nFiltering outliers finished in " << duration.count() / 1000. << " seconds\n";
            cout << "Points remaining: " << get_points_size(buildpoint) << endl;

        }
        makePointFile("points_after_outliers_filter.txt");
    }   
    void makeClusteredFiles(bool smoothed, string folder) {
        PointViewPtr view = *buildpoint.begin();

            //removing existing
        filesystem::path path = folder;
        for (const auto& entry : filesystem::directory_iterator(path)) {
            if (!filesystem::is_regular_file(entry.status()))
                continue; // skip subdirs etc.

            if (entry.path().extension() == ".txt") {
                filesystem::remove(entry.path());
            }
        }

        unordered_set<int> clusterIDs = getRoofsIDs(view);
        size_t numClusters = clusterIDs.size();


        //auto it = max_element(clusterIDs.begin(), clusterIDs.end());
        //size_t maxVal = *it;
        //cout << "maxVal " << maxVal << endl;
        //size_t numClusters = 30;
        //size_t maxVal = 30;


        vector<ofstream> files;
        files.reserve(numClusters);

        for (int i = 0; i < numClusters; i++) {
            files.emplace_back(folder + "/roof" + to_string(i + 1) + ".txt");
            files[i] << fixed << setprecision(3);
        }


        for (PointId i = 0; i < view->size(); ++i) {
            int id = view->getFieldAs<int>(Dimension::Id::ClusterID, i);
            double x = view->getFieldAs<double>(Dimension::Id::X, i);
            double y = view->getFieldAs<double>(Dimension::Id::Y, i);
            double z = view->getFieldAs<double>(Dimension::Id::Z, i);
            
            if (id > 0) {
                files[id - 1] << x << "," << y << "," << z << endl;
            }
        }

        for (int i = 0; i < numClusters; i++) {
            files[i].close();
        }

        cout << "Number of clusters: " << numClusters << endl;

        
    }
    void makePointViewSetRoofs() {
        PointViewPtr old_view = *buildpoint.begin();

        
        PointLayoutPtr layout = rooftable.layout();
        layout->registerDim(Dimension::Id::X);
        layout->registerDim(Dimension::Id::Y);
        layout->registerDim(Dimension::Id::Z);
        layout->registerDim(Dimension::Id::NormalX);
        layout->registerDim(Dimension::Id::NormalY);
        layout->registerDim(Dimension::Id::NormalZ);
        layout->registerDim(Dimension::Id::Curvature);




        unordered_set<int> clusterIDs = getRoofsIDs(old_view);
        for (int i = 0; i < clusterIDs.size(); i++) {


            PointViewPtr view(new PointView(rooftable));
            for (int j = 0; j < old_view->size(); j++) {

                size_t clusterID = old_view->getFieldAs<double>(Dimension::Id::ClusterID, j);
                if (clusterID == (i + 1)) {
                    double x = old_view->getFieldAs<double>(Dimension::Id::X, j);
                    double y = old_view->getFieldAs<double>(Dimension::Id::Y, j);
                    double z = old_view->getFieldAs<double>(Dimension::Id::Z, j);
                    double nx = old_view->getFieldAs<double>(Dimension::Id::NormalX, j);
                    double ny = old_view->getFieldAs<double>(Dimension::Id::NormalY, j);
                    double nz = old_view->getFieldAs<double>(Dimension::Id::NormalZ, j);
                    double curv = old_view->getFieldAs<double>(Dimension::Id::Curvature, j);


                    view->setField(Dimension::Id::X, j, x);
                    view->setField(Dimension::Id::Y, j, y);
                    view->setField(Dimension::Id::Z, j, z);
                    view->setField(Dimension::Id::NormalX, j, nx);
                    view->setField(Dimension::Id::NormalY, j, ny);
                    view->setField(Dimension::Id::NormalZ, j, nz);
                    view->setField(Dimension::Id::Curvature, j, curv);
                }
            }
            roofs.insert(view);


        }


        

        for (PointId i = 0; i < old_view->size(); ++i) {

            /*view->setField(Dimension::Id::X, i, x);
            view->setField(Dimension::Id::Y, i, y);
            view->setField(Dimension::Id::Z, i, z);*/
        }
    }
    void clusterRoofs(float tolerance, size_t min_points) {
        StageFactory factory;
        PointViewSet clusteredRoofs;


        PointViewSet temp;
        for (PointViewPtr view : roofs)
        {
            if (!view || view->empty())
                continue;

            BufferReader reader;
            reader.addView(view);

            Stage* dbscanStage = factory.createStage("filters.dbscan");
            if (!dbscanStage)
                throw pdal_error("Failed to create filters.dbscan stage");

            Options opts;
            opts.add("min_points", 10);
            opts.add("eps", 2.0);
            opts.add("dimensions", "X,Y,Z");

            dbscanStage->setOptions(opts);
            dbscanStage->setInput(reader);

            // Use the SAME table that was used to create the views:
            dbscanStage->prepare(rooftable);
            PointViewSet clusteredViews = dbscanStage->execute(rooftable);

            temp.insert(*clusteredViews.begin());
        }

        //for (PointViewPtr srcView : roofs)
        //{
        //    if (!srcView || srcView->empty())
        //        continue;

        //    // 1) Fresh table JUST for this roof

        //    // 2) Register the dimensions DBSCAN will use
        //    PointLayoutPtr layout = rooftable.layout();
        //    layout->registerDim(Dimension::Id::X);
        //    layout->registerDim(Dimension::Id::Y);
        //    layout->registerDim(Dimension::Id::Z);

        //    // 3) Set up BufferReader + DBSCAN
        //    BufferReader reader;

        //    Stage* dbscanStage = factory.createStage("filters.dbscan");
        //    if (!dbscanStage)
        //        throw pdal_error("Failed to create filters.dbscan stage");

        //    Options opts;
        //    opts.add("min_points", 10);
        //    opts.add("eps", 2.0);
        //    opts.add("dimensions", "X,Y,Z");   // these now exist in layout

        //    dbscanStage->setOptions(opts);
        //    dbscanStage->setInput(reader);

        //    // 4) Let DBSCAN add ClusterID and finalize the layout
        //    dbscanStage->prepare(rooftable);

        //    // 5) Create a new view on this table and copy XYZ from the original roof view
        //    PointViewPtr view(new PointView(rooftable));

        //    for (PointId i = 0; i < srcView->size(); ++i)
        //    {
        //        PointId j = view->size();

        //        double x = srcView->getFieldAs<double>(Dimension::Id::X, i);
        //        double y = srcView->getFieldAs<double>(Dimension::Id::Y, i);
        //        double z = srcView->getFieldAs<double>(Dimension::Id::Z, i);

        //        view->setField(Dimension::Id::X, j, x);
        //        view->setField(Dimension::Id::Y, j, y);
        //        view->setField(Dimension::Id::Z, j, z);
        //        // ClusterID will be filled by DBSCAN during execute().
        //    }

        //    reader.addView(view);

        //    // 6) Run DBSCAN – this writes ClusterID into `view`
        //    PointViewSet out = dbscanStage->execute(rooftable);

        //    // Use the returned views (normally contains `view`)
        //    clusteredRoofs.insert(out.begin(), out.end());
        //}

        //roofs.swap(clusteredRoofs);







        cout << "temp    " << temp.size() << endl;

        int hh = 0;
        for (PointViewPtr kkk : temp) {
            cout << "cluster " << hh << endl;

            for (int j = 0; j < kkk->size(); j++) {
                
                int ii = kkk->getFieldAs<double>(Dimension::Id::X, j);
                cout << " " << ii << endl;
            }
            hh++;
        }
    }
    void makeSmallRooofsFile() {
        string folderName = "folderujem";
        if (!filesystem::exists(folderName)) {
            if (filesystem::create_directories(folderName)) {
                std::cout << "Folder created: " << folderName << '\n';
            }
            else {
                std::cout << "Failed to create folder\n";
            }
        }



        for (PointViewPtr view:roofs) {





        }

    }
};




int main() {
    _putenv_s("PROJ_LIB", "C:\\vcpkg\\installed\\x64-windows\\share\\proj");
    _putenv_s("PROJ_DATA", "C:\\vcpkg\\installed\\x64-windows\\share\\proj");
    _putenv_s("GDAL_DATA", "C:\\vcpkg\\installed\\x64-windows\\share\\gdal");  
    auto start = chrono::high_resolution_clock::now();

    PointReader p;
    Stage* stage;
    stage = p.loadFile("LiDAR.laz");
    stage = p.filterClass6(stage);
    stage = p.computeNormals(stage, 16);
    stage = p.clusterPoints(stage, 1.5,20);
    stage = p.zsmooth(stage,1);
    p.execute(stage);


    p.filterWalls(0.01, 15);
    //p.filerByZValue();
    p.filterOutliers("statistical", 6, 0.5, 1, 4);


    //p.makeClusteredFiles(false, "roofs");
    p.makeClusteredFiles(true , "roofs_afterSmoothing");




    //p.makePointViewSetRoofs();
    //p.clusterRoofs(1.5, 20);

    chrono::duration<double> elapsedd = chrono::high_resolution_clock::now() - start;
    cout << "\n\nProgram ran in " << elapsedd.count() << " seconds.\n";
}
