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
    double ratio_curvature = 0.01;                       //curv ratio for filtering points
    double ratio_angle = 20;                            //degree ratio for filtering walls
    double ratio_Zvalue = 0;
    double delta = 0.1;                                //smoothing_parameter


    StageFactory factory;
    Stage* reader;
    Stage* range;
    Stage* cluster;
    Stage* normal;
    PointTable mainTable;



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

public:
    PointViewSet buildpoint;
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


    void loadFile() {
        auto start = chrono::high_resolution_clock::now();

        reader = factory.createStage("readers.las");
        PointTable table;
        Options opts;                                        ////options
        opts.add("filename", R"(LiDAR.laz)");
        reader->setOptions(opts);
        reader->prepare(table);
        PointViewSet all_points = reader->execute(table);

            //debug
        auto duration = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start);
        if(debug){
            cout << "\nLoading finished in " << duration.count()/1000. << " seconds\n";
            cout << "Points remaining: " << get_points_size(all_points) << endl;
        }
        
    }
    void filterClass6() {
        auto start = chrono::high_resolution_clock::now();

        range = factory.createStage("filters.range");
        Options opts;
        opts.add("limits", "Classification[6:6]");
        range->setOptions(opts);
        range->setInput(*reader);



            //debug
        auto duration = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start);
        if (debug) { cout << "\nFiltering class 6 finished in " << duration.count()/1000. << " miliseconds\n"; }
    }
    void computeNormals() {
        auto start = chrono::high_resolution_clock::now();

        normal = factory.createStage("filters.normal");
        PointTable table;
        Options opts;
        opts.add("knn", 16);                         // number of neighbors
        normal->setOptions(opts);
        normal->setInput(*range);        /////////input
        //normal->prepare(table);
        //buildpoint = normal->execute(table);



            //debug      
        if (debug) {
            auto duration = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start);
            cout << "\nComputing normals finished in " << duration.count() / 1000. << " miliseconds\n";

        }
        
    }
    void clusterPoints() {
        auto start = chrono::high_resolution_clock::now();
        cluster = factory.createStage("filters.cluster");
        Options opts;
        opts.add("tolerance", 1.5);
        opts.add("min_points", 20);
        cluster->setOptions(opts);
        cluster->setInput(*normal);
        
      
        //debug      
        if (debug) {
            auto duration = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start);
            cout << "\nClustering finished in " << duration.count() / 1000. << " miliseconds\n";
        }




        execute();
        makePointFile("points_after_filtering_class_6.txt");
    }
    void execute() {
        auto start2 = chrono::high_resolution_clock::now();
        cluster->prepare(mainTable);
        buildpoint = cluster->execute(mainTable);



        if (debug) {
            auto duration = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start2);
            cout << "\Executing point table finished in " << duration.count() / 1000. << " seconds\n";
            cout << "Points remaining: " << get_points_size(buildpoint) << endl;
        }
    }
    void filterWalls() {
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
                if (curv < ratio_curvature) {
                    outView->appendPoint(*view, i);
                }
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
    void filterOutliers() {
        auto start = chrono::high_resolution_clock::now();
        PointViewPtr view = *buildpoint.begin();


        BufferReader breader;
        breader.addView(view);
        PipelineManager manager;


        Stage& outlier = manager.makeFilter("filters.outlier", breader);
        Options opts;
        opts.add("method", "statistical");
        opts.add("mean_k", 6);
        opts.add("multiplier", 0.5);
        outlier.setOptions(opts);

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
    void smoothPoints(int iterations) {
        auto start = chrono::high_resolution_clock::now();
        PointViewPtr view = *buildpoint.begin();


        for (int j = 0; j < iterations; j++) {
            //double curv = view->getFieldAs<double>(Dimension::Id::Curvature, i);

            for (int i = 0; i < view->size(); i++) {
                double curv = view->getFieldAs<double>(Dimension::Id::Curvature, i);
                double x = view->getFieldAs<double>(Dimension::Id::X, i);
                double y = view->getFieldAs<double>(Dimension::Id::Y, i);
                double z = view->getFieldAs<double>(Dimension::Id::Z, i);
                double nx = view->getFieldAs<double>(Dimension::Id::NormalX, i);
                double ny = view->getFieldAs<double>(Dimension::Id::NormalY, i);
                double nz = view->getFieldAs<double>(Dimension::Id::NormalZ, i);

                x = x - nx * delta * curv;
                y = y - ny * delta * curv;
                z = z - nz * delta * curv;

                view->setField(Dimension::Id::X, i , x);
                view->setField(Dimension::Id::Y, i, y);
                view->setField(Dimension::Id::Z, i, z);

            }
            recomputeNormals();
        }

        


            //debug
        if (debug) {
            auto duration = chrono::duration_cast<chrono::milliseconds>(chrono::high_resolution_clock::now() - start);
            cout << "\nSmoothing points finished in " << duration.count() / 1000. << " seconds\n";
        }
        makePointFile("points_after_smoothing.txt");
    }
    void makeClusteredFiles(string folder) {
        PointViewPtr view = *buildpoint.begin();
        unordered_set<int> clusterIDs;


        for (PointId i = 0; i < view->size(); ++i)
        {
            int id = view->getFieldAs<int>(Dimension::Id::ClusterID, i);

            if (id >= 0)               // ignore noise
                clusterIDs.insert(id);
        }
        size_t numClusters = clusterIDs.size();
        auto it = max_element(clusterIDs.begin(), clusterIDs.end());
        size_t maxVal = *it;



        vector<ofstream> files;
        files.reserve(maxVal);
        for (int i = 0; i < maxVal; i++) {
            files.emplace_back(folder + "/roof" + to_string(i + 1) + ".txt");
            files[i] << fixed << setprecision(2);

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
};




int main() {
    _putenv_s("PROJ_LIB", "C:\\vcpkg\\installed\\x64-windows\\share\\proj");
    _putenv_s("PROJ_DATA", "C:\\vcpkg\\installed\\x64-windows\\share\\proj");
    _putenv_s("GDAL_DATA", "C:\\vcpkg\\installed\\x64-windows\\share\\gdal");  
    auto start = chrono::high_resolution_clock::now();

    PointReader p;
    p.loadFile();
    p.filterClass6();
    p.computeNormals();
    p.clusterPoints();
    p.filterWalls();
    //p.filerByZValue();
    p.filterOutliers();
    //p.makeClusteredFiles("roofs");
    //p.smoothPoints(40);
    //p.makeClusteredFiles("roofs_afterSmoothing");
    


    chrono::duration<double> elapsedd = chrono::high_resolution_clock::now() - start;
    cout << "\n\nProgram ran in " << elapsedd.count() << " seconds.\n";
}
