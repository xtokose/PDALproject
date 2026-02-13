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


#include <pdal/PointView.hpp>
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

const double PI = 3.14159265359;
using namespace std;
using namespace pdal;

struct BestFitPlane
{
    
    Eigen::Vector3d normal;   // unit-length (A,B,C)
    double d;                 // D in Ax + By + Cz + D = 0
    Eigen::Vector3d centroid; // point on plane
};

struct RoofData {
    pdal::PointTable table;
    pdal::PointViewSet view; // contains the PointViewPtrs made from `table`
};



class PointReader {
private:
    bool debug = true;
    double ratio_Zvalue = 0;



    StageFactory factory;
    PointTable mainTable;

    
    //vector<PointTable> rooftable;
    size_t numClusters;

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
    unordered_set<double> getRoofsIDs(PointViewPtr view){


        unordered_set<double> set;
        for (PointId i = 0; i < view->size(); ++i)
        {
            double id = view->getFieldAs<double>(Dimension::Id::ClusterID, i);


            if (id >= 0)               // ignore noise
                set.insert(id);
        }

        return set;
    }


public:
    PointViewSet buildpoint;
    //vector<PointViewSet> roofs;
   // std::deque<pdal::PointTable> rooftable;
    vector<BestFitPlane> planes;
    std::vector<std::unique_ptr<RoofData>> roofs;


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
    void makeClusteredFiles(const PointViewPtr& view , string folder) {
        
        /*int cc = 1;
        cout << "cluster " << cc << endl;
        cc++;*/

    

        filesystem::create_directories(folder);

            //removing existing
        filesystem::path path = folder;

        for (const auto& entry : filesystem::directory_iterator(path)) {
            if (!filesystem::is_regular_file(entry.status()))
                continue; // skip subdirs etc.

            if (entry.path().extension() == ".txt") {
                filesystem::remove(entry.path());
            }
        }

        

        unordered_set<double> clusterIDs = getRoofsIDs(view);

        numClusters = clusterIDs.size();
        cout << "numClusters " << numClusters << endl;

        


        vector<ofstream> files;
        files.reserve(numClusters);

        for (int i = 0; i < numClusters; i++) {
            files.emplace_back(folder + "/roof" + to_string(i + 1) + ".txt");
            files[i] << fixed << setprecision(3);
        }
        
        cout << "sum DONE " << endl;
        for (PointId i = 0; i < view->size(); ++i) {
            double id = view->getFieldAs<double>(Dimension::Id::ClusterID, i);
            double x = view->getFieldAs<double>(Dimension::Id::X, i);
            double y = view->getFieldAs<double>(Dimension::Id::Y, i);
            double z = view->getFieldAs<double>(Dimension::Id::Z, i);
            


            if ((id > 0) && (id < 10000)) {
                files[id - 1] << x << "," << y << "," << z << endl;
            }
        }
        cout << "DONE " << endl;
        for (int i = 0; i < numClusters; i++) {
            files[i].close();
        }

        cout << "Number of clusters: " << numClusters << endl;

        
    }
    void makePointViewSetRoofs() {
        PointViewPtr old_view = *buildpoint.begin();

		int sum = 0;

        for (int k = 0; k < numClusters; k++) {
            //rooftable.emplace_back();
            //PointTable& t = rooftable.back();

           // RoofData rd;

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
            cout << "view->size() " << view->size() << endl;
			sum += view->size();
            PointViewSet temp;
            temp.insert(view);

            rd->view.insert(view);
            roofs.push_back(std::move(rd));



            /////////////
            
    //        unordered_set<int> clusterIDs = getRoofsIDs(old_view);
    //        cout << "clusterIDs.size() " << clusterIDs.size() << endl;
    //        for (int i = 0; i < clusterIDs.size(); i++) {

    //            int count = 0;
    //            PointViewPtr view(new PointView(rooftable[k]));
    //            for (int j = 0; j < old_view->size(); j++) {

    //                size_t clusterID = old_view->getFieldAs<double>(Dimension::Id::ClusterID, j);
    //                if (clusterID == (i + 1)) {


    //                    double x = old_view->getFieldAs<double>(Dimension::Id::X, j);
    //                    double y = old_view->getFieldAs<double>(Dimension::Id::Y, j);
    //                    double z = old_view->getFieldAs<double>(Dimension::Id::Z, j);
    //                    double nx = old_view->getFieldAs<double>(Dimension::Id::NormalX, j);
    //                    double ny = old_view->getFieldAs<double>(Dimension::Id::NormalY, j);
    //                    double nz = old_view->getFieldAs<double>(Dimension::Id::NormalZ, j);
    //                    double curv = old_view->getFieldAs<double>(Dimension::Id::Curvature, j);


    //                    view->setField(Dimension::Id::X, count, x);
    //                    view->setField(Dimension::Id::Y, count, y);
    //                    view->setField(Dimension::Id::Z, count, z);
    //                    view->setField(Dimension::Id::NormalX, count, nx);
    //                    view->setField(Dimension::Id::NormalY, count, ny);
    //                    view->setField(Dimension::Id::NormalZ, count, nz);
    //                    view->setField(Dimension::Id::Curvature, count, curv);
    //                    count++;
    //                }

    //            }

    //            //cout << "view->size() " << view->size() << endl;
    //            PointViewSet temp;
				//temp.insert(view);
				//roofs.push_back(temp);
    //        }


        }
   

       


        


    }
    void clusterRoofs(float tolerance, size_t min_points) {
        StageFactory factory;
        PointViewSet clusteredRoofs;

        for (int k = 0; k < numClusters; k++) {
            auto& rd = *roofs[k];

            
            BufferReader reader;
            reader.addView(*rd.view.begin());


            Stage* cluster = factory.createStage("filters.cluster");
            Options opts;
            opts.add("tolerance", tolerance);                     //max distance point to be added to the cluster
            opts.add("min_points", min_points);                    //minimum number of points in a cluster
            cluster->setOptions(opts);
            cluster->setInput(reader);
            cluster->prepare(rd.table);
   

            PointViewSet temp = cluster->execute(rd.table);

            clusteredRoofs.insert(temp.begin(), temp.end());
        }
        



        //BufferReader reader;
        //reader.addView(*roofs.begin());

        /*size_t idx = 0;
        for (const auto& v : roofs)
        {
            auto cid = v->layout()->findDim("ClusterID");
            std::cout << "view " << idx++ << ": "
                << (cid == pdal::Dimension::Id::Unknown ? "MISSING ClusterID" : "has ClusterID")
                << "\n";
        }*/


        //Stage* cluster = factory.createStage("filters.cluster");
        //Options opts;
        //opts.add("tolerance", tolerance);                     //max distance point to be added to the cluster
        //opts.add("min_points", min_points);                    //minimum number of points in a cluster
        //cluster->setOptions(opts);
        //cluster->setInput(reader);

        //cluster->prepare(rooftable);
        //PointViewSet temp = cluster->execute(rooftable);
        //clusteredRoofs.insert(temp.begin(), temp.end());


		

        //for (const PointViewPtr& v : roofs) {

        //    cout << "\n\n\n " << endl;
        //    printSchema();


        //    BufferReader reader;
        //    reader.addView(v);

        //    //pdal::PointTableRef t = v->table();

        //    Stage* cluster = factory.createStage("filters.cluster");
        //    Options opts;
        //    opts.add("tolerance", tolerance);                     //max distance point to be added to the cluster
        //    opts.add("min_points", min_points);                    //minimum number of points in a cluster
        //    cluster->setOptions(opts);
        //    cluster->setInput(reader);

        //    
        //    cluster->prepare(rooftable);
        //    PointViewSet temp = cluster->execute(rooftable);
        //    clusteredRoofs.insert(temp.begin(), temp.end());
        //}



        //int hh = 1;
        //for (const PointViewPtr& v : clusteredRoofs) {

        //    cout << "CLUSTER " << hh << endl;
        //    for (int i = 0; i < v->size(); i++) {
        //         int id = v->getFieldAs<double>(Dimension::Id::ClusterID, i);
        //         cout << "ID " << id << endl;
        //    }
        //    hh++;   
        //}
    }
    void makeClusteredRoofsFiles(string folder) {

        /*string folderName = "folderujem";
        if (!filesystem::exists(folderName)) {
            if (filesystem::create_directories(folderName)) {
                std::cout << "Folder created: " << folderName << '\n';
            }
            else {
                std::cout << "Failed to create folder\n";
            }
        }*/


        std::error_code ec;
        if (filesystem::exists(folder, ec))
        {
            filesystem::remove_all(folder, ec);   // recursive
            if (ec) throw filesystem::filesystem_error("remove_all failed", folder, ec);
        }

       // filesystem::create_directories(folder);
        filesystem::create_directories(folder, ec);
        if (ec) throw filesystem::filesystem_error("create_directories failed", folder, ec);

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

            auto& rd = *roofs[i];
            PointViewPtr ff = *rd.view.begin();
            cout << "view->size()  " << ff->size() << endl;

            if (ff->size() != 0) {
                makeClusteredFiles(*rd.view.begin(), folderName);
                count++;
            }
            
            

        }


        //for (const auto& v : roofs) {
        //    

        //    string folderName = folder + "/roof" + to_string(count);
        //    //makeClusteredFiles(v,folderName);
        //    count++;

        //}



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


        Eigen::Vector3d b1(0, 0, 1);
        Eigen::Vector3d b2(2, 0, 2);
        Eigen::Vector3d b3(1, 3, 0);
        Eigen::Vector3d b4(1, 0, 1.5);
        Eigen::Vector3d b5(1, 3, 1);
		vector<Eigen::Vector3d> testPoints = { b1, b2, b3, b4, b5 };
        BestFitPlane plane = countPlane(testPoints);

        cout << "x " << plane.normal.x() << endl;
        cout << "y " << plane.normal.y() << endl;
        cout << "z " << plane.normal.z() << endl;
        cout << "d " << plane.d << endl;

        /*int cc = 1;
        for (PointViewPtr view : roofs) {

			vector<Eigen::Vector3d> points;
            for (PointId i = 0; i < view->size(); ++i) {
                double x = view->getFieldAs<double>(Dimension::Id::X, i);
                double y = view->getFieldAs<double>(Dimension::Id::Y, i);
                double z = view->getFieldAs<double>(Dimension::Id::Z, i);
				Eigen::Vector3d point(x, y, z);
				points.push_back(point);
            }
            BestFitPlane plane = countPlane(points);
			planes.push_back(plane);

            cout << "\n\nPLANE " << cc << endl;
            cout << "x " << plane.normal.x() << endl;
            cout << "y " << plane.normal.y() << endl;
            cout << "z " << plane.normal.z() << endl;
            cout << "d " << plane.d << endl;
            cc++;

        }*/


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
    p.makeClusteredFiles(*p.buildpoint.begin(), "kotlety");


    


    p.makePointViewSetRoofs();
   

    p.clusterRoofs(0.5, 5);
    
	p.makeClusteredRoofsFiles("subroofs");


     //p.computatePlanes();
    //Eigen::Matrix3d C;
    //C << 2.0, 0.0, 1.0,
    //    0.0, 6.0, -3.0,
    //    1.0, -3.0, 2.0;
    //std::cout << C << endl;

    //Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(C);


    // //Eigen returns them in ascending order for SelfAdjointEigenSolver
    // std::cout << "Eigenvalues:\n" << solver.eigenvalues() << "\n";
    // std::cout << "Eigenvectors:\n" << solver.eigenvectors() << "\n";



    

    chrono::duration<double> elapsedd = chrono::high_resolution_clock::now() - start;
    cout << "\n\nProgram ran in " << elapsedd.count() << " seconds.\n";
}
