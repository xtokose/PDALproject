#include <pdal/StageFactory.hpp>
#include <pdal/Options.hpp>
#include <pdal/PointTable.hpp>
#include <pdal/PointView.hpp>
#include <pdal/Dimension.hpp>
#include <pdal/Stage.hpp>
#include <pdal/io/BufferReader.hpp>
#include <pdal/PointRef.hpp>

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

const double PI = 3.14159265359;
using namespace std;

struct BuildingPoint {
    double x, y, z;
    int intensity;
    int classification;
    double nx = 0, ny = 0, nz = 0;
    double curvature;
};



//void Stage::setInput(Stage& prev);

void printSchema(const pdal::PointTable& table)
{
    auto layout = table.layout();
    std::cout << "Dimensions (" << layout->dims().size() << "):\n";
    for (auto id : layout->dims())
        std::cout << "  - " << pdal::Dimension::name(id) << "\n";
}

// --- dump first N points (all dimensions) ---
void dumpFirstN(const pdal::PointTable& table,
    const pdal::PointView& v, std::size_t N = 10)
{
    auto layout = table.layout();
    const pdal::PointId n = std::min<pdal::PointId>(N, v.size());

    for (pdal::PointId i = 0; i < n; ++i)
    {
        std::cout << "Point " << i << ":\n";
        for (auto id : layout->dims())
        {
            // getFieldAs<double> is convenient for numeric dims (applies scale/offset)
            double val = v.getFieldAs<double>(id, i);
            std::cout << "  " << std::left << std::setw(18)
                << pdal::Dimension::name(id) << " = " << val << "\n";
        }
    }
}


int main() {
    _putenv_s("PROJ_LIB", "C:\\vcpkg\\installed\\x64-windows\\share\\proj");
    _putenv_s("PROJ_DATA", "C:\\vcpkg\\installed\\x64-windows\\share\\proj");
    _putenv_s("GDAL_DATA", "C:\\vcpkg\\installed\\x64-windows\\share\\gdal");
    auto start = std::chrono::high_resolution_clock::now();


    pdal::StageFactory factory;         ////stage


        //reader
    pdal::Stage* reader = factory.createStage("readers.las");
    pdal::Options opts;                 ////options
    opts.add("filename", R"(LiDAR.laz)");
    reader->setOptions(opts);

        //filtering class 6
    pdal::Stage* range = factory.createStage("filters.range");
    pdal::Options rng;
    rng.add("limits", "Classification[6:6]");
    range->setOptions(rng);
    range->setInput(*reader);

        // Normal filter
    pdal::Stage* normal = factory.createStage("filters.normal");
    pdal::Options no;
    no.add("knn", 16);                         // number of neighbors
    normal->setOptions(no);
    normal->setInput(*range);        /////////input

        //execute
    pdal::PointTable table;
    normal->prepare(table);
    reader->prepare(table);
    pdal::PointViewSet views = normal->execute(table);
    pdal::PointViewSet views_all = reader->execute(table);


        //printg total points
    std::size_t total2 = 0;
    for (auto v : views_all) total2 += v->size();
    std::cout << "Total points (LAZ): " << total2 << "\n";
    for (auto const& v : views) {
        std::cout << "Class-6 points: " << v->size() << "\n";
    }



    

    
    

        //filter points based on angle and curvature (only class-6 points)
    double raito_curvature = 0.1;                       //curv ratio
    double ratio_angle = 20;                            //degree ratio

    vector<BuildingPoint> buildpoint;
    size_t total = 0;
    double min, max, average;
    vector<double> excludedZ;
    ofstream file("buildpoints/building_points_without_" + to_string((int)(ratio_angle)) + "degrees_walls.txt");
    file << fixed << setprecision(3);                          
    double ratio = sin(ratio_angle * PI / 180.0);
    for (auto const& v : views)
    {
        total += v->size();
        for (pdal::PointId i = 0; i < v->size(); ++i)
        {
            BuildingPoint a;
            double x = v->getFieldAs<double>(pdal::Dimension::Id::X, i);
            double y = v->getFieldAs<double>(pdal::Dimension::Id::Y, i);
            double z = v->getFieldAs<double>(pdal::Dimension::Id::Z, i);
            double nx = v->getFieldAs<double>(pdal::Dimension::Id::NormalX, i);
            double ny = v->getFieldAs<double>(pdal::Dimension::Id::NormalY, i);
            double nz = v->getFieldAs<double>(pdal::Dimension::Id::NormalZ, i);
            double curv = v->getFieldAs<double>(pdal::Dimension::Id::Curvature, i);
            a.x = x;
            a.y = y;
            a.z = z;
            a.nx = nx;
            a.ny = ny;
            a.nz = nz;
            a.curvature = curv;

            if (nz > ratio) {
                if (curv < raito_curvature) {
                    buildpoint.push_back(a);
                    file << a.x << "," << a.y << "," << a.z << std::endl;
                }             
            }
            else {
                excludedZ.push_back(abs(nz));
            }
        }
    }
    file.close();

    min = *min_element(excludedZ.begin(), excludedZ.end());
    max = *max_element(excludedZ.begin(), excludedZ.end());
    average = accumulate(excludedZ.begin(), excludedZ.end(), 0.0) / excludedZ.size();
    cout << "min " << min << endl;
    cout << "max " << max << endl;
    cout << "average " << average << endl;




        //creating groups based on similar normals
    ratio_angle = 5;    // degree ratio
    vector<BuildingPoint> temporary_buildpoint = buildpoint; //temporary vector from which we will pop out points with same group
    vector<vector<BuildingPoint>> roofs;
    ratio = sin(ratio_angle * PI / 180.0);
    size_t current_temp_size = temporary_buildpoint.size();
    int size_of_group = 0, id = 1;
    while (temporary_buildpoint.size() > 1) {
        ofstream file("roofs/roof" + to_string(id) + ".txt");
        file << fixed << setprecision(3);



        vector<BuildingPoint> group;
        group.push_back(temporary_buildpoint[0]);
        for (int i = 1; i < current_temp_size - size_of_group; i++) {
            
            double dot = temporary_buildpoint[i].nx * temporary_buildpoint[0].nx + temporary_buildpoint[i].ny * temporary_buildpoint[0].ny + temporary_buildpoint[i].nz * temporary_buildpoint[0].nz;
            double magA = sqrt(temporary_buildpoint[i].nx * temporary_buildpoint[i].nx + temporary_buildpoint[i].ny * temporary_buildpoint[i].ny + temporary_buildpoint[i].nz * temporary_buildpoint[i].nz);
            double magB = sqrt(temporary_buildpoint[0].nx * temporary_buildpoint[0].nx + temporary_buildpoint[0].ny * temporary_buildpoint[0].ny + temporary_buildpoint[0].nz * temporary_buildpoint[0].nz);

            double cosTheta = dot / (magA * magB);
            if (cosTheta > 1.0) cosTheta = 1.0;
            if (cosTheta < -1.0) cosTheta = -1.0;

            double angleRad = acos(cosTheta);
            double angleDeg = angleRad * 180.0 / PI;

            if (angleDeg < ratio_angle) {

                group.push_back(temporary_buildpoint[i]);
                file << temporary_buildpoint[i].x << "," << temporary_buildpoint[i].y << "," << temporary_buildpoint[i].z << std::endl;

                temporary_buildpoint.erase(temporary_buildpoint.begin() + i);
                i--;
                size_of_group++;
            }
        }
        file.close();
        cout << "size_of_roof " << id << " is " << size_of_group << endl;
        id++;

        
        size_of_group = 0;
        temporary_buildpoint.erase(temporary_buildpoint.begin());
        current_temp_size = temporary_buildpoint.size();
        roofs.push_back(group);


    }
    if (temporary_buildpoint.size() == 1) {
        vector<BuildingPoint> group;
        group.push_back(temporary_buildpoint[0]);
        roofs.push_back(group);
    }
    





    size_t amount = 0;
    cout << "pocet group " << roofs.size() << endl;
    for (int i = 0; i < roofs.size(); i++) {

        amount = amount + roofs[i].size();

    }
    cout << "celkovy pocet bodov " << amount << endl;


    int group_number = 5;
    for (int i = 0; i < roofs[group_number].size(); i++) {

        cout << "nx   " << roofs[group_number][i].nx << " ny   " << roofs[group_number][i].ny << " nz   " << roofs[group_number][i].nz << endl;
        cout << " " << endl;
    }


    std::cout << "Class-6 points after filtering null normals in x direction: " << buildpoint.size() << "\n\n\n";
        //printing first 10 points
    for (int i = 0; i < 10; i++) {
        /*std::cout << "point " << i << endl;
        std::cout << "x: " << buildpoint[i].x << "   y: " << buildpoint[i].y << "   z: " << buildpoint[i].z << endl;
        std::cout << "nx: " << buildpoint[i].nx << "   ny: " << buildpoint[i].ny << "   nz: " << buildpoint[i].nz << endl;
        std::cout << "curvature " << buildpoint[i].curvature << endl << endl;*/

    }












    //auto v = *views.begin();                 // first view
    //printSchema(table);
    //dumpFirstN(table, *v, 10);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "\n\nProgram ran in " << elapsed.count() << " seconds.\n";
}
