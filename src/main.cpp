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



    

    
    

        // count normals (only class-6 points)
    std::vector<BuildingPoint> buildpoint;
    std::size_t total = 0;
    std::ofstream file("points.txt");
    file << std::fixed << std::setprecision(3);
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

            file << a.x << "," << a.y << "," << a.z << std::endl;


            double ratio = 0.01;
            if (nx > ratio) {
                buildpoint.push_back(a);
            }

        }
    }
    file.close();

    std::cout << "Class-6 points after filtering null normals in x direction: " << buildpoint.size() << "\n\n\n";


        //printing first 10 points
    for (int i = 0; i < 10; i++) {
        std::cout << "point " << i << endl;
        std::cout << "x: " << buildpoint[i].x << "   y: " << buildpoint[i].y << "   z: " << buildpoint[i].z << endl;
        std::cout << "nx: " << buildpoint[i].nx << "   ny: " << buildpoint[i].ny << "   nz: " << buildpoint[i].nz << endl;
        std::cout << "curvature " << buildpoint[i].curvature << endl << endl;

    }


    auto v = *views.begin();                 // first view
    //printSchema(table);
    //dumpFirstN(table, *v, 10);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "\n\nProgram ran in " << elapsed.count() << " seconds.\n";
}
