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

    std::cout << "test";

    pdal::StageFactory factory;     ////stage


    //reader
    pdal::Stage* reader = factory.createStage("readers.las");
    pdal::Options opts;     ////options
    opts.add("filename", R"(LiDAR.laz)"); // <-- .laz file
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
    no.add("knn", 16);                         // number of neighbors (default 8)
    // Optional:
    // no.add("always_up", true);              // flip toward +Z (default true)
    // no.add("refine", true);                 // MST propagation for more consistency
    // no.add("viewpoint", "POINT Z (0 0 10)");// flip toward a 3D viewpoint (WKT/GeoJSON)
    normal->setOptions(no);
    normal->setInput(*range); /////////input

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


    std::ofstream file("points.txt");


    // count normals (only class-6 points exist now)
    std::vector<BuildingPoint> buildpoint;
    std::size_t total = 0, ok = 0;
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
            a.x = x;
            a.y = y;
            a.z = z;
            a.nx = nx;
            a.ny = ny;
            a.nz = nz;

            file << x << "," << y << "," << z << endl;

            double ratio = 0.01;
            if (nx > ratio) {
                buildpoint.push_back(a);
            }


            //std::cout << "\npoint " << i << std::endl;
           // std::cout << "nx " << nx << std::endl;

            const double len2 = nx * nx + ny * ny + nz * nz;
            if (std::isfinite(nx) && std::isfinite(ny) && std::isfinite(nz) && len2 > 0.0) ++ok;
        }
    }
    file.close();

    //std::cout << "Class-6 points: " << total << "\n"
    //    << "Class-6 normals valid: " << ok << "\n";


    std::cout << " " << buildpoint.size() << std::endl;
    for (int i = 0; i < 100; i++) {
        std::cout << "point " << i << " ---    x " << buildpoint[i].x << std::endl;
        std::cout << " nx " << buildpoint[i].nx << endl << endl;

    }



    //loading all points from .las
    {
        ///////////

    //pdal::StageFactory factory;

    //auto reader = factory.createStage("readers.las");
    //pdal::Options opts;
    //opts.add("filename", "Mracno_bodov_vyrez.las");   // put your LAS path here
    //reader->setOptions(opts);




    //pdal::PointTable table;
    //reader->prepare(table);     //metadata    
    //pdal::PointViewSet views = reader->execute(table);  //reading data


    ////suma
    //std::size_t total = 0;
    //for (auto v : views) total += v->size();
    //std::cout << "Total points: " << total << "\n";


    /////////////////
    }


    auto v = *views.begin();                 // first view
    //printSchema(table);
    //dumpFirstN(table, *v, 10);

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Program ran in " << elapsed.count() << " seconds.\n";
}
