//#include <pdal/StageFactory.hpp>
//#include <pdal/Options.hpp>
//#include <pdal/PointTable.hpp>
//#include <pdal/PointView.hpp>
//#include <pdal/Dimension.hpp>
//#include <pdal/Stage.hpp>
//#include <pdal/io/BufferReader.hpp>
//#include <pdal/PointRef.hpp>
//
//#include <cstdlib>
//#include <iomanip>
//
//#include <vector>
//#include <string>
//#include <iostream>
//#include <cmath>
//
//auto fail = [](const std::string& msg) { std::cerr << "[ERR] " << msg << "\n"; std::exit(1); };
//
//struct BuildingPoint {
//    double x, y, z;
//    int intensity;
//    int classification;
//    double nx = 0, ny = 0, nz = 0;
//};
//
//
//
////void Stage::setInput(Stage& prev);
//
//void printSchema(const pdal::PointTable& table)
//{
//    auto layout = table.layout();
//    std::cout << "Dimensions (" << layout->dims().size() << "):\n";
//    for (auto id : layout->dims())
//        std::cout << "  - " << pdal::Dimension::name(id) << "\n";
//}
//
//// --- dump first N points (all dimensions) ---
//void dumpFirstN(const pdal::PointTable& table,
//    const pdal::PointView& v, std::size_t N = 10)
//{
//    auto layout = table.layout();
//    const pdal::PointId n = std::min<pdal::PointId>(N, v.size());
//
//    for (pdal::PointId i = 0; i < n; ++i)
//    {
//        std::cout << "Point " << i << ":\n";
//        for (auto id : layout->dims())
//        {
//            // getFieldAs<double> is convenient for numeric dims (applies scale/offset)
//            double val = v.getFieldAs<double>(id, i);
//            std::cout << "  " << std::left << std::setw(18)
//                << pdal::Dimension::name(id) << " = " << val << "\n";
//        }
//    }
//}
//
//
//int main() {
//    _putenv_s("PROJ_LIB", "C:\\vcpkg\\installed\\x64-windows\\share\\proj");
//    _putenv_s("PROJ_DATA", "C:\\vcpkg\\installed\\x64-windows\\share\\proj");
//    _putenv_s("GDAL_DATA", "C:\\vcpkg\\installed\\x64-windows\\share\\gdal");
//    pdal::StageFactory factory;     ////stage
//
//
//        //reader
//    auto reader = factory.createStage("readers.las");
//    pdal::Options opts;     ////options
//    opts.add("filename", R"(LiDAR.laz)"); // <-- .laz file
//    reader->setOptions(opts);
//
//    
//
//        //loading .las
//    pdal::PointTable table;         ////pointTable
//    reader->prepare(table);
//    pdal::PointViewSet views = reader->execute(table);
//    
//        //printg total points
//    std::size_t total2 = 0;
//    for (auto v : views) total2 += v->size();
//    std::cout << "Total points (LAZ): " << total2 << "\n";
//    
//
//    auto layout = table.layout(); // rovnaké rozloženie
//
//
//    pdal::PointViewPtr viewClass6 = std::make_shared<pdal::PointView>(table);
//    pdal::Dimension::Id classDim = pdal::Dimension::Id::Classification;
//
//
//    for (auto const& view : views)
//    {
//        for (pdal::PointId i = 0; i < view->size(); ++i)
//        {
//            uint8_t classification = view->getFieldAs<uint8_t>(classDim, i);
//            if (classification == 6)
//            {
//                viewClass6->appendPoint(*view, i);
//            }
//        }
//    }
//    std::cout << "Poèet bodov s classification = 6: " << viewClass6->size() << "\n";
//
//    
//
//
//        // Normal filter
//    pdal::BufferReader bufReader;
//    bufReader.addView(viewClass6);
//    auto normal = factory.createStage("filters.normal");
//    pdal::Options nopts;
//    nopts.add("knn", 16);
//    // alternatíva: nopts.add("radius", 0.5);  // v jednotkách dát (napr. metre)
//    // nopts.add("viewpoint", "0,0,100"); // volite¾né – orientuje normály smerom k danému bodu (x,y,z)
//    // nopts.add("always_up", true);      // volite¾né – drží Z zložku normály “nahor” pri ambiguite
//    normal->setOptions(nopts);
//    normal->setInput(bufReader);
//
//
//    pdal::PointTable tableNorm;
//
//    try {
//        normal->prepare(tableNorm);
//        auto viewsNorm = normal->execute(tableNorm);
//        if (viewsNorm.empty()) fail("filters.normal returned no views");
//
//        auto out = *viewsNorm.begin();
//        std::cout << "viewsNorm.size() = " << viewsNorm.size()
//            << ", out->size() = " << out->size() << "\n";
//
//        bool hasX = out->hasDim(pdal::Dimension::Id::X);
//        bool hasY = out->hasDim(pdal::Dimension::Id::Y);
//        bool hasZ = out->hasDim(pdal::Dimension::Id::Z);
//        bool hasNx = out->hasDim(pdal::Dimension::Id::NormalX);
//        bool hasNy = out->hasDim(pdal::Dimension::Id::NormalY);
//        bool hasNz = out->hasDim(pdal::Dimension::Id::NormalZ);
//
//        std::cout << "Dims  X:" << hasX << " Y:" << hasY << " Z:" << hasZ
//            << " | Nx:" << hasNx << " Ny:" << hasNy << " Nz:" << hasNz << "\n";
//
//        if (!hasNx || !hasNy || !hasNz)
//            fail("filters.normal did not register NormalX/Y/Z — stage likely not running with Eigen");
//
//        // vypíš prvé body
//        std::size_t toPrint = std::min<std::size_t>(10, out->size());
//        std::size_t nzCount = 0;
//        for (pdal::PointId i = 0; i < out->size(); ++i)
//        {
//            float nx = out->getFieldAs<float>(pdal::Dimension::Id::NormalX, i);
//            float ny = out->getFieldAs<float>(pdal::Dimension::Id::NormalY, i);
//            float nz = out->getFieldAs<float>(pdal::Dimension::Id::NormalZ, i);
//            if (nx != 0.f || ny != 0.f || nz != 0.f) ++nzCount;
//
//            if (i < toPrint)
//            {
//                double x = out->getFieldAs<double>(pdal::Dimension::Id::X, i);
//                double y = out->getFieldAs<double>(pdal::Dimension::Id::Y, i);
//                double z = out->getFieldAs<double>(pdal::Dimension::Id::Z, i);
//                std::cout << "i=" << i << "  XYZ=(" << x << "," << y << "," << z
//                    << ")  N=(" << nx << "," << ny << "," << nz << ")\n";
//            }
//        }
//        std::cout << "Nenulové normály: " << nzCount << " / " << out->size() << "\n";
//    }
//    catch (const pdal::pdal_error& e)
//    {
//        std::cerr << "PDAL error: " << e.what() << "\n";
//        std::exit(1);
//    }
//
//    ///////////////
//    //normal->prepare(tableNorm);
//    //pdal::PointViewSet viewsNorm = normal->execute(tableNorm);
//
//    //std::cout << "viewsNorm.size() = " << viewsNorm.size() << "\n";
//    //auto out = *viewsNorm.begin();
//    //std::cout << "out->size() = " << out->size() << "\n";
//
//    //// musia existova NormalX/Y/Z
//    //bool hasNx = out->hasDim(pdal::Dimension::Id::NormalX);
//    //bool hasNy = out->hasDim(pdal::Dimension::Id::NormalY);
//    //bool hasNz = out->hasDim(pdal::Dimension::Id::NormalZ);
//    //std::cout << "Normals dims: Nx:" << hasNx << " Ny:" << hasNy << " Nz:" << hasNz << "\n";
//
//
//    //pdal::PointViewPtr outs = *viewsNorm.begin();
//   
//    //////
//
//    //// (b) vypíš prvých pár bodov: XYZ + normály a normu
//    //std::size_t toPrint = std::min<std::size_t>(10, outs->size());
//    //std::size_t nonZero = 0;
//    //for (pdal::PointId i = 0; i < toPrint; ++i)
//    //{
//    //    double x = outs->getFieldAs<double>(pdal::Dimension::Id::X, i);
//    //    double y = outs->getFieldAs<double>(pdal::Dimension::Id::Y, i);
//    //    double z = outs->getFieldAs<double>(pdal::Dimension::Id::Z, i);
//
//    //    float nx = outs->getFieldAs<float>(pdal::Dimension::Id::NormalX, i);
//    //    float ny = outs->getFieldAs<float>(pdal::Dimension::Id::NormalY, i);
//    //    float nz = outs->getFieldAs<float>(pdal::Dimension::Id::NormalZ, i);
//
//    //    double len = std::sqrt(nx * nx + ny * ny + nz * nz);
//    //    if (len > 0.0) ++nonZero;
//
//    //    std::cout << "i=" << i
//    //        << "  XYZ=(" << x << "," << y << "," << z << ")"
//    //        << "  N=(" << nx << "," << ny << "," << nz << ")"
//    //        << "  |N|=" << len << "\n";
//    //}
//
//    //// (c) spoèítaj poèet nenulových normál v celom view
//    //std::size_t nzCount = 0;
//    ////////
//
//
//    //for (pdal::PointId i = 0; i < out->size(); ++i)
//    //{
//    //    float nx = out->getFieldAs<float>(pdal::Dimension::Id::NormalX, i);
//    //    float ny = out->getFieldAs<float>(pdal::Dimension::Id::NormalY, i);
//    //    float nz = out->getFieldAs<float>(pdal::Dimension::Id::NormalZ, i);
//    //    float x = out->getFieldAs<float>(pdal::Dimension::Id::X, i);
//    //    float y = out->getFieldAs<float>(pdal::Dimension::Id::Y, i);
//    //    float z = out->getFieldAs<float>(pdal::Dimension::Id::Z, i);
//
//    //    if (nx != 0.f || ny != 0.f || nz != 0.f) ++nzCount;
//
//
//    //    if ((nx !=0) || (ny != 0) || (nz != 0)){
//    //        std::cout << "id  " << i << std::endl;
//    //        std::cout << nx << std::endl;
//    //        std::cout << ny << std::endl;
//    //        std::cout << nz << std::endl << std::endl;
//    //    }
//    //    if (i<10) {
//    //        std::cout << "id  " << i << std::endl;
//    //        std::cout << x << std::endl;
//    //        std::cout << y << std::endl;
//    //        std::cout << z << std::endl;
//    //        std::cout << nx << std::endl;
//    //        std::cout << ny << std::endl;
//    //        std::cout << nz << std::endl << std::endl;
//    //    }
//    //    // ... urob èo potrebuješ (napr. uloži, spracova, vizualizova)
//    //}
//
//    //std::cout << "Nenulové normály: " << nzCount << " / " << out->size() << "\n";
//
//
//   /* std::vector<BuildingPoint> points;
//    auto layout = table.layout();
//    const bool hasIntensity =  layout->findDim("Intensity") != pdal::Dimension::Id::Unknown;
//    const bool hasClass = layout->findDim("Classification") != pdal::Dimension::Id::Unknown;*/
//
//
//        //filtering class 6
//    /*int sort_clasification = 6;
//    std::size_t kept = 0;
//    for (auto v : views) {
//        for (pdal::PointId i = 0; i < v->size(); ++i) {
//            int cls = hasClass ? v->getFieldAs<int>(pdal::Dimension::Id::Classification, i) : 0;
//            if (cls != sort_clasification) continue;
//
//            BuildingPoint p;
//            p.x = v->getFieldAs<double>(pdal::Dimension::Id::X, i);
//            p.y = v->getFieldAs<double>(pdal::Dimension::Id::Y, i);
//            p.z = v->getFieldAs<double>(pdal::Dimension::Id::Z, i);
//            p.intensity = hasIntensity ? v->getFieldAs<int>(pdal::Dimension::Id::Intensity, i) : 0;
//            p.classification = cls;
//
//            points.push_back(p);
//            ++kept;
//        }
//    }*/
//
//        //printing info about filtered data (class 6) 
//    /*std::cout << "Kept " << kept << " points with Classification== " << sort_clasification << "\n";
//    if (!points.empty()) {
//        const auto& p = points.front();
//        std::cout << "First kept: " << p.x << "," << p.y << "," << p.z
//            << " I=" << p.intensity << " cls=" << p.classification << "\n";
//    }*/
//
//
//
//        //loading all points from .las
//    {
//        ///////////
//
//    //pdal::StageFactory factory;
//
//    //auto reader = factory.createStage("readers.las");
//    //pdal::Options opts;
//    //opts.add("filename", "Mracno_bodov_vyrez.las");   // put your LAS path here
//    //reader->setOptions(opts);
//
//
//
//
//    //pdal::PointTable table;
//    //reader->prepare(table);     //metadata    
//    //pdal::PointViewSet views = reader->execute(table);  //reading data
//
//
//    ////suma
//    //std::size_t total = 0;
//    //for (auto v : views) total += v->size();
//    //std::cout << "Total points: " << total << "\n";
//
//
//    /////////////////
//    }
//
//
//
//
//    //auto v = *views.begin();                 // first view
//    //printSchema(table);
//    //dumpFirstN(table, *v, 10);
//
//        
//}


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

auto fail = [](const std::string& msg) { std::cerr << "[ERR] " << msg << "\n"; std::exit(1); };

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


    pdal::StageFactory factory;     ////stage


    //reader
//auto reader = factory.createStage("readers.las");
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
