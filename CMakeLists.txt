cmake_minimum_required(VERSION 3.20)
project(MyPdalApp LANGUAGES CXX)

# C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find PDAL (uses the PDALConfig.cmake provided by your install)
# You can pin a version if you want: find_package(PDAL 2.6 REQUIRED CONFIG)
find_package(PDAL REQUIRED CONFIG)

# If PDAL provides legacy vars, use them; if it provides an imported target, use that.
add_executable(MyPdalApp
    src/main.cpp
)

# Prefer modern imported target if it exists
if(TARGET PDAL::pdalcpp)
    target_link_libraries(MyPdalApp PRIVATE PDAL::pdalcpp)
else()
    # Back-compat: variables exported by PDAL’s package config
    target_link_libraries(MyPdalApp PRIVATE ${PDAL_LIBRARIES})
    target_include_directories(MyPdalApp PRIVATE ${PDAL_INCLUDE_DIRS})
    add_definitions(${PDAL_DEFINITIONS})
endif()
