cmake_minimum_required(VERSION 3.10)
project(pvd)

# Set the C++ standard to at least C++11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# To change to OPTIMIZE=OFF use: cmake -DOPTIMIZE=OFF ..
option(OPTIMIZE "Build with optimization [default: on]" ON)

if(OPTIMIZE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
else()
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type" FORCE)
endif()

# Put the executables in a bin directory
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_CURRENT_SOURCE_DIR}/bin)

set(SOURCES
  src/main.cpp
  src/variant_site.cpp
  src/pg_utils.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

# add zlib
find_package(ZLIB REQUIRED)
target_link_libraries(${PROJECT_NAME} ZLIB::ZLIB)
