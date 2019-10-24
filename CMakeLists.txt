cmake_minimum_required(VERSION 2.8.9)


project (MCnetwork)
include_directories("./src/")

set(CMAKE_CXX_FLAGS "-std=c++11 -O0 -DDEBUG -g")
# set(CMAKE_CXX_FLAGS "-std=c++11")
# set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -DDEBUG" ${CMAKE_CXX_FLAGS_DEBUG})


find_package(HDF5 REQUIRED)
FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)

include_directories(${HDF5_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
set(_hdf5_libs hdf5 hdf5_cpp)

file(GLOB MCnetwork_SRC
  src/*.cpp
  src/*/*.cpp
)
add_executable(MCnetwork ${MCnetwork_SRC} ${sources})

add_library(enhance SHARED lib/enhance.cpp)

target_link_libraries(MCnetwork enhance)
target_link_libraries(MCnetwork ${_hdf5_libs})
target_link_libraries(MCnetwork ${Boost_LIBRARIES})

install(TARGETS MCnetwork RUNTIME DESTINATION bin)
