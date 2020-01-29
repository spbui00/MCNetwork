cmake_minimum_required(VERSION 2.8.9)


project (MCnetwork)
include_directories("./src/")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -pthread -O0 -g -Wl,--start-group")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSWAPTRACKER -std=c++17 -pthread -O0 -g -Wl,--start-group")

find_package(Boost COMPONENTS program_options REQUIRED)
find_package (HDF5 REQUIRED COMPONENTS CXX)
set(_hdf5_libs hdf5 hdf5_cpp)

file(GLOB MCnetwork_SRC
  src/*.cpp
  src/*/*.cpp
)
add_executable(MCnetwork ${MCnetwork_SRC} ${sources})

add_library(enhance SHARED lib/enhance.cpp)
add_library(finiteElemente SHARED lib/finiteElemente/finiteElemente.cpp)

target_link_libraries(MCnetwork enhance)
target_link_libraries(MCnetwork finiteElemente)
target_link_libraries(MCnetwork ${_hdf5_libs})
target_link_libraries(MCnetwork ${Boost_LIBRARIES})
target_link_libraries(MCnetwork mfem)


# install(TARGETS MCnetwork RUNTIME DESTINATION bin)
