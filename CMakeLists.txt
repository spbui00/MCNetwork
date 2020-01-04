cmake_minimum_required(VERSION 2.8.9)


project (MCnetwork)
include_directories("./src/")

set(CMAKE_CXX_FLAGS "-std=c++11 -O0 -DDEBUG -g -Wl,--start-group")
# set(CMAKE_CXX_FLAGS "-std=c++11")
# set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -DDEBUG" ${CMAKE_CXX_FLAGS_DEBUG})

find_package(Boost COMPONENTS program_options REQUIRED)

link_directories(/usr/local/hdf5/lib)

include_directories(/usr/local/hdf5/include ${Boost_INCLUDE_DIRS})


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
target_link_libraries (MCnetwork mfem)


# install(TARGETS MCnetwork RUNTIME DESTINATION bin)
