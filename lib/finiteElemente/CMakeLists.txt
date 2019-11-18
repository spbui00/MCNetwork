cmake_minimum_required(VERSION 2.8.9)


project (finEle)
include_directories("./src/")

set(CMAKE_CXX_FLAGS "-std=c++11 -O0 -DDEBUG -g -Wl,--start-group")
# set(CMAKE_CXX_FLAGS "-std=c++11")
# set(CMAKE_CXX_FLAGS_DEBUG "-std=c++11 -DDEBUG" ${CMAKE_CXX_FLAGS_DEBUG})



add_executable(finEle main.cpp finiteElemente.cpp)


target_link_libraries (finEle mfem)


