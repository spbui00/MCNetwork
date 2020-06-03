cmake_minimum_required(VERSION 3.1)


project (MCnetwork LANGUAGES CXX)
include_directories("./src/")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)
 
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

find_package(Boost COMPONENTS program_options REQUIRED)

find_package (HDF5 COMPONENTS CXX REQUIRED)
find_package (MFEM REQUIRED)

set(ENHANCE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/lib/enhance)
set(FINITE_ELEMENTE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/lib/finiteElemente)

add_library(enhance STATIC lib/enhance/enhance.cpp)
add_library(finiteElemente STATIC lib/finiteElemente/finiteElemente.cpp)

target_include_directories(finiteElemente PUBLIC ${MFEM_INCLUDE_DIRS})
target_link_libraries(finiteElemente mfem)

file(GLOB MCnetwork_SRC
  src/*.cpp
  src/*/*.cpp
)

add_executable(MCnetwork ${MCnetwork_SRC})
target_include_directories(MCnetwork PRIVATE ${HDF5_INCLUDE_DIRS} ${ENHANCE_INCLUDE_DIRECTORIES} ${FINITE_ELEMENTE_INCLUDE_DIRECTORIES})
target_link_libraries(MCnetwork Threads::Threads)
target_link_libraries(MCnetwork enhance)
target_link_libraries(MCnetwork finiteElemente)
target_link_libraries(MCnetwork ${HDF5_CXX_LIBRARIES})
target_link_libraries(MCnetwork ${Boost_LIBRARIES})

if(${CMAKE_VERSION} VERSION_LESS "3.8.0")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
else()
    set_property(TARGET MCnetwork PROPERTY CXX_STANDARD 17)
endif()

