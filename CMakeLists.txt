cmake_minimum_required(VERSION 3.5)
project(MCnetwork LANGUAGES CXX)

include_directories("./src/")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(BOOST_ROOT "/opt/homebrew/Cellar/boost/1.84.0_1")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

add_compile_definitions(BOOST_ALLOW_DEPRECATED_HEADERS)

find_package(Threads REQUIRED)

find_package(Boost 1.67 REQUIRED COMPONENTS program_options)
if(NOT Boost_FOUND)
  message(FATAL_ERROR "Boost not found")
endif()

find_package(HDF5 COMPONENTS CXX REQUIRED)
if(NOT HDF5_FOUND)
  message(FATAL_ERROR "HDF5 not found")
endif()

find_package(MFEM REQUIRED)
if(NOT MFEM_FOUND)
  message(FATAL_ERROR "MFEM not found")
endif()

file(GLOB MCnetwork_SRC
  src/*.cpp
  src/*/*.cpp
)

if(SWAPTRACKER)
  add_definitions(-DSWAPTRACKER)
endif()

add_subdirectory(lib/)

add_executable(MCnetwork ${MCnetwork_SRC})
target_include_directories(MCnetwork PRIVATE
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/src/system
    ${PROJECT_SOURCE_DIR}/lib/enhance
    ${PROJECT_SOURCE_DIR}/lib/finiteElemente
    ${Boost_INCLUDE_DIRS}
    ${HDF5_INCLUDE_DIRS}
    ${ENHANCE_INCLUDE_DIRECTORIES}
    ${FINITEELEMENTE_INCLUDE_DIRECTORIES}
)
target_link_libraries(MCnetwork 
    Threads::Threads 
    ${Boost_LIBRARIES}
    ${HDF5_CXX_LIBRARIES} 
    enhance 
    finiteElemente
)
