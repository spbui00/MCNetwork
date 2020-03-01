cmake_minimum_required(VERSION 2.8.9)


project (MCnetwork)
include_directories("./src/")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif(NOT CMAKE_BUILD_TYPE)
 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -pthread -O0 -g -Wl,--start-group")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSWAPTRACKER -std=c++17 -pthread -O0 -g -Wl,--start-group")

find_package(Boost COMPONENTS program_options REQUIRED)



file(GLOB MCnetwork_SRC
  src/*.cpp
  src/*/*.cpp
)

find_package (HDF5 COMPONENTS CXX)
find_package (MFEM)


if (NOT (${HDF5_FOUND} OR ${MFEM_FOUND}))
  message(" -->> could not find MFEM or HDF5(incl CXX) in default path, maybe on bagheera? searching locally...  <<-- ")

  find_package (MFEM REQUIRED PATHS /home/m_beck52/local)
  find_package (HDF5 REQUIRED COMPONENTS CXX PATHS /home/m_beck52/local NO_DEFAULT_PATH)

  message(" -->> found HDF5 and MFEM ! <<-- ")

  link_directories(/home/m_beck52/local/lib/)
  include_directories(/home/m_beck52/local/include/)
  file(GLOB BAGHEERA_SRC /home/m_beck52/local/include/*.hpp)


endif (NOT (${HDF5_FOUND} OR ${MFEM_FOUND}))

# include_directories(${HDF5_INCLUDE_DIRS})


add_executable(MCnetwork ${MCnetwork_SRC} ${sources} ${BAGHEERA_SRC})


add_library(enhance SHARED lib/enhance.cpp)
add_library(finiteElemente SHARED lib/finiteElemente/finiteElemente.cpp)


target_link_libraries(MCnetwork enhance)
target_link_libraries(MCnetwork finiteElemente)
target_link_libraries(MCnetwork hdf5 hdf5_cpp)
target_link_libraries(MCnetwork ${Boost_LIBRARIES})
target_link_libraries(MCnetwork mfem)
