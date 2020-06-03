set(ENHANCE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/enhance PARENT_SCOPE)
set(FINITE_ELEMENTE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_SOURCE_DIR}/finiteElemente PARENT_SCOPE)

add_library(enhance STATIC enhance/enhance.cpp)
add_library(finiteElemente STATIC finiteElemente/finiteElemente.cpp)

target_include_directories(finiteElemente PUBLIC ${MFEM_INCLUDE_DIRS})
target_link_libraries(finiteElemente mfem)

