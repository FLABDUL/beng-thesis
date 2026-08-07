cmake_minimum_required(VERSION 3.16)

project(beng_thesis_geometry_filter VERSION 1.1.0 LANGUAGES CXX)

option(BENG_THESIS_ENABLE_OPENMP "Enable parallel medial-axis processing" ON)
option(BENG_THESIS_VERBOSE "Print processing timings" OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(Eigen3 REQUIRED)
find_package(PCL 1.8 REQUIRED COMPONENTS common search features)
find_package(ZLIB REQUIRED)

if(BENG_THESIS_ENABLE_OPENMP)
  find_package(OpenMP QUIET)
endif()

add_library(cnpy STATIC tmp/cnpy/cnpy.cpp)
target_include_directories(cnpy PUBLIC tmp/cnpy ${ZLIB_INCLUDE_DIRS})
target_link_libraries(cnpy PUBLIC ZLIB::ZLIB)
target_compile_features(cnpy PUBLIC cxx_std_11)

add_library(masbcpp STATIC
  src/io.cpp
  src/compute_normals_processing.cpp
  src/compute_ma_processing.cpp
  src/simplify_processing.cpp
)
target_compile_features(masbcpp PUBLIC cxx_std_17)
target_include_directories(masbcpp
  PUBLIC
    src
    thirdparty
    ${PCL_INCLUDE_DIRS}
    ${EIGEN3_INCLUDE_DIR}
)
target_link_directories(masbcpp PUBLIC ${PCL_LIBRARY_DIRS})
target_link_libraries(masbcpp
  PUBLIC
    cnpy
    ${PCL_COMMON_LIBRARIES}
    ${PCL_SEARCH_LIBRARIES}
    ${PCL_FEATURES_LIBRARIES}
)
target_compile_definitions(masbcpp PUBLIC PCL_NO_PRECOMPILE)

foreach(definition IN LISTS PCL_DEFINITIONS)
  string(REGEX REPLACE "^-D" "" definition "${definition}")
  target_compile_definitions(masbcpp PUBLIC "${definition}")
endforeach()

if(OpenMP_CXX_FOUND)
  target_link_libraries(masbcpp PUBLIC OpenMP::OpenMP_CXX)
  target_compile_definitions(masbcpp PUBLIC WITH_OPENMP)
endif()

if(BENG_THESIS_VERBOSE)
  target_compile_definitions(masbcpp PUBLIC VERBOSEPRINT)
endif()

if(MSVC)
  target_compile_options(masbcpp PRIVATE /W4 /permissive-)
else()
  target_compile_options(masbcpp PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(compute_ma src/compute_ma.cpp)
add_executable(compute_normals src/compute_normals.cpp)
add_executable(simplify src/simplify.cpp)

foreach(executable compute_ma compute_normals simplify)
  target_link_libraries(${executable} PRIVATE masbcpp)
endforeach()

include(CTest)
if(BUILD_TESTING)
  add_test(NAME compute_ma_help COMMAND compute_ma --help)
  add_test(NAME compute_normals_help COMMAND compute_normals --help)
  add_test(NAME simplify_help COMMAND simplify --help)
endif()

install(TARGETS compute_ma compute_normals simplify RUNTIME DESTINATION bin)
