cmake_minimum_required(VERSION 3.14)
project(DecisionMesh LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(decision_mesh
  gated_main.cpp
  pipeline.cpp
  output.cpp
  run_config.cpp
  dataset.cpp
  diagnostics.cpp
  hierarchical_design.cpp
  lfdr.cpp
  mesh.cpp
  vertex.cpp
  edge.cpp
  face.cpp
  tree.cpp)

target_compile_options(decision_mesh PRIVATE
  $<$<CXX_COMPILER_ID:GNU,Clang>:-O2;-Wall;-Wextra;-Wpedantic>
  $<$<CXX_COMPILER_ID:MSVC>:/O2;/W4>)

include(CTest)
if(BUILD_TESTING)
  add_executable(test_lindsey_guard tests/test_lindsey_guard.cpp lfdr.cpp)
  target_include_directories(test_lindsey_guard PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_options(test_lindsey_guard PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-O2;-Wall;-Wextra;-Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/O2;/W4>)
  add_test(NAME lindsey_guard COMMAND test_lindsey_guard)
endif()
