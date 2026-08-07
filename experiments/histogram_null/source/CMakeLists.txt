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
