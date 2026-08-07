cmake_minimum_required(VERSION 3.14)
project(DecisionMesh LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(DMESH_NATIVE "Enable architecture-specific release tuning" ON)
option(DMESH_LTO "Enable interprocedural optimization in Release builds" ON)

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
  $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang>>:-O3>
  $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wpedantic>
  $<$<CXX_COMPILER_ID:MSVC>:/O2;/W4>)
if(DMESH_NATIVE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
  target_compile_options(decision_mesh PRIVATE $<$<CONFIG:Release>:-march=native;-mtune=native>)
endif()
if(DMESH_LTO)
  include(CheckIPOSupported)
  check_ipo_supported(RESULT DMESH_IPO_OK OUTPUT DMESH_IPO_ERROR)
  if(DMESH_IPO_OK)
    set_property(TARGET decision_mesh PROPERTY INTERPROCEDURAL_OPTIMIZATION_RELEASE TRUE)
  else()
    message(WARNING "IPO/LTO unavailable: ${DMESH_IPO_ERROR}")
  endif()
endif()

include(CTest)
if(BUILD_TESTING)
  add_executable(test_lindsey_guard tests/test_lindsey_guard.cpp lfdr.cpp)
  target_include_directories(test_lindsey_guard PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_options(test_lindsey_guard PRIVATE
    $<$<AND:$<CONFIG:Release>,$<CXX_COMPILER_ID:GNU,Clang>>:-O3>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wpedantic>
    $<$<CXX_COMPILER_ID:MSVC>:/O2;/W4>)
  add_test(NAME lindsey_guard COMMAND test_lindsey_guard)
endif()
