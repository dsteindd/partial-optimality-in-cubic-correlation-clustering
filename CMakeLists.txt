cmake_minimum_required(VERSION 3.1)
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)
#set(CMAKE_BUILD_TYPE DEBUG)

project(cubic-multicut)
enable_testing()
#set(CMAKE_CXX_STANDARD 20)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

# tests namespace andres
add_executable(test-multicut-cubic-problem src/test/andres/graph/multicut-cubic/problem.cxx)
target_include_directories(test-multicut-cubic-problem PUBLIC include)
add_test(test-multicut-cubic-problem test-multicut-cubic-problem)

add_executable(test-partition src/test/partition/partition.cxx)
target_include_directories(test-partition PUBLIC include)
add_test(test-partition test-partition)

# tests namespace dssc
add_executable(test-math src/test/dssc/math.cxx)
target_include_directories(test-math PUBLIC include)
add_test(test-math test-math)

add_executable(test-triangles-adaptor-costs src/test/dssc/adaptors/triangles/costs.cxx)
target_include_directories(test-triangles-adaptor-costs PUBLIC include)
add_test(test-triangles-adaptor-costs test-triangles-adaptor-costs)

add_executable(test-qpbo-adaptor src/test/andres/graph/multicut-cubic/qpbo-adaptor.cxx)
target_include_directories(test-qpbo-adaptor PUBLIC include)
add_test(test-qpbo-adaptor test-qpbo-adaptor)

add_executable(test-qpbo-problem src/test/andres/graph/qpbo-problem.cxx)
target_include_directories(test-qpbo-problem PUBLIC include)
add_test(test-qpbo-problem test-qpbo-problem)

add_executable(test-submodular-qpbo-solver src/test/andres/graph/submodular-qpbo-maxflow-solver.cxx)
target_include_directories(test-submodular-qpbo-solver PUBLIC include)
add_test(test-submodular-qpbo-solver test-submodular-qpbo-solver)

add_executable(test-edge-cut src/test/andres/graph/multicut-cubic/persistency-edge-cut.cxx)
target_include_directories(test-edge-cut PUBLIC include)
add_test(test-edge-cut test-edge-cut)

add_executable(test-edge-join src/test/andres/graph/multicut-cubic/persistency-edge-join.cxx)
target_include_directories(test-edge-join PUBLIC include)
add_test(test-edge-join test-edge-join)

add_executable(test-triplet-cut src/test/andres/graph/multicut-cubic/persistency-triplet-cut.cxx)
target_include_directories(test-triplet-cut PUBLIC include)
add_test(test-triplet-cut test-triplet-cut)

add_executable(test-independent-subproblems src/test/andres/graph/multicut-cubic/persistency-independent-subproblems.cxx)
target_include_directories(test-independent-subproblems PUBLIC include)
add_test(test-independent-subproblems test-independent-subproblems)

add_executable(test-triplet-join src/test/andres/graph/multicut-cubic/persistency-triplet-join.cxx)
target_include_directories(test-triplet-join PUBLIC include)
add_test(test-triplet-join test-triplet-join)

add_executable(test-edge-subgraph src/test/andres/graph/multicut-cubic/persistency-edge-subgraph.cxx)
target_include_directories(test-edge-subgraph PUBLIC include)
add_test(test-edge-subgraph test-edge-subgraph)

add_executable(test-triplet-subgraph src/test/andres/graph/multicut-cubic/persistency-triplet-subgraph.cxx)
target_include_directories(test-triplet-subgraph PUBLIC include)
add_test(test-triplet-subgraph test-triplet-subgraph)

add_executable(test-triplet-edge-join src/test/andres/graph/multicut-cubic/persistency-triplet-edge-join.cxx)
target_include_directories(test-triplet-edge-join PUBLIC include)
add_test(test-triplet-edge-join test-triplet-edge-join)

add_executable(test-cost-rescaling src/test/dssc/adaptors/problem-rescaling-adaptor.cxx)
target_include_directories(test-cost-rescaling PUBLIC include)
add_test(test-cost-rescaling test-cost-rescaling)

add_executable(test-subset-join src/test/andres/graph/multicut-cubic/persistency-subset-join.cxx)
target_include_directories(test-subset-join PUBLIC include)
add_test(test-subset-join test-subset-join)

# experiments
add_executable(equilateral-triangles-persistency src/equilateral-triangles/persistency.cxx)
target_include_directories(equilateral-triangles-persistency PUBLIC include)

add_executable(predefined-partition-persistency src/predefined-partition/persistency.cxx)
target_include_directories(predefined-partition-persistency PUBLIC include)
