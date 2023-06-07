#include <stdexcept>
#include <iostream>

#include "andres/graph/multicut-cubic/persistency.hxx"

inline void test(bool const condition) {
    if(!condition) throw std::logic_error("test failed");
}

typedef andres::graph::multicut_cubic::Problem<double> Problem;
typedef andres::graph::multicut_cubic::Persistency<double> Persistency;


void testTripletSubgraphCriterionUnfulfilledTriplet(){
    Problem problem(3);

    problem.costOfTriple({0, 1, 2}) = -3.1;
    problem.costOfEdge({0, 1}) = -1;
    problem.costOfEdge({0, 2}) = -1;
    problem.costOfEdge({1, 2}) = -1;

    Persistency  persistency(problem);
    persistency.edgeSubgraphCriterion();

    std::set<std::array<size_t , 2>> edges = persistency.edges1();

    test(persistency.edges1().empty());
}

void testTripletSubgraphCriterionFulfilledTripletSingleMerge(){
    Problem problem(3);

    problem.costOfTriple({0, 1, 2}) = -1.9;
    problem.costOfEdge({0, 1}) = -2;
    problem.costOfEdge({0, 2}) = 1;
    problem.costOfEdge({1, 2}) = 1;

    Persistency  persistency(problem);
    persistency.edgeSubgraphCriterion();

    std::set<std::array<size_t , 2>> edges = persistency.edges1();

    test(persistency.edges1().size() == 1);


    // here, all of them should be contained
    std::array<std::array<size_t, 2>, 1> expectedEdges{};
    expectedEdges[0] = {1, 0};

    for (auto & edge : expectedEdges){
        test(edges.find(edge) != edges.end());
    }
}

void testEdgeSubgraphCriterionFulfilledTripletDoubleMerge(){
    Problem problem(3);

    problem.costOfTriple({0, 1, 2}) = -2.1;
    problem.costOfEdge({0, 1}) = -2.2;
    problem.costOfEdge({0, 2}) = 1;
    problem.costOfEdge({1, 2}) = 1;

    Persistency  persistency(problem);
    persistency.edgeSubgraphCriterion();

    std::set<std::array<size_t , 2>> edges = persistency.edges1();

    test(persistency.edges1().size() == 3);


    // here, all of them should be contained
    std::array<std::array<size_t, 2>, 3> expectedEdges{};
    expectedEdges[0] = {1, 0};
    expectedEdges[1] = {2, 0};
    expectedEdges[2] = {2, 1};

    for (auto & edge : expectedEdges){
        test(edges.find(edge) != edges.end());
    }
}

void testEdgeSubgraphCriterionFulfilledDoublet(){
    Problem problem(2);

    problem.costOfEdge({0, 1}) = -1;

    Persistency  persistency(problem);
    persistency.edgeSubgraphCriterion();

    std::set<std::array<size_t , 2>> edges = persistency.edges1();

    test(persistency.edges1().size() == 1);


    // here, all of them should be contained
    std::array<std::array<size_t, 2>, 1> expectedEdges{};
    expectedEdges[0] = {1, 0};

    for (auto & edge : expectedEdges){
        test(edges.find(edge) != edges.end());
    }
}

void testEdgeSubgraphCriterionUnfullfilledDoublet(){
    Problem problem(2);

    problem.costOfEdge({0, 1}) = 1;

    Persistency  persistency(problem);
    persistency.edgeSubgraphCriterion();

    std::set<std::array<size_t , 2>> edges = persistency.edges1();

    test(persistency.edges1().size() == 0);

}

void testTripletSubgraphCriterionFulfilledArbitraryNumberOfNodes(int numberOfNodes = 5){
    Problem problem(numberOfNodes);

    for (size_t i = 0; i < numberOfNodes; ++i){
        for (size_t j = 0; j < i; ++j){
            problem.costOfEdge({i, j}) = -1;

            for (size_t k = 0; k < j; ++k){
                problem.costOfTriple({i, j, k}) = -1;
            }
        }
    }

    problem.costOfEdge({1, 0}) = -numberOfNodes*(numberOfNodes - 2) - 0.1;

    Persistency  persistency(problem);
    persistency.edgeSubgraphCriterion();

    std::set<std::array<size_t , 2>> edges = persistency.edges1();

    test(persistency.edges1().size() == 1);


    // here, all of them should be contained
    std::array<std::array<size_t, 2>, 1> expectedEdges{};
    expectedEdges[0] = {1, 0};

    for (auto & edge : expectedEdges){
        test(edges.find(edge) != edges.end());
    }
}

void testEdgeSubgraphCriterionUnfulfilledArbitraryNumberOfNodesNodes(int numberOfNodes = 5){
    Problem problem(numberOfNodes);

    for (size_t i = 0; i < numberOfNodes; ++i){
        for (size_t j = 0; j < i; ++j){
            problem.costOfEdge({i, j}) = -1;

            for (size_t k = 0; k < j; ++k){
                problem.costOfTriple({i, j, k}) = -1;
            }
        }
    }

    problem.costOfEdge({1, 0}) = -numberOfNodes*(numberOfNodes - 2) + 0.1;

    Persistency  persistency(problem);
    persistency.edgeSubgraphCriterion();

    std::set<std::array<size_t , 2>> edges = persistency.edges1();

    test(persistency.edges1().empty());
}

void testEdgeSubgraphCriterionFulfilledArbitraryNumberOfNodesNodesWithTripletBias(int numberOfNodes = 5){
    Problem problem(numberOfNodes);

    for (size_t i = 0; i < numberOfNodes; ++i){
        for (size_t j = 0; j < i; ++j){
            problem.costOfEdge({i, j}) = -1;

            for (size_t k = 0; k < j; ++k){
                problem.costOfTriple({i, j, k}) = -1;
            }
        }
    }

    problem.costOfEdge({1, 0}) = -numberOfNodes*(numberOfNodes - 2) -0.1 + 2;
    problem.costOfTriple({0, 2, 3}) = 0.1;
    problem.costOfTriple({1, 2, 3}) = 0.1;

    Persistency  persistency(problem);
    persistency.edgeSubgraphCriterion();

    std::set<std::array<size_t , 2>> edges = persistency.edges1();

    test(persistency.edges1().size() == 1);

    // here, all of them should be contained
    std::array<std::array<size_t, 2>, 1> expectedEdges{};
    expectedEdges[0] = {1, 0};

    for (auto & edge : expectedEdges){
        test(edges.find(edge) != edges.end());
    }
}

void testEdgeSubgraphCriterionUnfulfilledArbitraryNumberOfNodesNodesWithTripletBias(int numberOfNodes = 5){
    Problem problem(numberOfNodes);

    for (size_t i = 0; i < numberOfNodes; ++i){
        for (size_t j = 0; j < i; ++j){
            problem.costOfEdge({i, j}) = -1;

            for (size_t k = 0; k < j; ++k){
                problem.costOfTriple({i, j, k}) = -1;
            }
        }
    }

    problem.costOfEdge({1, 0}) = -numberOfNodes*(numberOfNodes - 2) +0.1 + 2;
    problem.costOfTriple({0, 2, 3}) = 0.1;
    problem.costOfTriple({1, 2, 3}) = 0.1;

    Persistency  persistency(problem);
    persistency.edgeSubgraphCriterion();

    std::set<std::array<size_t , 2>> edges = persistency.edges1();

    test(persistency.edges1().empty());

}


int main() {
    testEdgeSubgraphCriterionUnfullfilledDoublet();
    testEdgeSubgraphCriterionFulfilledDoublet();

    testTripletSubgraphCriterionUnfulfilledTriplet();
    testTripletSubgraphCriterionFulfilledTripletSingleMerge();
    testEdgeSubgraphCriterionFulfilledTripletDoubleMerge();

    testTripletSubgraphCriterionFulfilledArbitraryNumberOfNodes(5);
    testTripletSubgraphCriterionFulfilledArbitraryNumberOfNodes(10);
    testTripletSubgraphCriterionFulfilledArbitraryNumberOfNodes(20);
    testEdgeSubgraphCriterionUnfulfilledArbitraryNumberOfNodesNodes(5);
    testEdgeSubgraphCriterionUnfulfilledArbitraryNumberOfNodesNodes(10);
    testEdgeSubgraphCriterionUnfulfilledArbitraryNumberOfNodesNodes(20);

    testEdgeSubgraphCriterionFulfilledArbitraryNumberOfNodesNodesWithTripletBias(5);
    testEdgeSubgraphCriterionFulfilledArbitraryNumberOfNodesNodesWithTripletBias(10);
    testEdgeSubgraphCriterionFulfilledArbitraryNumberOfNodesNodesWithTripletBias(20);

    testEdgeSubgraphCriterionUnfulfilledArbitraryNumberOfNodesNodesWithTripletBias(5);
    testEdgeSubgraphCriterionUnfulfilledArbitraryNumberOfNodesNodesWithTripletBias(10);
    testEdgeSubgraphCriterionUnfulfilledArbitraryNumberOfNodesNodesWithTripletBias(20);
    return 0;
}
