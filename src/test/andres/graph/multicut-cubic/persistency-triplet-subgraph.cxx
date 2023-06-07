#include <stdexcept>
#include <iostream>

#include "andres/graph/multicut-cubic/persistency.hxx"

inline void test(bool const condition) {
    if(!condition) throw std::logic_error("test failed");
}

typedef andres::graph::multicut_cubic::Problem<double> Problem;
typedef andres::graph::multicut_cubic::Persistency<double> Persistency;


void testTripletSubgraphCriterionUnfulfilledTriplet(){
    {
        Problem problem(3);

        problem.costOfTriple({0, 1, 2}) = -3.1;
        problem.costOfEdge({0, 1}) = +1;
        problem.costOfEdge({0, 2}) = +1;
        problem.costOfEdge({1, 2}) = +1;

        Persistency persistency(problem);
        persistency.tripletSubgraphCriterion();

        auto edges = persistency.edges1();

        test(persistency.edges1().empty());
    }

    {
        Problem problem(3);

        problem.costOfTriple({0, 1, 2}) = 8;
        problem.costOfEdge({0, 1}) = -1;
        problem.costOfEdge({0, 2}) = -1;
        problem.costOfEdge({1, 2}) = -1;

        Persistency persistency(problem);
        persistency.tripletSubgraphCriterion();

        auto edges = persistency.edges1();

        test(persistency.edges1().empty());
    }

    {
        Problem problem(3);

        problem.costOfTriple({0, 1, 2}) = 4;
        problem.costOfEdge({0, 1}) = -1;
        problem.costOfEdge({0, 2}) = -1;
        problem.costOfEdge({1, 2}) = -1;

        Persistency persistency(problem);
        persistency.tripletSubgraphCriterion();

        auto edges = persistency.edges1();

        test(persistency.edges1().empty());
    }

    {
        Problem problem(3);

        problem.costOfTriple({0, 1, 2}) = 9.1;
        problem.costOfEdge({0, 1}) = -1;
        problem.costOfEdge({0, 2}) = -1;
        problem.costOfEdge({1, 2}) = -8;

        Persistency persistency(problem);
        persistency.tripletSubgraphCriterion();

        auto edges = persistency.edges1();

        test(persistency.edges1().empty());
    }
}

void testTripletSubgraphCriterionFulfilledTripletSingleMerge(){
    {
        Problem problem(3);

        problem.costOfTriple({0, 1, 2}) = 2.1;
        problem.costOfEdge({0, 1}) = -1;
        problem.costOfEdge({0, 2}) = -1;
        problem.costOfEdge({1, 2}) = -8;

        Persistency persistency(problem);
        persistency.tripletSubgraphCriterion();

        std::set<std::array<size_t, 2>> edges = persistency.edges1();

        test(persistency.edges1().size() == 1);

        test(edges.find({2, 1}) != edges.end());
    }

    {
        Problem problem(3);

        problem.costOfTriple({0, 1, 2}) = 8.9;
        problem.costOfEdge({0, 1}) = -1;
        problem.costOfEdge({0, 2}) = -1;
        problem.costOfEdge({1, 2}) = -8;

        Persistency persistency(problem);
        persistency.tripletSubgraphCriterion();

        std::set<std::array<size_t, 2>> edges = persistency.edges1();

        test(persistency.edges1().size() == 1);

        test(edges.find({2, 1}) != edges.end());
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
    {
        Problem problem(numberOfNodes);

        for (size_t i = 0; i < numberOfNodes; ++i) {
            for (size_t j = 0; j < i; ++j) {
                problem.costOfEdge({i, j}) = -1;

                for (size_t k = 0; k < j; ++k) {
                    problem.costOfTriple({i, j, k}) = -1;
                }
            }
        }

        problem.costOfEdge({1, 0}) = -1;
        problem.costOfEdge({2, 0}) = -2;
        problem.costOfEdge({2, 1}) = -3;
        problem.costOfTriple({2, 1, 0}) = -numberOfNodes * (numberOfNodes - 3) * 3 / 2 + 3.1;

        Persistency persistency(problem);
        persistency.tripletSubgraphCriterion();

        auto edges = persistency.edges1();

        test(edges.find({2, 1}) != edges.end());
    }

    {
        Problem problem(numberOfNodes);

        for (size_t i = 0; i < numberOfNodes; ++i) {
            for (size_t j = 0; j < i; ++j) {
                problem.costOfEdge({i, j}) = -1;

                for (size_t k = 0; k < j; ++k) {
                    problem.costOfTriple({i, j, k}) = -1;
                }
            }
        }

        problem.costOfEdge({1, 0}) = -1;
        problem.costOfEdge({2, 0}) = -2;
        problem.costOfEdge({2, 1}) = -3;
        problem.costOfTriple({2, 1, 0}) = -numberOfNodes * (numberOfNodes - 3) * 3 / 2 + 3.9;

        Persistency persistency(problem);
        persistency.tripletSubgraphCriterion();

        auto edges = persistency.edges1();

        test(edges.find({2, 1}) != edges.end());
    }

    {
        Problem problem(numberOfNodes);

        for (size_t i = 0; i < numberOfNodes; ++i) {
            for (size_t j = 0; j < i; ++j) {
                problem.costOfEdge({i, j}) = -1;

                for (size_t k = 0; k < j; ++k) {
                    problem.costOfTriple({i, j, k}) = -1;
                }
            }
        }

        problem.costOfEdge({1, 0}) = -1;
        problem.costOfEdge({2, 0}) = -1;
        problem.costOfEdge({2, 1}) = -numberOfNodes * (numberOfNodes - 3) * 3 / 2 + 1.9;
        problem.costOfTriple({2, 1, 0}) = -1;

        Persistency persistency(problem);
        persistency.tripletSubgraphCriterion();

        auto edges = persistency.edges1();

        test(edges.size() == 1);

        test(edges.find({2, 1}) != edges.end());
    }
}


void testTripletSubgraphCriterionFulfilledArbitraryNumberOfNodesPositiveTripletWeights(int numberOfNodes = 5){
    {
        Problem problem(numberOfNodes);

        for (size_t i = 0; i < numberOfNodes; ++i) {
            for (size_t j = 0; j < i; ++j) {
                problem.costOfEdge({i, j}) = -1;

                for (size_t k = 0; k < j; ++k) {
                    problem.costOfTriple({i, j, k}) = +1;
                }
            }
        }

        problem.costOfEdge({1, 0}) = -1;
        problem.costOfEdge({2, 0}) = -1;
        problem.costOfEdge({2, 1}) = -3*(numberOfNodes - 3) - 0.1;
        problem.costOfTriple({2, 1, 0}) = +1;

        Persistency persistency(problem);
        persistency.tripletSubgraphCriterion();

        auto edges = persistency.edges1();

        test(edges.size() == 1);

        test(edges.find({2, 1}) != edges.end());
    }
    {
        Problem problem(numberOfNodes);

        for (size_t i = 0; i < numberOfNodes; ++i) {
            for (size_t j = 0; j < i; ++j) {
                problem.costOfEdge({i, j}) = -1;

                for (size_t k = 0; k < j; ++k) {
                    problem.costOfTriple({i, j, k}) = +1;
                }
            }
        }

        problem.costOfEdge({1, 0}) = -1;
        problem.costOfEdge({2, 0}) = -1;
        problem.costOfEdge({2, 1}) = -3*(numberOfNodes - 3);
        problem.costOfTriple({2, 1, 0}) = +1;

        Persistency persistency(problem);
        persistency.tripletSubgraphCriterion();

        auto edges = persistency.edges1();

        test(edges.size() == 1);

        test(edges.find({2, 1}) != edges.end());
    }
}


void testTripletSubgraphCriterionUnfulfilledArbitraryNumberOfNodesPositiveTripletWeights(int numberOfNodes = 5){
    {
        Problem problem(numberOfNodes);

        for (size_t i = 0; i < numberOfNodes; ++i) {
            for (size_t j = 0; j < i; ++j) {
                problem.costOfEdge({i, j}) = -1;

                for (size_t k = 0; k < j; ++k) {
                    problem.costOfTriple({i, j, k}) = +1;
                }
            }
        }

        problem.costOfEdge({1, 0}) = -1;
        problem.costOfEdge({2, 0}) = -1;
        problem.costOfEdge({2, 1}) = -3*(numberOfNodes - 3) + 0.1;
        problem.costOfTriple({2, 1, 0}) = +1;

        Persistency persistency(problem);
        persistency.tripletSubgraphCriterion();

        auto edges = persistency.edges1();

        test(edges.empty() == 1);
    }
}

void testTripletSubgraphCriterionUnfulfilledArbitraryNumberOfNodes(int numberOfNodes = 5){
    {
        Problem problem(numberOfNodes);

        for (size_t i = 0; i < numberOfNodes; ++i) {
            for (size_t j = 0; j < i; ++j) {
                problem.costOfEdge({i, j}) = -1;

                for (size_t k = 0; k < j; ++k) {
                    problem.costOfTriple({i, j, k}) = -1;
                }
            }
        }

        problem.costOfEdge({1, 0}) = -1;
        problem.costOfEdge({2, 0}) = -2;
        problem.costOfEdge({2, 1}) = -3;
        problem.costOfTriple({2, 1, 0}) = -numberOfNodes * (numberOfNodes - 3) * 3 / 2 + 4.1;

        Persistency persistency(problem);
        persistency.tripletSubgraphCriterion();

        auto edges = persistency.edges1();

        test(edges.empty());
    }

    {
        Problem problem(numberOfNodes);

        for (size_t i = 0; i < numberOfNodes; ++i) {
            for (size_t j = 0; j < i; ++j) {
                problem.costOfEdge({i, j}) = -1;

                for (size_t k = 0; k < j; ++k) {
                    problem.costOfTriple({i, j, k}) = -1;
                }
            }
        }

        problem.costOfEdge({1, 0}) = -1;
        problem.costOfEdge({2, 0}) = -1;
        problem.costOfEdge({2, 1}) = -numberOfNodes * (numberOfNodes - 3) * 3 / 2 + 2.1;
        problem.costOfTriple({2, 1, 0}) = -1;

        Persistency persistency(problem);
        persistency.tripletSubgraphCriterion();

        auto edges = persistency.edges1();

        test(edges.empty());
    }
}



int main() {
    testTripletSubgraphCriterionUnfulfilledTriplet();
    testTripletSubgraphCriterionFulfilledTripletSingleMerge();

    testTripletSubgraphCriterionFulfilledArbitraryNumberOfNodes(5);
    testTripletSubgraphCriterionFulfilledArbitraryNumberOfNodes(10);
    testTripletSubgraphCriterionFulfilledArbitraryNumberOfNodes(20);

    testTripletSubgraphCriterionFulfilledArbitraryNumberOfNodesPositiveTripletWeights(5);
    testTripletSubgraphCriterionFulfilledArbitraryNumberOfNodesPositiveTripletWeights(10);
    testTripletSubgraphCriterionFulfilledArbitraryNumberOfNodesPositiveTripletWeights(20);

    testTripletSubgraphCriterionUnfulfilledArbitraryNumberOfNodesPositiveTripletWeights(5);
    testTripletSubgraphCriterionUnfulfilledArbitraryNumberOfNodesPositiveTripletWeights(10);
    testTripletSubgraphCriterionUnfulfilledArbitraryNumberOfNodesPositiveTripletWeights(20);

    testTripletSubgraphCriterionUnfulfilledArbitraryNumberOfNodes(5);
    testTripletSubgraphCriterionUnfulfilledArbitraryNumberOfNodes(10);
    testTripletSubgraphCriterionUnfulfilledArbitraryNumberOfNodes(20);
    return 0;
}
