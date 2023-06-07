#include <stdexcept>
#include <iostream>

#include "andres/graph/multicut-cubic/persistency.hxx"

inline void test(bool const condition) {
    if(!condition) throw std::logic_error("test failed");
}

typedef andres::graph::multicut_cubic::Problem<double> Problem;
typedef andres::graph::multicut_cubic::Persistency<double> Persistency;


void testTripletEdgeJoinCriterionFulfilledSingleTriplet(){
    Problem problem(3);

    problem.costOfTriple({0, 1, 2}) = -2.1;

    Persistency  persistency(problem);
    persistency.tripletEdgeJoinCriterion();

    test(persistency.edges1().size() == 1);

    std::set<std::array<size_t, 2>> expectedEdges;

    // todo: add assert for correct edge
}



void testTripletEdgeJoinCriterionFulfilledOneEdge(){
    Problem problem(5);

    problem.costOfTriple({0, 1, 2}) = -5;
    problem.costOfTriple({0, 1, 3}) = 1;
    problem.costOfTriple({0, 1, 4}) = 1;
    problem.costOfTriple({0, 2, 3}) = 1;
    problem.costOfTriple({0, 2, 4}) = 1;
    problem.costOfTriple({1, 2, 3}) = 1;
    problem.costOfTriple({1, 2, 4}) = 1;
    problem.costOfTriple({2, 3, 4}) = 0.6;
    problem.costOfTriple({0, 3, 4}) = 0.6;
    problem.costOfTriple({1, 3, 4}) = 1.5;

    Persistency  persistency(problem);
    persistency.tripletEdgeJoinCriterion();

    test(persistency.edges1().size() == 1);

    std::set<std::array<size_t, 2>> expectedEdges;
    expectedEdges.insert({2, 0});

    for (auto const edge : expectedEdges){
        test(persistency.edges1().find(edge) != persistency.edges1().end());
    }
}





void testTripletEdgeJoinCriterionOneEdgeFulfilled(){

    Problem problem(4);

    problem.costOfTriple({0, 1, 2}) = -2.1;
    problem.costOfTriple({0, 1, 3}) = 1;
    problem.costOfTriple({0, 2, 3}) = 1;
    problem.costOfTriple({1, 2, 3}) = 1;

    Persistency  persistency(problem);
    persistency.tripletEdgeJoinCriterion();

    std::set<std::array<size_t , 2>> edges = persistency.edges1();

    test(persistency.edges1().size() == 1);

    std::array<std::array<size_t, 2>, 3> expectedEdges{};
    expectedEdges[0] = {1, 0};
    expectedEdges[1] = {2, 0};
    expectedEdges[2] = {2, 1};

    std::array<bool, 3> isEdgeContained{};

    for (size_t index = 0; index < expectedEdges.size(); index++){
        isEdgeContained[index] = (edges.find(expectedEdges[index]) != edges.end());
    }

    test(std::any_of(isEdgeContained.begin(), isEdgeContained.end(),[](bool isContained){return isContained;}));
}

void testTripletEdgeJoinCriterionTwoEdgesFulfilled(){

    Problem problem(4);

    problem.costOfTriple({0, 1, 2}) = -3.1;
    problem.costOfTriple({0, 1, 3}) = 1;
    problem.costOfTriple({0, 2, 3}) = 1;
    problem.costOfTriple({1, 2, 3}) = 1;

    Persistency  persistency(problem);
    persistency.tripletEdgeJoinCriterion();

    std::set<std::array<size_t , 2>> edges = persistency.edges1();

    test(edges.size() == 3);


    // here, all of them should be contained
    std::array<std::array<size_t, 2>, 3> expectedEdges{};
    expectedEdges[0] = {1, 0};
    expectedEdges[1] = {2, 0};
    expectedEdges[2] = {2, 1};

    for (auto & edge : expectedEdges){
        test(edges.find(edge) != edges.end());
    }
}


void testTripletEdgeJoinCriterionUnfulfilled(){
    Problem problem(4);

    problem.costOfTriple({0, 1, 2}) = -1.9;
    problem.costOfTriple({0, 1, 3}) = 1;
    problem.costOfTriple({0, 2, 3}) = 1;
    problem.costOfTriple({1, 2, 3}) = 1;

    Persistency  persistency(problem);
    persistency.tripletEdgeJoinCriterion();

    test(persistency.edges1().empty());
}

void testTripletEdgeJoinCriterionUnfulfilledFiveNodes(){
    int numberOfNodes = 5;

    Problem problem(numberOfNodes);

    for (size_t i = 0; i < numberOfNodes; ++i){
        for (size_t j = 0; j < i; ++j){
            problem.costOfEdge({i, j}) = 0;

            for (size_t k = 0; k < j; k++){
                problem.costOfTriple({i, j, k}) = -1;
            }
        }
    }

    problem.costOfTriple({2, 1, 0}) = -2.9;


    Persistency  persistency(problem);
    persistency.tripletEdgeJoinCriterion();

    test(persistency.edges1().empty());
}

void testTripletEdgeJoinCriterionFulfilledArbitraryNumberOfNodes(int numberOfNodes = 5){
    Problem problem(numberOfNodes);

    for (size_t i = 0; i < numberOfNodes; ++i){
        for (size_t j = 0; j < i; ++j){
            problem.costOfEdge({i, j}) = +1;

            for (size_t k = 0; k < j; k++){
                problem.costOfTriple({i, j, k}) = +1;
            }
        }
    }

    problem.costOfEdge({2, 1}) = -(numberOfNodes + 1)*(numberOfNodes -2)/2 - 0.1;


    Persistency  persistency(problem);
    persistency.tripletEdgeJoinCriterion();

    auto edges = persistency.edges1();

    test(edges.size() == 1);
    test(edges.find({2, 1}) != edges.end());
}

void testTripletEdgeJoinCriterionUnfulfilledArbitraryNumberOfNodes(int numberOfNodes = 5){
    Problem problem(numberOfNodes);

    for (size_t i = 0; i < numberOfNodes; ++i){
        for (size_t j = 0; j < i; ++j){
            problem.costOfEdge({i, j}) = +1;

            for (size_t k = 0; k < j; k++){
                problem.costOfTriple({i, j, k}) = +1;
            }
        }
    }

    problem.costOfEdge({2, 1}) = -(numberOfNodes + 1)*(numberOfNodes -2)/2 +0.1;


    Persistency  persistency(problem);
    persistency.tripletEdgeJoinCriterion();

    auto edges = persistency.edges1();

    test(edges.empty());
    test(edges.find({2, 1}) == edges.end());
}

void testTripletEdgeJoinCriterionUnfulfilledArbitraryNumberOfNodes2(int numberOfNodes = 5){
    Problem problem(numberOfNodes);

    for (size_t i = 0; i < numberOfNodes; ++i){
        for (size_t j = 0; j < i; ++j){
            problem.costOfEdge({i, j}) = +1;

            for (size_t k = 0; k < j; k++){
                problem.costOfTriple({i, j, k}) = -1;
            }
        }
    }

    problem.costOfEdge({2, 1}) = -3.0*(numberOfNodes -3)*(numberOfNodes -4)/2 - 1 + 0.1;


    Persistency  persistency(problem);
    persistency.tripletEdgeJoinCriterion();

    auto edges = persistency.edges1();

    test(edges.empty());
    test(edges.find({2, 1}) == edges.end());
}

void testTripletEdgeJoinCriterionFulfilledArbitraryNumberOfNodes2(int numberOfNodes = 5){
    Problem problem(numberOfNodes);

    for (size_t i = 0; i < numberOfNodes; ++i){
        for (size_t j = 0; j < i; ++j){
            problem.costOfEdge({i, j}) = +1;

            for (size_t k = 0; k < j; k++){
                problem.costOfTriple({i, j, k}) = -1;
            }
        }
    }

    problem.costOfEdge({2, 1}) = -3.0*(numberOfNodes -3)*(numberOfNodes -4)/2 - 1 - 0.1;


    Persistency  persistency(problem);
    persistency.tripletEdgeJoinCriterion();

    auto edges = persistency.edges1();

    test(edges.size() == 1);
    test(edges.find({2, 1}) != edges.end());
}


int main() {
    testTripletEdgeJoinCriterionFulfilledSingleTriplet();
    testTripletEdgeJoinCriterionUnfulfilled();
    testTripletEdgeJoinCriterionFulfilledOneEdge();
    testTripletEdgeJoinCriterionOneEdgeFulfilled();
    testTripletEdgeJoinCriterionTwoEdgesFulfilled();
    testTripletEdgeJoinCriterionUnfulfilledFiveNodes();
    testTripletEdgeJoinCriterionFulfilledArbitraryNumberOfNodes(5);
    testTripletEdgeJoinCriterionFulfilledArbitraryNumberOfNodes(10);
    testTripletEdgeJoinCriterionFulfilledArbitraryNumberOfNodes(20);
    testTripletEdgeJoinCriterionUnfulfilledArbitraryNumberOfNodes(5);
    testTripletEdgeJoinCriterionUnfulfilledArbitraryNumberOfNodes(10);
    testTripletEdgeJoinCriterionUnfulfilledArbitraryNumberOfNodes(20);
    testTripletEdgeJoinCriterionUnfulfilledArbitraryNumberOfNodes2(5);
    testTripletEdgeJoinCriterionUnfulfilledArbitraryNumberOfNodes2(10);
    testTripletEdgeJoinCriterionUnfulfilledArbitraryNumberOfNodes2(20);
    testTripletEdgeJoinCriterionFulfilledArbitraryNumberOfNodes2(5);
    testTripletEdgeJoinCriterionFulfilledArbitraryNumberOfNodes2(10);
    testTripletEdgeJoinCriterionFulfilledArbitraryNumberOfNodes2(20);
    return 0;
}
