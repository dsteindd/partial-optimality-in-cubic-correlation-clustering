#include <stdexcept>
#include <iostream>

#include "andres/graph/multicut-cubic/persistency.hxx"

inline void test(bool const condition) {
    if(!condition) throw std::logic_error("test failed");
}

typedef andres::graph::multicut_cubic::Problem<double> Problem;
typedef andres::graph::multicut_cubic::Persistency<double> Persistency;

void testNotJoinedThreeNodes(){
    Problem problem(3);
    problem.costOfTriple({2, 1, 0}) = 3;
    problem.costOfEdge({2, 0}) = -2;
    problem.costOfEdge({1, 0}) = -2;
    problem.costOfEdge({2, 1}) = -2;

    Persistency persistency(problem);

    persistency.tripletCutCriterion();

    auto triples0 = persistency.triples0();

    test(triples0.find({2, 1, 0}) == triples0.end());
}

void testJoinedThreeNodes(){
    Problem  problem(3);
    problem.costOfTriple({2, 1, 0}) = 3;
    problem.costOfEdge({2, 0}) = -1;
    problem.costOfEdge({1, 0}) = -1;
    problem.costOfEdge({2, 1}) = 2.5;

    Persistency persistency(problem);

    persistency.tripletCutCriterion();

    auto triples0 = persistency.triples0();

    test(triples0.find({2, 1, 0}) != triples0.end());
}

void testNotCutFiveNodes(){
    Problem problem(5);
    const double epsilon = 1;
    const double delta = 1;

    for (size_t i = 0; i < 5; ++i){
        for (size_t j = 0; j < i; ++j){
            for (size_t k = 0; k < j; ++k){
                problem.costOfTriple({i, j, k}) = -epsilon;
            }
        }
    }

    for (size_t i = 0; i < 5; ++i){
        for (size_t j = 0; j < i; ++j){
            problem.costOfEdge({i, j}) = -epsilon;
        }
    }

    problem.costOfTriple({2, 1, 0}) = 2;
    problem.costOfEdge({2, 0}) = 2;
    problem.costOfEdge({1, 0}) = 2;

    problem.costOfTriple({4, 3, 0}) = -delta;

    Persistency persistency(problem);

    persistency.tripletCutCriterion();

    auto triples0 = persistency.triples0();

    test(triples0.find({2, 1, 0}) == triples0.end());
}

void testCutWithSingletonFiveNodes(){
    Problem problem(5);
    const double epsilon = 0.5;
    const double delta = 0.5;

    for (size_t i = 0; i < 5; ++i){
        for (size_t j = 0; j < i; ++j){
            for (size_t k = 0; k < j; ++k){
                problem.costOfTriple({i, j, k}) = -epsilon;
            }
        }
    }

    for (size_t i = 0; i < 5; ++i){
        for (size_t j = 0; j < i; ++j){
            problem.costOfEdge({i, j}) = -epsilon;
        }
    }

    problem.costOfTriple({2, 1, 0}) = 2;
    problem.costOfEdge({2, 0}) = 2;
    problem.costOfEdge({1, 0}) = 2;

    problem.costOfTriple({4, 3, 0}) = -delta;

    Persistency persistency(problem);

    persistency.tripletCutCriterion();

    auto triples0 = persistency.triples0();

    test(triples0.find({2, 1, 0}) != triples0.end());
}

void testCutWithTripletButNotWithSingletonFiveNodes(){
    Problem problem(5);
    const double epsilon = 0.25;
    const double delta = 5;

    for (size_t i = 0; i < 5; ++i){
        for (size_t j = 0; j < i; ++j){
            for (size_t k = 0; k < j; ++k){
                problem.costOfTriple({i, j, k}) = -epsilon;
            }
        }
    }

    for (size_t i = 0; i < 5; ++i){
        for (size_t j = 0; j < i; ++j){
            problem.costOfEdge({i, j}) = -epsilon;
        }
    }

    problem.costOfTriple({2, 1, 0}) = 2;
    problem.costOfEdge({2, 0}) = 2;
    problem.costOfEdge({1, 0}) = 2;

    problem.costOfTriple({4, 3, 0}) = -delta;

    Persistency persistency(problem);

    persistency.tripletCutCriterion();

    auto triples0 = persistency.triples0();

    test(triples0.find({2, 1, 0}) != triples0.end());
}



int main() {
    testNotJoinedThreeNodes();
    testJoinedThreeNodes();
    testNotCutFiveNodes();
    testCutWithSingletonFiveNodes();
    testCutWithTripletButNotWithSingletonFiveNodes();
}
