#include <stdexcept>
#include <iostream>

#include "andres/graph/multicut-cubic/persistency.hxx"

inline void test(bool const condition) {
    if(!condition) throw std::logic_error("test failed");
}

typedef andres::graph::multicut_cubic::Problem<double> Problem;
typedef andres::graph::multicut_cubic::Persistency<double> Persistency;

void testNotJoined(){
    Problem problem(3);
    problem.costOfTriple({2, 1, 0}) = -3;
    problem.costOfEdge({2, 0}) = 2.5;
    problem.costOfEdge({1, 0}) = -2;
    problem.costOfEdge({2, 1}) = 2.5;

    Persistency persistency(problem);

    persistency.edgeJoinCriterion();

    auto edges1 = persistency.edges1();

    test(edges1.size() == 0);
}

void testJoinedSingle(){
    Problem problem(3);
    problem.costOfTriple({2, 1, 0}) = -3;
    problem.costOfEdge({2, 0}) = 1.6;
    problem.costOfEdge({1, 0}) = -2;
    problem.costOfEdge({2, 1}) = 1.6;

    Persistency persistency(problem);

    persistency.edgeJoinCriterion();

    auto edges1 = persistency.edges1();

    test(edges1.size() == 1);
    test(edges1.find({1, 0}) != edges1.end());
}


void testJoinedDouble(){
    Problem problem(3);
    problem.costOfTriple({2, 1, 0}) = -3;
    problem.costOfEdge({2, 0}) = 1.4;
    problem.costOfEdge({1, 0}) = -2;
    problem.costOfEdge({2, 1}) = 1.4;

    Persistency persistency(problem);

    persistency.edgeJoinCriterion();

    auto edges1 = persistency.edges1();

    test(edges1.size() == 3);
    test(edges1.find({1, 0}) != edges1.end());
    test(edges1.find({2, 0}) != edges1.end());
    test(edges1.find({2, 1}) != edges1.end());
}

void testJoinedSingleZeroCost(){
    Problem problem(3);
    problem.costOfTriple({2, 1, 0}) = -1;
    problem.costOfEdge({2, 0}) = 0;
    problem.costOfEdge({1, 0}) = 0;
    problem.costOfEdge({2, 1}) = 1.5;

    Persistency persistency(problem);

    persistency.edgeJoinCriterion();

    auto edges1 = persistency.edges1();

    test(edges1.size() == 1);
    test(edges1.find({1, 0}) != edges1.end());
}


void testJoinedDoubleZeroCost(){
    Problem problem(3);
    problem.costOfTriple({2, 1, 0}) = -1;
    problem.costOfEdge({2, 0}) = 0;
    problem.costOfEdge({1, 0}) = 0.5;
    problem.costOfEdge({2, 1}) = 0;

    Persistency persistency(problem);

    persistency.edgeJoinCriterion();

    auto edges1 = persistency.edges1();

    test(edges1.size() == 3);
    test(edges1.find({1, 0}) != edges1.end());
    test(edges1.find({2, 0}) != edges1.end());
    test(edges1.find({2, 1}) != edges1.end());
}


int main() {
    testNotJoined();
    testJoinedSingle();
    testJoinedDouble();
    testJoinedSingleZeroCost();
    testJoinedDoubleZeroCost();

    return 0;
}
