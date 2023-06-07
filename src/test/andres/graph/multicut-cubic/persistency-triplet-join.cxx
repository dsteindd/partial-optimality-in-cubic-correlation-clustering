#include <stdexcept>
#include <iostream>

#include "andres/graph/multicut-cubic/persistency.hxx"

inline void test(bool const condition) {
    if(!condition) throw std::logic_error("test failed");
}

typedef andres::graph::multicut_cubic::Problem<double> Problem;
typedef andres::graph::multicut_cubic::Persistency<double> Persistency;

void testNotJoinedThreeNodes(){
    {
        Problem problem(3);
        problem.costOfTriple({2, 1, 0}) = -2;
        problem.costOfEdge({2, 0}) = -2;
        problem.costOfEdge({1, 0}) = -2;
        problem.costOfEdge({2, 1}) = +7;

        Persistency persistency(problem);

        persistency.tripletJoinCriterion();

        auto edges1 = persistency.edges1();

        test(edges1.find({1, 0}) == edges1.end());
        test(edges1.find({2, 0}) == edges1.end());
        test(edges1.find({2, 1}) == edges1.end());
    }

    {
        Problem problem(3);
        problem.costOfTriple({2, 1, 0}) = -2;
        problem.costOfEdge({2, 0}) = -2;
        problem.costOfEdge({1, 0}) = -2;
        problem.costOfEdge({2, 1}) = +5;

        Persistency persistency(problem);

        persistency.tripletJoinCriterion();

        auto edges1 = persistency.edges1();

        test(edges1.find({1, 0}) == edges1.end());
        test(edges1.find({2, 0}) == edges1.end());
        test(edges1.find({2, 1}) == edges1.end());
    }
}

void testJoinedThreeNodes(){
    Problem problem(3);
    problem.costOfTriple({2, 1, 0}) = -2;
    problem.costOfEdge({2, 0}) = -2;
    problem.costOfEdge({1, 0}) = -2;
    problem.costOfEdge({2, 1}) = +3;

    Persistency persistency(problem);

    persistency.tripletJoinCriterion();

    auto edges1 = persistency.edges1();

    test(edges1.find({1, 0}) != edges1.end());
    test(edges1.find({2, 0}) != edges1.end());
    test(edges1.find({2, 1}) != edges1.end());
}

void testJoinedFourNodes(){
    Problem problem(4);
    problem.costOfTriple({2, 1, 0}) = -4;
    problem.costOfTriple({3, 1, 0}) = +2;
    problem.costOfTriple({3, 2, 0}) = +2;
    problem.costOfTriple({3, 2, 1}) = +2;
    problem.costOfEdge({2, 0}) = -4;
    problem.costOfEdge({1, 0}) = -4;
    problem.costOfEdge({2, 1}) = -4;
    problem.costOfEdge({3, 0}) = +2;
    problem.costOfEdge({3, 1}) = +2;
    problem.costOfEdge({3, 2}) = +2;

    Persistency persistency(problem);

    persistency.tripletJoinCriterion();

    auto edges1 = persistency.edges1();

    test(edges1.find({1, 0}) != edges1.end());
    test(edges1.find({2, 0}) != edges1.end());
    test(edges1.find({2, 1}) != edges1.end());
}


void testNotJoinedFourNodes(){
    Problem problem(4);
    problem.costOfTriple({2, 1, 0}) = -2;
    problem.costOfTriple({3, 1, 0}) = +2;
    problem.costOfTriple({3, 2, 0}) = +2;
    problem.costOfTriple({3, 2, 1}) = +2;
    problem.costOfEdge({2, 0}) = -2;
    problem.costOfEdge({1, 0}) = -2;
    problem.costOfEdge({2, 1}) = -2;
    problem.costOfEdge({3, 0}) = +2;
    problem.costOfEdge({3, 1}) = +2;
    problem.costOfEdge({3, 2}) = +2;

    Persistency persistency(problem);

    persistency.tripletJoinCriterion();

    auto edges1 = persistency.edges1();

    test(edges1.find({1, 0}) == edges1.end());
    test(edges1.find({2, 0}) == edges1.end());
    test(edges1.find({2, 1}) == edges1.end());
}

void testNotJoinedBySingletonButByDoubletFourNodes(){
    Problem problem(4);
    problem.costOfTriple({2, 1, 0}) = -2.3;
    problem.costOfTriple({3, 1, 0}) = +2;
    problem.costOfTriple({3, 2, 0}) = +2;
    problem.costOfTriple({3, 2, 1}) = -1;
    problem.costOfEdge({2, 0}) = -2.3;
    problem.costOfEdge({1, 0}) = -2.3;
    problem.costOfEdge({2, 1}) = -2.3;
    problem.costOfEdge({3, 0}) = -6;
    problem.costOfEdge({3, 1}) = -1;
    problem.costOfEdge({3, 2}) = -1;

    Persistency persistency(problem);

    persistency.tripletJoinCriterion();

    auto edges1 = persistency.edges1();

    test(edges1.find({1, 0}) != edges1.end());
    test(edges1.find({2, 0}) != edges1.end());
    test(edges1.find({2, 1}) != edges1.end());
}



int main() {
    testNotJoinedThreeNodes();
    testJoinedThreeNodes();
    testJoinedFourNodes();
    testNotJoinedFourNodes();
    testNotJoinedBySingletonButByDoubletFourNodes();
}
