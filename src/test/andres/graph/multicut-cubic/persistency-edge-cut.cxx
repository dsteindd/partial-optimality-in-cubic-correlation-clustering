#include <stdexcept>

#include "andres/graph/multicut-cubic/persistency.hxx"

inline void test(bool const condition) {
    if(!condition) throw std::logic_error("test failed");
}

typedef andres::graph::multicut_cubic::Problem<double> Problem;
typedef andres::graph::multicut_cubic::Persistency<double> Persistency;

void testCutOnlyForSingleton(){
    Problem problem(3);
    problem.costOfTriple({2, 1, 0}) = -1;
    problem.costOfEdge({2, 0}) = -1;
    problem.costOfEdge({1, 0}) = 3;
    problem.costOfEdge({2, 1}) = 1;

    Persistency persistency(problem);

    persistency.edgeCutCriterion();

    test(persistency.edges0().find({1, 0}) != persistency.edges0().end());
    test(persistency.remainingNodes() == 3);
    test(persistency.remainingVariables() <= 2);
}

void testNotCut(){
    Problem problem(3);
    problem.costOfTriple({2, 1, 0}) = -1;
    problem.costOfEdge({2, 0}) = -1;
    problem.costOfEdge({1, 0}) = 1;
    problem.costOfEdge({2, 1}) = -2;

    Persistency persistency(problem);

    persistency.edgeCutCriterion();

    auto const & edges0 = persistency.edges0();

    test(edges0.find({1, 0}) == edges0.end());
    test(persistency.remainingNodes() == 3);
    test(persistency.remainingVariables() == 3);
}


void testCutForSingletonAndDoublet(){
    Problem problem(3);
    problem.costOfTriple({2, 1, 0}) = -1;
    problem.costOfEdge({2, 0}) = -1;
    problem.costOfEdge({1, 0}) = 3;
    problem.costOfEdge({2, 1}) = +1;

    Persistency persistency(problem);

    persistency.edgeCutCriterion();

    auto const & edges0 = persistency.edges0();

    test(edges0.find({1, 0}) != edges0.end());
    test(persistency.remainingNodes() == 3);
    test(persistency.remainingVariables() <= 2);
}

void testNotCutForSingletonNorBiggerSet(){
    Problem problem(4);
    problem.costOfTriple({3, 1, 0}) = -5;
    problem.costOfTriple({2, 1, 0}) = -5;
    problem.costOfTriple({3, 2, 0}) = -1;
    problem.costOfEdge({1, 0}) = 3;
    problem.costOfEdge({2, 1}) = -2;
    problem.costOfEdge({3, 1}) = -2;

    Persistency persistency(problem);

    persistency.edgeCutCriterion();
    auto const edges0 = persistency.edges0();


    test(edges0.find({1, 0}) == edges0.end());
    test(persistency.remainingNodes() == 4);
    test(persistency.remainingVariables() == 6);
}

void testCutBySingletonButNotByTriplet(){
    Problem problem(4);
    problem.costOfTriple({3, 1, 0}) = -1;
    problem.costOfTriple({2, 1, 0}) = -1;
    problem.costOfTriple({3, 2, 0}) = -1;
    problem.costOfEdge({1, 0}) = 5;
    problem.costOfEdge({2, 1}) = -2;
    problem.costOfEdge({3, 1}) = -2;

    Persistency persistency(problem);

    persistency.edgeCutCriterion();
    auto const edges0 = persistency.edges0();


    test(edges0.find({1, 0}) != edges0.end());
    test(persistency.remainingNodes() == 4);
    test(persistency.remainingVariables() <= 5);
}
void testCutWithZeroEdgeCostBySingleton(){
    Problem problem(3);
    problem.costOfTriple({2, 1, 0}) = 1;
    problem.costOfEdge({1, 0}) = 0;
    problem.costOfEdge({2, 1}) = -3;
    problem.costOfEdge({2, 0}) = +2;

    Persistency persistency(problem);

    persistency.edgeCutCriterion();
    auto const edges0 = persistency.edges0();


    test(edges0.find({1, 0}) != edges0.end());
    test(persistency.remainingNodes() == 3);
    test(persistency.remainingVariables() <= 2);
}

void testCutWithZeroEdgeCostByDoublet(){
    Problem problem(3);
    problem.costOfTriple({2, 1, 0}) = 1;
    problem.costOfEdge({1, 0}) = 0;
    problem.costOfEdge({2, 1}) = +1;
    problem.costOfEdge({2, 0}) = -1;

    Persistency persistency(problem);

    persistency.edgeCutCriterion();

    auto const edges0 = persistency.edges0();

    test(edges0.find({1, 0}) != edges0.end());
    test(persistency.remainingNodes() == 3);
    test(persistency.remainingVariables() <= 2);
}


int main() {
    testCutOnlyForSingleton();
    testNotCut();
    testCutForSingletonAndDoublet();
    testNotCutForSingletonNorBiggerSet();
    testCutBySingletonButNotByTriplet();
    testCutBySingletonButNotByTriplet();
    testCutWithZeroEdgeCostBySingleton();
    testCutWithZeroEdgeCostByDoublet();

    return 0;
}
