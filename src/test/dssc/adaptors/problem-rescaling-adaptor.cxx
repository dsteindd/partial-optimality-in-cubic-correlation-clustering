#include <array>

#include <andres/graph/multicut-cubic/persistency.hxx>
#include "dssc/cost-adjustment-adaptors/problem-rescaling-adaptor.hxx"

typedef andres::graph::multicut_cubic::Problem<double> Problem;
typedef dssc::RescalingCostType CostType;
typedef dssc::CostRescalingAdaptor<double, CostType::PAIRS_AND_TRIPLES> ProblemAdaptor;

inline void test(bool const condition) {
    if(!condition) throw std::logic_error("test failed.");
}

inline bool isClose(double value, double compareValue, double tolerance = 0.001){
    return (abs(value - compareValue) <= tolerance);
}

void testStandardCases(){
    double minCost = -3;
    double maxPositiveCost = 3;

    Problem problem(4);

    problem.costOfEdge({1, 0}) = -1;
    problem.costOfEdge({2, 0}) = -5;
    problem.costOfEdge({2, 1}) = 5;
    problem.costOfEdge({3, 0}) = 2;
    problem.costOfEdge({3, 1}) = 3;
    problem.costOfEdge({3, 2}) = 1;
    problem.costOfTriple({2, 1, 0}) = -4;
    problem.costOfTriple({3, 1, 0}) = -2;
    problem.costOfTriple({3, 2, 0}) = 2;
    problem.costOfTriple({3, 2, 1}) = 5;

    ProblemAdaptor adaptor(problem, minCost, maxPositiveCost);

    test(isClose(adaptor.costOfEdge({1, 0}), -3.0/5));
    test(isClose(adaptor.costOfEdge({2, 0}), -3));
    test(isClose(adaptor.costOfEdge({2, 1}), 3));
    test(isClose(adaptor.costOfEdge({3, 0}), 6.0/5));
    test(isClose(adaptor.costOfEdge({3, 1}), 9.0/5));
    test(isClose(adaptor.costOfEdge({3, 2}), 3.0/5));
    test(isClose(adaptor.costOfTriple({2, 1, 0}), -3));
    test(isClose(adaptor.costOfTriple({3, 1, 0}), -1.5));
    test(isClose(adaptor.costOfTriple({3, 2, 0}),6.0/5));
    test(isClose(adaptor.costOfTriple({3, 2, 1}), 3));
}

void testAllBelowThreshold(){
    double minCost = -3;
    double maxPositiveCost = 3;

    Problem problem(4);

    problem.costOfEdge({1, 0}) = -1;
    problem.costOfEdge({2, 0}) = -2;
    problem.costOfEdge({2, 1}) = -2;
    problem.costOfEdge({3, 0}) = -2;
    problem.costOfEdge({3, 1}) = -2;
    problem.costOfEdge({3, 2}) = -1;
    problem.costOfTriple({2, 1, 0}) = -2;
    problem.costOfTriple({3, 1, 0}) = -2;
    problem.costOfTriple({3, 2, 0}) = -2;
    problem.costOfTriple({3, 2, 1}) = -2;

    ProblemAdaptor adaptor(problem, minCost, maxPositiveCost);

    test(isClose(adaptor.costOfEdge({1, 0}), -1.5));
    test(isClose(adaptor.costOfEdge({2, 0}), -3));
    test(isClose(adaptor.costOfEdge({2, 1}), -3));
    test(isClose(adaptor.costOfEdge({3, 0}), -3));
    test(isClose(adaptor.costOfEdge({3, 1}), -3));
    test(isClose(adaptor.costOfEdge({3, 2}), -1.5));
    test(isClose(adaptor.costOfTriple({2, 1, 0}), -3));
    test(isClose(adaptor.costOfTriple({3, 1, 0}), -3));
    test(isClose(adaptor.costOfTriple({3, 2, 0}),-3));
    test(isClose(adaptor.costOfTriple({3, 2, 1}), -3));
}

void testAllAboveThreshold(){
    double minCost = -3;
    double maxPositiveCost = 3;

    Problem problem(4);

    problem.costOfEdge({1, 0}) = 1;
    problem.costOfEdge({2, 0}) = 5;
    problem.costOfEdge({2, 1}) = 5;
    problem.costOfEdge({3, 0}) = 2;
    problem.costOfEdge({3, 1}) = 3;
    problem.costOfEdge({3, 2}) = 1;
    problem.costOfTriple({2, 1, 0}) = 4;
    problem.costOfTriple({3, 1, 0}) = 2;
    problem.costOfTriple({3, 2, 0}) = 2;
    problem.costOfTriple({3, 2, 1}) = 5;

    ProblemAdaptor adaptor(problem, minCost, maxPositiveCost);

    test(isClose(adaptor.costOfEdge({1, 0}), 3.0/5));
    test(isClose(adaptor.costOfEdge({2, 0}), 3));
    test(isClose(adaptor.costOfEdge({2, 1}), 3));
    test(isClose(adaptor.costOfEdge({3, 0}), 6.0/5));
    test(isClose(adaptor.costOfEdge({3, 1}), 9.0/5));
    test(isClose(adaptor.costOfEdge({3, 2}), 3.0/5));
    test(isClose(adaptor.costOfTriple({2, 1, 0}), 12.0/5));
    test(isClose(adaptor.costOfTriple({3, 1, 0}), 6.0/5));
    test(isClose(adaptor.costOfTriple({3, 2, 0}),6.0/5));
    test(isClose(adaptor.costOfTriple({3, 2, 1}), 3));
}



int main() {
    testStandardCases();
    testAllBelowThreshold();
    testAllAboveThreshold();


    return 0;
}
