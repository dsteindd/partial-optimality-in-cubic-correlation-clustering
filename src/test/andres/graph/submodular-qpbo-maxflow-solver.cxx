#include <stdexcept>

#include "andres/graph/qpbo.hxx"

typedef andres::graph::qpbo::MaxFlowSolver MaxFlowSolver;

inline void test(bool const condition) {
    if(!condition) throw std::logic_error("test failed");
}


template<MaxFlowSolver Solver>
inline
void testQPBOProblemSolverByMaxFlow(){
    typedef andres::graph::QPBOProblem<double> Problem;

    // note: this is the same test case as the first in maxflow tests
    Problem problem(4);

    problem.vertexCost(0) = -10;
    problem.vertexCost(1) = -10;
    problem.vertexCost(2) = 13;
    problem.vertexCost(3) = 10;

    problem.edgeCost(2, 0) = -10;
    problem.edgeCost(3, 2) = -5;
    problem.edgeCost(3, 1) = -5;

    double optimalValue;

    andres::graph::qpbo::solveSubmodularQPBO<double, Solver>(problem, optimalValue);

    test(optimalValue == -20);
}


int main() {
    testQPBOProblemSolverByMaxFlow<MaxFlowSolver::BOOST_PUSH_RELABEL>();
    return 0;
}
