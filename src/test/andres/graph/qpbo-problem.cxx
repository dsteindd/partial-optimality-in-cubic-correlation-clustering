#include <random>
#include <vector>
#include <stdexcept>

#include "andres/graph/qpbo-problem.hxx"

inline void test(const bool condition) {
    if(!condition) throw std::logic_error("test failed.");
}

void testCostSettingForProblem(size_t numberOfVariables) {
    typedef andres::graph::QPBOProblem<double> Problem;

    Problem problem(numberOfVariables);
    problem.constant() = 1;

    for (size_t v = 0; v < numberOfVariables; ++v){
        problem.vertexCost(v) = v+2;
    }

    for (size_t u = 0; u < numberOfVariables; ++u){
        for (size_t v = 0; v < u; ++v){
            problem.edgeCost(u, v) = 10*u + v + 2;
        }
    }

    test(problem.numberOfVertices() == numberOfVariables);
    test(problem.constant() == 1);

    for (size_t v = 0; v < numberOfVariables; ++v){
        test(problem.vertexCost(v) == v+2);
    }

    for (size_t u = 0; u < numberOfVariables; ++u){
        for (size_t v = 0; v < u; ++v){
            test(problem.edgeCost(u, v) == 10*u + v + 2);
        }
    }
}

int main() {
    testCostSettingForProblem(10);

    return 0;
}
