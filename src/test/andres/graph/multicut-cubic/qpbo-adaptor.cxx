#include <stdexcept>

#include "andres/graph/multicut-cubic/persistency.hxx"
#include "andres/graph/multicut-cubic/qpbo-cut-problem-adaptor.hxx"

inline void test(bool const condition) {
    if(!condition) throw std::logic_error("test failed");
}



void testCostAdaption(size_t const numberOfElements){
    typedef andres::graph::multicut_cubic::Problem<int> Problem;
    typedef andres::graph::multicut_cubic::QPBOCutProblemAdaptor<int> Adaptor;

    Problem problem(numberOfElements);
    problem.costConstant() = 4;

    for (size_t i = 0; i < numberOfElements; ++i){
        for (size_t j = 0; j < i; ++j){
            for (size_t k = 0; k < j; ++k){
                problem.costOfTriple({i, j, k}) = 1;
            }
        }
    }
    for (size_t i = 0; i < numberOfElements; ++i){
        for (size_t j = 0; j < i; ++j){
            problem.costOfEdge({i, j}) = 1;
        }
    }

    Adaptor adaptor(problem);

    // vertex cost should be |V| - 1 + (|V|-1)(|V|-2)/2 = |V|(|V|-1)/2
    // for |V| = 4 -> 4*3/2 = 6
    for (size_t i = 0; i < numberOfElements; ++i){
        test(adaptor.vertexCost(i) == numberOfElements*(numberOfElements - 1)/2);
    }

    // edge cost should be -2 + -|V|+2 = -|V|
    for (size_t i = 0; i < numberOfElements; ++i){
        for (size_t j = 0; j < i; ++j)
            test(adaptor.edgeCost(i, j) == -numberOfElements);
    }

    test(adaptor.constant() == 4);
}

void testCostAdaption2(size_t const numberOfElements){
    typedef andres::graph::multicut_cubic::Problem<int> Problem;
    typedef andres::graph::multicut_cubic::QPBOCutProblemAdaptor<int> Adaptor;

    Problem problem(numberOfElements);
    problem.costConstant() = 4;

    const int n = static_cast<int>(numberOfElements);

    for (size_t i = 0; i < numberOfElements; ++i){
        for (size_t j = 0; j < i; ++j){
            for (size_t k = 0; k < j; ++k){
                problem.costOfTriple({i, j, k}) = i + j + k;
            }
        }
    }
    for (size_t i = 0; i < numberOfElements; ++i){
        for (size_t j = 0; j < i; ++j){
            problem.costOfEdge({i, j}) = i + j;
        }
    }

    Adaptor adaptor(problem);

    // vertex cost for i should be:
    for (size_t i = 0; i < numberOfElements; ++i){
        const int value = static_cast<int>(i)*(n-1)*(n-2)/2 + n*(n-1)*(n-1)/2;
        test(adaptor.vertexCost(i) == value);
    }

    // edge cost should be -2 + -|V|+2 = -|V|
    for (size_t i = 0; i < numberOfElements; ++i){
        for (size_t j = 0; j < i; ++j) {
            const int value = -n * (n - 1) / 2 - (n - 1) * (static_cast<int>(i) + static_cast<int>(j));
            test(adaptor.edgeCost(i, j) == value);
        }
    }

    test(adaptor.constant() == 4);
}

void testCostAdaptionWithZeroConstraint(size_t const numberOfElements){
    typedef andres::graph::multicut_cubic::Problem<int> Problem;
    typedef andres::graph::multicut_cubic::QPBOCutProblemAdaptor<int> Adaptor;

    Problem problem(numberOfElements);
    problem.costConstant() = 4;

    const int n = static_cast<int>(numberOfElements);

    for (size_t i = 0; i < numberOfElements; ++i){
        for (size_t j = 0; j < i; ++j){
            for (size_t k = 0; k < j; ++k){
                problem.costOfTriple({i, j, k}) = i + j + k;
            }
        }
    }
    for (size_t i = 0; i < numberOfElements; ++i){
        for (size_t j = 0; j < i; ++j){
            problem.costOfEdge({i, j}) = i + j;
        }
    }

    const size_t zeroVertex = numberOfElements / 2;

    Adaptor adaptor(problem, std::vector<size_t>(), std::vector<size_t>({zeroVertex}));

    test(adaptor.numberOfVertices() == numberOfElements - 1);

    // vertex cost for i should be:
    for (size_t i = 0; i < adaptor.numberOfVertices(); ++i){
        if (i == zeroVertex){
            const size_t last = adaptor.numberOfVertices();
            const int value = static_cast<int>(last)*(n-1)*(n-2)/2 + n*(n-1)*(n-1)/2;
            test(adaptor.vertexCost(i) == value);
        } else {
            const int value = static_cast<int>(i)*(n-1)*(n-2)/2 + n*(n-1)*(n-1)/2;
            test(adaptor.vertexCost(i) == value);
        }
    }

    // edge cost should be -2 + -|V|+2 = -|V|
    for (size_t i = 0; i < adaptor.numberOfVertices(); ++i){
        size_t i1 = i;
        if (i == zeroVertex){
            i1 = adaptor.numberOfVertices();
        }
        for (size_t j = 0; j < i; ++j) {
            size_t j1 = j;
            if (j == zeroVertex){
                j1 = adaptor.numberOfVertices();
            }

            const int value = -n * (n - 1) / 2 - (n - 1) * (static_cast<int>(i1) + static_cast<int>(j1));
            test(adaptor.edgeCost(i, j) == value);
        }
    }

    test(adaptor.constant() == 4);
}

void testCostAdaptionWithOneConstraint(size_t const numberOfElements){
    typedef andres::graph::multicut_cubic::Problem<int> Problem;
    typedef andres::graph::multicut_cubic::QPBOCutProblemAdaptor<int> Adaptor;

    Problem problem(numberOfElements);
    problem.costConstant() = 4;

    const int n = static_cast<int>(numberOfElements);

    for (size_t i = 0; i < numberOfElements; ++i){
        for (size_t j = 0; j < i; ++j){
            for (size_t k = 0; k < j; ++k){
                problem.costOfTriple({i, j, k}) = i + j + k;
            }
        }
    }
    for (size_t i = 0; i < numberOfElements; ++i){
        for (size_t j = 0; j < i; ++j){
            problem.costOfEdge({i, j}) = i + j;
        }
    }

    const size_t oneVertex = numberOfElements / 2;

    Adaptor adaptor(problem, std::vector<size_t>({oneVertex}), std::vector<size_t>());

    test(adaptor.numberOfVertices() == numberOfElements - 1);

    // vertex cost for i should be:
    for (size_t i = 0; i < adaptor.numberOfVertices(); ++i){
        size_t i1 = i;
        if (i == oneVertex){
            i1 = adaptor.numberOfVertices();
        }


        const int value = static_cast<int>(i1)*(n-1)*(n-2)/2 + n*(n-1)*(n-1)/2 -n * (n - 1) / 2 - (n - 1) * (static_cast<int>(i1) + static_cast<int>(oneVertex));
        test(adaptor.vertexCost(i) == value);
    }

//     edge cost should be -2 + -|V|+2 = -|V|
    for (size_t i = 0; i < adaptor.numberOfVertices(); ++i){
        size_t i1 = i;
        if (i == oneVertex){
            i1 = adaptor.numberOfVertices();
        }

        for (size_t j = 0; j < i; ++j) {
            size_t j1 = j;
            if (j == oneVertex){
                j1 = adaptor.numberOfVertices();
            }

            const int value = -n * (n - 1) / 2 - (n - 1) * (static_cast<int>(i1) + static_cast<int>(j1));
            test(adaptor.edgeCost(i, j) == value);
        }
    }
    const int expectedConstant = 4 + static_cast<int>(oneVertex)*(n-1)*(n-2)/2 + n*(n-1)*(n-1)/2;
    test(adaptor.constant() == expectedConstant);
}


int main() {
    testCostAdaption(5);
    testCostAdaption(10);
    testCostAdaption(20);
    testCostAdaption2(5);
    testCostAdaption2(10);
    testCostAdaption2(20);
    testCostAdaptionWithZeroConstraint(5);
    testCostAdaptionWithZeroConstraint(10);
    testCostAdaptionWithZeroConstraint(20);
    testCostAdaptionWithOneConstraint(5);
    testCostAdaptionWithOneConstraint(10);
    testCostAdaptionWithOneConstraint(20);
    return 0;
}
