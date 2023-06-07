#include <stdexcept>
#include <iostream>

#include "andres/graph/multicut-cubic/persistency.hxx"

inline void test(bool const condition) {
    if(!condition) throw std::logic_error("test failed");
}

typedef andres::graph::multicut_cubic::Problem<double> Problem;
typedef andres::graph::multicut_cubic::Persistency<double> Persistency;


void testCriterion3FulfilledSingletNode(){
    Problem problem(4);
    // should separate the vertex 0
    problem.costOfTriple({0, 1, 2}) = 1;
    problem.costOfTriple({0, 1, 3}) = 1;
    problem.costOfTriple({0, 2, 3}) = 1;
    problem.costOfTriple({1, 2, 3}) = -1;

    Persistency persistency(problem);
    persistency.findIndependentSubProblems();

    test(persistency.subProblems().size() == 2);
    test(persistency.subPersistencies().size() == 2);

    std::set<std::set<size_t>> expectedSubsets{
            std::set<size_t >{0},
            std::set<size_t >{1, 2, 3}
    };

    for (auto const& expected : expectedSubsets){
        test(persistency.subProblems().find(expected) != persistency.subProblems().end());
    }

    auto const & edges0 = persistency.edges0();

    test(edges0.find({1, 0}) != edges0.end());
    test(edges0.find({2, 0}) != edges0.end());
    test(edges0.find({3, 0}) != edges0.end());
    test(edges0.size() == 3);


}


void testCriterion3Unfulfilled(){
    Problem problem(5);
    // should separate the vertex 0
    for (size_t i = 0; i < 5; i++){
        for (size_t j = 0; j < i; j++){
            for (size_t k = 0; k < j; k++){
                problem.costOfTriple({i, j, k}) = +1;
            }
        }
    }

    problem.costOfTriple({0, 1, 2}) = -1;
    problem.costOfTriple({0, 3, 4}) = -1;

    Persistency persistency(problem);
    persistency.findIndependentSubProblems();

    test(persistency.subProblems().size() == 1);
    test(persistency.subPersistencies().empty());

    std::set<std::set<size_t>> expectedSubsets{
            std::set<size_t >{0, 1, 2, 3, 4},
    };

    for (auto const& expected : expectedSubsets){
        test(persistency.subProblems().find(expected) != persistency.subProblems().end());
    }

    auto const & edges0 = persistency.edges0();

    test(edges0.empty());
}

void testCriterion3FulfilledTripletNodes(){
    Problem problem(6);
    // should separate the vertex 0
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 3; j < 6; j++) {
            for (size_t k = 0; k < 6; k++) {
                if (i == j || j == k || i == k) continue;

                if (i > j && j > k){
                    problem.costOfTriple({i, j, k}) = +1;
                } else if (j > i && i > k){
                    problem.costOfTriple({j, i, k}) = +1;
                }
                else if (k > i && i > j){
                    problem.costOfTriple({k, i, j}) = 1;
                }
                else {
                    problem.costOfTriple({i, j, k}) = 1;
                }
            }
        }
    }
    problem.costOfTriple({0, 1, 2}) = -1;
    problem.costOfTriple({3, 4, 5}) = -1;


    Persistency persistency(problem);
    persistency.findIndependentSubProblems();

    test(persistency.subProblems().size() == 2);
    test(persistency.subPersistencies().size() == 2);

    std::set<std::set<size_t>> expectedSubsets{
            std::set<size_t >{0, 1, 2},
            std::set<size_t >{3, 4, 5}
    };

    for (auto const& expected : expectedSubsets){
        test(persistency.subProblems().find(expected) != persistency.subProblems().end());
    }

    for (auto const & subPersistency : persistency.subPersistencies()){
        test(subPersistency.edges0().empty());
        test(subPersistency.edges1().empty());
        test(subPersistency.remainingNodes() == 3);
        test(subPersistency.remainingVariables() == 3);
        test(subPersistency.remainingTriples() == 1);
    }

    auto const & edges0 = persistency.edges0();

    test(edges0.find({3, 0}) != edges0.end());
    test(edges0.find({3, 1}) != edges0.end());
    test(edges0.find({3, 2}) != edges0.end());
    test(edges0.find({4, 0}) != edges0.end());
    test(edges0.find({4, 1}) != edges0.end());
    test(edges0.find({4, 2}) != edges0.end());
    test(edges0.find({5, 0}) != edges0.end());
    test(edges0.find({5, 1}) != edges0.end());
    test(edges0.find({5, 2}) != edges0.end());
    test(edges0.size() == 9);

}


void testOnlyEdgeCosts(){
    Problem problem(6);
    // should separate the vertex 0
    for (size_t i = 0; i < problem.numberOfElements(); ++i){
        for (size_t j = 0; j < i; ++j){
            problem.costOfEdge({i, j}) = -1;

            for (size_t k = 0; k < j; ++k){
                problem.costOfTriple({i, j, k}) = 0;
            }
        }
    }


    Persistency persistency(problem);
    persistency.findIndependentSubProblems();

    test(persistency.subProblems().size() == 1);
    test(persistency.subPersistencies().empty());
    test(persistency.edges0().empty());
}

int main() {
    testCriterion3FulfilledSingletNode();
    testCriterion3FulfilledTripletNodes();
    testCriterion3Unfulfilled();
    testOnlyEdgeCosts();
    return 0;
}
