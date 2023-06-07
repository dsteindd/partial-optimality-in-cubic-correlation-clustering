#include <stdexcept>
#include <iostream>

#include "andres/graph/multicut-cubic/persistency.hxx"

inline void test(bool const condition) {
    if(!condition) throw std::logic_error("test failed");
}

typedef andres::graph::multicut_cubic::Problem<double> Problem;
typedef andres::graph::multicut_cubic::Persistency<double> Persistency;


void testPositivelySeparatedClusters(std::vector<size_t> clusterSizes){
    size_t numClusters = clusterSizes.size();
    size_t numberOfElements = 0;
    for (auto & clusterSize: clusterSizes){
        numberOfElements += clusterSize;
    }

    std::vector<size_t > labels(numberOfElements);
    size_t offset = 0;
    for (size_t clusterIndex = 0; clusterIndex < numClusters; ++clusterIndex){
        size_t clusterSize = clusterSizes[clusterIndex];
        for (size_t index = offset; index < offset + clusterSize; ++index){
            labels[index] = clusterIndex;
        }
        offset += clusterSize;
    }

    Problem problem(numberOfElements);

    // set costs
    for (size_t i = 0; i < numberOfElements; ++i){
        for (size_t j = 0; j < i; ++j){
            if (labels[i] == labels[j]){
                problem.costOfEdge({i, j}) = -1;
            } else {
                problem.costOfEdge({i, j}) = +1;
            }


            for (size_t k = 0; k < j; ++k){
                if (labels[i] == labels[j] && labels[j] == labels[k]){
                    problem.costOfTriple({i, j, k}) = -1;
                } else {
                    problem.costOfTriple({i, j, k}) = +1;
                }
            }
        }
    }

    Persistency persistency(problem);
    persistency.subsetJoinCriterion();

    test(persistency.remainingNodes() == clusterSizes.size());
    test(persistency.remainingVariables() == numClusters*(numClusters - 1) / 2);

    std::set<Persistency::IndexEdge > const edges1 = persistency.edges1();

    size_t expectedNumEdges1 = 0;
    for (auto & clusterSize: clusterSizes){
        expectedNumEdges1 += clusterSize * (clusterSize - 1) / 2;
    }
    test(edges1.size() == expectedNumEdges1);

    offset = 0;
    for (size_t clusterIndex = 0; clusterIndex < numClusters; ++clusterIndex){
        size_t clusterSize = clusterSizes[clusterIndex];

        for (size_t i = offset; i < offset + clusterSize; ++i){
            for (size_t j = offset; j < i; ++j){
                test(edges1.find({i, j}) != edges1.end());
            }
        }

        offset += clusterSize;
    }
}

void testClusterSeparationWithNegativeCostsBetweenClusters(std::vector<size_t> clusterSizes){
    size_t numClusters = clusterSizes.size();
    assert(numClusters > 1);
    for (auto & clusterSize : clusterSizes){
        assert(clusterSize > 1);
    }


    size_t numberOfElements = 0;
    size_t minClusterSize = std::numeric_limits<unsigned long>().max();
    for (auto & clusterSize: clusterSizes){
        numberOfElements += clusterSize;
        minClusterSize = std::min(clusterSize, minClusterSize);
    }

    // should be the exact threshold
    double epsilon = (static_cast<double>(minClusterSize) - 1) / (static_cast<double>(numberOfElements) - static_cast<double>(minClusterSize)) / 2;

    // to be safe
    epsilon = epsilon / 2;

    std::vector<size_t > labels(numberOfElements);
    size_t offset = 0;
    for (size_t clusterIndex = 0; clusterIndex < numClusters; ++clusterIndex){
        size_t clusterSize = clusterSizes[clusterIndex];
        for (size_t index = offset; index < offset + clusterSize; ++index){
            labels[index] = clusterIndex;
        }
        offset += clusterSize;
    }

    Problem problem(numberOfElements);

    // set costs
    for (size_t i = 0; i < numberOfElements; ++i){
        for (size_t j = 0; j < i; ++j){
            if (labels[i] == labels[j]){
                problem.costOfEdge({i, j}) = -1;
            } else {
                problem.costOfEdge({i, j}) = -epsilon;
            }


            for (size_t k = 0; k < j; ++k){
                if (labels[i] == labels[j] && labels[j] == labels[k]){
                    problem.costOfTriple({i, j, k}) = -1;
                } else {
                    problem.costOfTriple({i, j, k}) = +1;
                }
            }
        }
    }

    Persistency persistency(problem);
    persistency.subsetJoinCriterion();

    test(persistency.remainingNodes() == clusterSizes.size());
    test(persistency.remainingVariables() == numClusters*(numClusters - 1) / 2);

    std::set<Persistency::IndexEdge > const edges1 = persistency.edges1();

    size_t expectedNumEdges1 = 0;
    for (auto & clusterSize: clusterSizes){
        expectedNumEdges1 += clusterSize * (clusterSize - 1) / 2;
    }
    test(edges1.size() == expectedNumEdges1);

    offset = 0;
    for (size_t clusterIndex = 0; clusterIndex < numClusters; ++clusterIndex){
        size_t clusterSize = clusterSizes[clusterIndex];

        for (size_t i = offset; i < offset + clusterSize; ++i){
            for (size_t j = offset; j < i; ++j){
                test(edges1.find({i, j}) != edges1.end());
            }
        }

        offset += clusterSize;
    }
}


void testClustersCannotBeMergedDueToTooHighInterClusterNegativeCosts(std::vector<size_t> clusterSizes){
    size_t numClusters = clusterSizes.size();
    assert(numClusters > 1);
    for (auto & clusterSize : clusterSizes){
        assert(clusterSize > 1);
    }


    size_t numberOfElements = 0;
    size_t maxClusterSize = 0;
    for (auto & clusterSize: clusterSizes){
        numberOfElements += clusterSize;
        maxClusterSize = std::max(clusterSize, maxClusterSize);
    }

    // should be the exact threshold
    double epsilon = (static_cast<double>(maxClusterSize) - 1) / (static_cast<double>(numberOfElements) - static_cast<double>(maxClusterSize)) / 2;

    epsilon = epsilon + 0.1;

    std::vector<size_t > labels(numberOfElements);
    size_t offset = 0;
    for (size_t clusterIndex = 0; clusterIndex < numClusters; ++clusterIndex){
        size_t clusterSize = clusterSizes[clusterIndex];
        for (size_t index = offset; index < offset + clusterSize; ++index){
            labels[index] = clusterIndex;
        }
        offset += clusterSize;
    }

    Problem problem(numberOfElements);

    // set costs
    for (size_t i = 0; i < numberOfElements; ++i){
        for (size_t j = 0; j < i; ++j){
            if (labels[i] == labels[j]){
                problem.costOfEdge({i, j}) = -1;
            } else {
                problem.costOfEdge({i, j}) = -epsilon;
            }


            for (size_t k = 0; k < j; ++k){
                if (labels[i] == labels[j] && labels[j] == labels[k]){
                    problem.costOfTriple({i, j, k}) = -1;
                } else {
                    problem.costOfTriple({i, j, k}) = +1;
                }
            }
        }
    }

    Persistency persistency(problem);
    persistency.subsetJoinCriterion();

    test(persistency.remainingNodes() == problem.numberOfElements());
    test(persistency.remainingVariables() == problem.numberOfElements()*(problem.numberOfElements() - 1) / 2);

    std::set<Persistency::IndexEdge > const edges1 = persistency.edges1();
    test(edges1.empty());
}




int main() {
    testPositivelySeparatedClusters({10, 10, 10});
    testPositivelySeparatedClusters({5, 15, 20});
    testPositivelySeparatedClusters({1, 1, 1});
    testClusterSeparationWithNegativeCostsBetweenClusters({10, 10, 10});
    testClusterSeparationWithNegativeCostsBetweenClusters({5, 15, 20});
    testClusterSeparationWithNegativeCostsBetweenClusters({2, 2, 2});
    testClustersCannotBeMergedDueToTooHighInterClusterNegativeCosts({10, 10, 10});
    testClustersCannotBeMergedDueToTooHighInterClusterNegativeCosts({5, 15, 20});
    testClustersCannotBeMergedDueToTooHighInterClusterNegativeCosts({2, 2, 2});

}
