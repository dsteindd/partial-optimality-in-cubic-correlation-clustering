#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_CUBIC_PERSISTENCY_HXX
#define ANDRES_GRAPH_MULTICUT_CUBIC_PERSISTENCY_HXX

#include <cstddef>
#include <vector>
#include <set>
#include <array>

#include "problem.hxx"
#include "andres/graph/multicut-cubic/qpbo-cut-problem-adaptor.hxx"
#include "andres/graph/qpbo.hxx"
#include "chrono"
#include "boost/graph/stoer_wagner_min_cut.hpp"
#include "boost/graph/graph_traits.hpp"
#include "boost/graph/adjacency_list.hpp"
#include "boost/property_map/property_map.hpp"
#include "andres/graph/multicut-cubic/absolute-cost-problem-adaptor.hxx"
#include "andres/graph/multicut-cubic/negative-part-cost-problem-adaptor.hxx"


namespace andres {
    namespace graph {
        namespace multicut_cubic {

            template<class T = double>
            class Persistency {
            public:
                typedef T Rational;
                typedef Problem<Rational> ProblemType;
                typedef AbsoluteCostProblemAdaptor<Rational> AbsoluteCostAdaptor;
                typedef NegativePartCostProblemAdaptor<Rational> NegativePartCostAdaptor;
                typedef std::array<size_t, 3> IndexTriple;
                typedef std::set<IndexTriple> IndexTripleSet;
                typedef std::array<size_t, 2> IndexEdge;
                typedef std::set<IndexEdge> IndexEdgeSet;
                typedef std::set<std::set<size_t>> PartitionSubsets;
                typedef std::chrono::duration<double> Duration;

                // requirements for PROBLEM_VIEW
                // - typedef Rational
                // - size_t numberOfElements() const
                // - Rational cost(size_t const, size_t const, size_t const) const
                template<class PROBLEM_VIEW>
                Persistency(PROBLEM_VIEW const &);

                template<class PROBLEM_VIEW>
                Persistency(PROBLEM_VIEW const &, std::vector<size_t> const &);

                IndexTripleSet const &triples0() const;

                IndexEdgeSet edges1() const;

                IndexEdgeSet const edges0() const;

                PartitionSubsets const &subProblems() const;

                std::vector<Persistency<T>> const &subPersistencies() const;

                // tested
                bool findIndependentSubProblems();

                Duration const &totalDurationFindIndependentProblems() const;

                // tested
                void edgeCutCriterion();

                Duration const &totalDurationEdgeCutCriterion() const;

                // tested
                void tripletCutCriterion();

                Duration const &totalDurationTripletCutCriterion() const;

                // tested
                void edgeJoinCriterion();

                Duration const &totalDurationEdgeJoinCriterion() const;

                // tested
                void tripletJoinCriterion();

                Duration const &totalDurationTripletJoinCriterion() const;

                // tested
                void tripletEdgeJoinCriterion();

                Duration const &totalDurationTripletEdgeJoinCriterion() const;

                // tested
                void edgeSubgraphCriterion();

                Duration const &totalDurationEdgeSubgraphCriterion() const;

                // tested
                void tripletSubgraphCriterion();

                Duration const &totalDurationTripletSubgraphCriterion() const;

                // tested
                void subsetJoinCriterion();

                Duration const &totalDurationSubsetJoinCriterion() const;

                size_t remainingNodes() const;

                size_t remainingVariables() const;

                size_t remainingTriples() const;

            private:
                Duration findIndependentProblemsDuration_;
                Duration edgeCutDuration_;
                Duration tripletCutDuration_;
                Duration edgeJoinDuration_;
                Duration tripletJoinDuration_;
                Duration tripletEdgeJoinDuration_;
                Duration edgeSubgraphDuration_;
                Duration tripletSubgraphDuration_;
                Duration subsetJoinDuration_;


                // findIndependentSubsets Helper
                // tries to find a node which if inserted to 'component' conflicts the conditions in the criterion, i.e.  introduces some negative cost (either negative edge or triplet cost)
                // if it finds such a node it returns 'False' (!) as opposed to the standard 'try'-API
                bool tryFindNegativeNode(const std::set<size_t> &component, size_t *nextNode);

                // triplet edge join criterion
                bool tripletEdgeJoinCriterion(IndexEdge *mergeEdge);

                bool tripletEdgeJoinCriterion(IndexTriple triple, IndexEdge edge);

                Rational calculateTripletEdgeJoinCriterionThirdRightHandSide(IndexTriple const);

                Rational calculateTripletEdgeJoinCriterionFirstOrSecondLeftHandSide(IndexEdge, IndexEdge);

                // edge join criterion
                bool edgeJoinCriterion(IndexEdge *);

                // triplet join criterion
                bool tripletJoinCriterion(IndexTriple *triple);

                Rational calculateTripletJoinCriterionLeftHandSide(IndexTriple const, IndexEdge const);

                // edge subgraph criterion
                bool checkEdgeSubgraph(IndexEdge *mergeEdge);

                // triplet subgraph criterion
                bool checkTripletSubgraph(IndexTriple triplet, IndexEdge edge);

                bool checkTripletSubgraph(IndexEdge *mergeEdge);


                // subset join criterion
                // checks whether a node can be added to a subset such that the resulting subset is negatively connected
                // returns also the negative cost offset
                std::pair<bool, Rational> canMergeSubsetJoin(std::set<size_t> const &, size_t const &);

                bool subsetJoinCriterion(std::set<size_t> *subset);


                // single node seed
                std::set<size_t> findMaximalComponentWithOnlyNegativeWeights(size_t const);

                // arbirary seed set
                std::set<size_t> findMaximalComponentWithOnlyNegativeWeights(std::set<size_t> const &);

                Rational subsetJoinCalculateLeftHandSide(std::set<size_t> const &);

                Rational subsetJoinCalculateRightHandSide(std::set<size_t> const &);

                // general helper methods
                Rational positivePart(const Rational);

                Rational negativePart(const Rational);

                template<class PROBLEM_VIEW>
                Rational
                getMinimalSeparationCost(PROBLEM_VIEW const &, std::vector<size_t> const &,
                                         std::vector<size_t> const &);

                // initializer props
                ProblemType problem_;

                // props to hold state of persistent variables
                IndexTripleSet triples0_;
                IndexEdgeSet edges0_;
                IndexEdgeSet edges1_;
                PartitionSubsets subProblems_;
                std::vector<Persistency> subProblemPersistencies_;
                // exclusively used for communicating to the client the correct node labels
                std::vector<size_t> nodeNames_;
            };

            template<class T>
            const std::vector<Persistency<T>> &Persistency<T>::subPersistencies() const {
                return subProblemPersistencies_;
            }


            template<class T>
            template<class PROBLEM_VIEW>
            Persistency<T>::Persistency(
                    PROBLEM_VIEW const &problemView
            )
                    :   problem_(problemView.numberOfElements()),
                        triples0_(),
                        edges1_(),
                        edges0_(),
                        findIndependentProblemsDuration_(0),
                        edgeCutDuration_(0),
                        tripletCutDuration_(0),
                        edgeJoinDuration_(0),
                        tripletJoinDuration_(0),
                        tripletEdgeJoinDuration_(0),
                        edgeSubgraphDuration_(0),
                        tripletSubgraphDuration_(0),
                        subsetJoinDuration_(0),
                        nodeNames_(problemView.numberOfElements()) {
                // read triplet costs
                for (size_t j = 0; j < problem_.numberOfElements(); ++j)
                    for (size_t k = 0; k < j; ++k)
                        for (size_t l = 0; l < k; ++l) {
                            problem_.costOfTriple({j, k, l}) = problemView.costOfTriple({j, k, l});
                        }

                // init edge costs
                for (size_t j = 0; j < problem_.numberOfElements(); ++j) {
                    for (size_t k = 0; k < j; k++) {
                        problem_.costOfEdge({j, k}) = problemView.costOfEdge({j, k});
                    }
                }

                for (size_t i = 0; i < problem_.numberOfElements(); ++i) {
                    nodeNames_[i] = i;
                }
            }

            template<class T>
            template<class PROBLEM_VIEW>
            Persistency<T>::Persistency(
                    PROBLEM_VIEW const &problemView,
                    std::vector<size_t> const &nodeNames
            )
                    :
                    Persistency(problemView) {
                nodeNames_ = nodeNames;
            }

            template<class T>
            inline
            typename Persistency<T>::IndexTripleSet const &
            Persistency<T>::triples0() const {
                return triples0_;
            }

            template<class T>
            inline
            typename Persistency<T>::IndexEdgeSet const
            Persistency<T>::edges0() const {
                // translate for client with provided names
                IndexEdgeSet returnEdges0 = {};
                for (IndexEdge const &edge: edges0_) {
                    returnEdges0.insert({nodeNames_[edge[0]], nodeNames_[edge[1]]});
                }

                // need to add all from independent sub persistencies
                if (!subProblemPersistencies_.empty()) {
                    for (auto &subPersistency: subProblemPersistencies_) {
                        IndexEdgeSet const subEdges0 = subPersistency.edges0();
                        returnEdges0.insert(subEdges0.begin(), subEdges0.end());
                    }
                }

                return returnEdges0;
            }

            template<class T>
            inline
            typename Persistency<T>::IndexEdgeSet
            Persistency<T>::edges1() const {

                // calculate from current partition
                std::vector<int> r(static_cast<std::size_t>(problem_.partition().numberOfElements()));

                problem_.partition().elementLabeling(r.begin());

                // get subsets
                std::map<size_t, std::set<size_t>> cLabelToSubset;
                size_t element = 0;
                for (size_t cLabel: r) {
                    if (cLabelToSubset.find(cLabel) == cLabelToSubset.end()) {
                        // key does not yet exist, assign empty set
                        cLabelToSubset[cLabel] = std::set<size_t>();
                    }
                    cLabelToSubset[cLabel].insert(element);
                    element++;
                }

                IndexEdgeSet returnEdges1 = {};

                for (auto const &[_, val]: cLabelToSubset) {
                    for (auto &i: val) {
                        for (auto &j: val) {
                            if (i <= j) continue;

                            // translate internal node labels to client side node labels
                            returnEdges1.insert({nodeNames_[i], nodeNames_[j]});

                        }
                    }
                }

                // add the 1-edges from the sub problems
                if (!subProblemPersistencies_.empty()) {
                    for (auto &subPersistency: subProblemPersistencies_) {
                        IndexEdgeSet const subEdges1 = subPersistency.edges1();
                        returnEdges1.insert(subEdges1.begin(), subEdges1.end());
                    }
                }

                return returnEdges1;
            }

            template<class T>
            inline
            typename Persistency<T>::PartitionSubsets const &
            Persistency<T>::subProblems() const {
                return subProblems_;
            }

            template<class T>
            bool
            Persistency<T>::tryFindNegativeNode(const std::set<size_t> &component, size_t *nextNode) {

                // check edges
                for (const size_t c: component) {
                    // edge (c, i)
                    for (size_t i = 0; i < c; ++i) {
                        if (component.find(i) != component.end()) {
                            // (c, i) is in component
                            continue;
                        }

                        if (problem_.costOfEdge({c, i}) < 0) {
                            *nextNode = i;
                            return false;
                        }
                    }

                    // edges (i, c)
                    for (size_t i = c + 1; i < problem_.numberOfElements(); ++i) {
                        if (component.find(i) != component.end()) {
                            continue;
                        }

                        if (problem_.costOfEdge({i, c}) < 0) {
                            *nextNode = i;
                            return false;
                        }
                    }
                }

                // iterate over all c's in component
                for (const size_t c: component) {
                    // all (c, i, j)
                    for (size_t i = 0; i < c; ++i) {

                        for (size_t j = 0; j < i; ++j) {
                            if (component.find(i) != component.end() && component.find(j) != component.end()) {
                                // (c,i,j) fully contained in component
                                continue;
                            }
                            if (problem_.costOfTriple({c, i, j}) < 0) {

                                if (component.find(i) != component.end()) {
                                    // j is not in component
                                    *nextNode = j;
                                } else if (component.find(j) != component.end()) {
                                    *nextNode = i;
                                } else {
                                    // break tie
                                    *nextNode = i;
                                }
                                return false;
                            }
                        }
                    }

                    // all (i, c, j)
                    for (size_t i = c + 1; i < problem_.numberOfElements(); ++i) {
                        for (size_t j = 0; j < c; ++j) {
                            if (component.find(i) != component.end() && component.find(j) != component.end()) {
                                // (c,i,j) fully contained in component
                                continue;
                            }
                            if (problem_.costOfTriple({i, c, j}) < 0) {

                                if (component.find(i) != component.end()) {
                                    // j is not in component
                                    *nextNode = j;
                                } else if (component.find(j) != component.end()) {
                                    *nextNode = i;
                                } else {
                                    // break tie
                                    *nextNode = i;
                                }
                                return false;
                            }
                        }
                    }

                    // all (i, j, c)
                    for (size_t i = c + 2; i < problem_.numberOfElements(); ++i) {
                        for (size_t j = c + 1; j < i; ++j) {
                            if (component.find(i) != component.end() && component.find(j) != component.end()) {
                                // (c,i,j) fully contained in component
                                continue;
                            }
                            if (problem_.costOfTriple({i, j, c}) < 0) {

                                if (component.find(i) != component.end()) {
                                    // j is not in component
                                    *nextNode = j;
                                } else if (component.find(j) != component.end()) {
                                    *nextNode = i;
                                } else {
                                    // break tie
                                    *nextNode = i;
                                }
                                return false;
                            }
                        }
                    }
                }

                return true;
            }

            // returns true if independent sets have been found
            template<class T>
            bool
            Persistency<T>::findIndependentSubProblems() {
                auto start = std::chrono::system_clock::now();

                bool foundIndependentSubProblems = false;

                // check if there are already sub persistencies
                std::set<std::set<size_t >> subProblems = {};

                if (subProblemPersistencies_.empty()) {
                    std::vector<size_t> points;
                    for (int i = 0; i < problem_.numberOfElements(); i++) {
                        points.push_back(i);
                    }

                    std::set<std::set<size_t>> partitions;

                    while (!points.empty()) {
                        std::set<size_t> subset;
                        const size_t seed = points[0];
                        subset.insert(seed);

                        size_t nextNode;
                        while (!tryFindNegativeNode(subset, &nextNode)) {
                            subset.insert(nextNode);
                        }

                        partitions.insert(subset);
                        subProblems.insert(subset);
                        for (size_t node: subset) {
                            if (std::find(points.begin(), points.end(), node) != points.end()) {
                                points.erase(std::remove(points.begin(), points.end(), node), points.end());
                            }
                        }
                    }

                    // adjust zero edges
                    for (auto const &subsetOne: subProblems) {
                        for (auto const &subsetTwo: subProblems) {
                            if (subsetOne == subsetTwo) continue;

                            for (auto const &i: subsetOne) {
                                for (auto const &j: subsetTwo) {
                                    if (i > j) {
                                        edges0_.insert({i, j});
                                    } else {
                                        edges0_.insert({j, i});
                                    }
                                }
                            }
                        }
                    }

                    if (subProblems.size() != 1) {
                        subProblemPersistencies_.clear();

                        // create sub persistencies
                        for (auto &subProblem: subProblems) {
                            std::vector<size_t> subsetArray(subProblem.size());
                            std::copy(subProblem.begin(), subProblem.end(), subsetArray.begin());
                            std::sort(subsetArray.begin(), subsetArray.end(), std::less<>());

                            ProblemType problem = problem_.orderedProjectTo(subsetArray);
                            // initialize new persistency class with respective nodeNames
                            Persistency<Rational> subPersistency(problem, subsetArray);
                            subProblemPersistencies_.push_back(subPersistency);
                        }


                        // more than one sub problem found
                        foundIndependentSubProblems = true;
                    } else {
                        // did not find any sub problems:
                        foundIndependentSubProblems = false;
                    }
                } else {
                    for (auto &subPersistency: subProblemPersistencies_) {
                        bool foundSubPersistencySubProblems = subPersistency.findIndependentSubProblems();
                        foundIndependentSubProblems = foundIndependentSubProblems || foundSubPersistencySubProblems;
                    }
                }

                auto end = std::chrono::system_clock::now();

                findIndependentProblemsDuration_ += end - start;

                return foundIndependentSubProblems;
            }


            template<class T>
            void Persistency<T>::edgeSubgraphCriterion() {
                auto start = std::chrono::system_clock::now();

                if (subProblemPersistencies_.empty()) {
                    IndexEdge mergeEdge;

                    while (checkEdgeSubgraph(&mergeEdge)) {
                        size_t s = mergeEdge[0] > mergeEdge[1] ? mergeEdge[0] : mergeEdge[1];
                        size_t t = mergeEdge[0] > mergeEdge[1] ? mergeEdge[1] : mergeEdge[0];
                        problem_.contract(s, t);
                    }
                } else {
                    for (auto &subPersistency: subProblemPersistencies_) {
                        subPersistency.edgeSubgraphCriterion();
                    }
                }


                auto end = std::chrono::system_clock::now();
                edgeSubgraphDuration_ += end - start;
            }

            template<class T>
            bool Persistency<T>::checkEdgeSubgraph(Persistency::IndexEdge *mergeEdge) {

                // iterate over all edges
                for (size_t j = 0; j < problem_.numberOfElements(); ++j) {
                    for (size_t k = 0; k < j; ++k) {
                        // Check criterion for edge (j, k)

                        // left hand side
                        const Rational edgeCost = problem_.costOfEdge({j, k});

                        // calculate right hand side
                        Rational rhs = 0;
                        // EDGE PART
                        // all edges with {j, l} with l neq k
                        for (size_t l = 0; l < j; ++l) {
                            if (l != k) {
                                const Rational cost = problem_.costOfEdge({j, l});
                                if (cost <= 0) {
                                    rhs += cost;
                                }

                                // triplet contribs with edge (j, l)
                                // (j, l, m)
                                for (size_t m = 0; m < l; ++m) {
                                    const Rational tripleCost = problem_.costOfTriple({j, l, m});
                                    if (tripleCost <= 0) {
                                        rhs += tripleCost / 2;
                                    }
                                }

                                // (j, m, l)
                                for (size_t m = l + 1; m < j; ++m) {
                                    const Rational tripleCost = problem_.costOfTriple({j, m, l});
                                    if (tripleCost <= 0) {
                                        rhs += tripleCost / 2;
                                    }
                                }

                                // (m, j, l)
                                for (size_t m = j + 1; m < problem_.numberOfElements(); ++m) {
                                    const Rational tripleCost = problem_.costOfTriple({m, j, l});
                                    if (tripleCost <= 0) {
                                        rhs += tripleCost / 2;
                                    }
                                }
                            }
                        }

                        // all edges with {l, j} with (l != k)
                        for (size_t l = j + 1; l < problem_.numberOfElements(); ++l) {
                            if (l != k) {
                                const Rational cost = problem_.costOfEdge({l, j});
                                if (cost <= 0) {
                                    rhs += cost;
                                }

                                // triple contributions with edge (l, j)
                                // (l, j, m)
                                for (size_t m = 0; m < j; ++m) {
                                    const Rational tripleCost = problem_.costOfTriple({l, j, m});
                                    if (tripleCost <= 0) {
                                        rhs += tripleCost / 2;
                                    }
                                }

                                // (l, m, j)
                                for (size_t m = j + 1; m < l; ++m) {
                                    const Rational tripleCost = problem_.costOfTriple({l, m, j});
                                    if (tripleCost <= 0) {
                                        rhs += tripleCost / 2;
                                    }
                                }

                                // (m, l, j)
                                for (size_t m = l + 1; m < problem_.numberOfElements(); ++m) {
                                    const Rational tripleCost = problem_.costOfTriple({m, l, j});
                                    if (tripleCost <= 0) {
                                        rhs += tripleCost / 2;
                                    }
                                }
                            }
                        }

                        // all edges with {k, l} with l \neq j
                        for (size_t l = 0; l < k; ++l) {
                            if (l != j) {
                                const Rational cost = problem_.costOfEdge({k, l});
                                if (cost <= 0) {
                                    rhs += cost;
                                }

                                // triple contributions
                                // (k, l, m)
                                for (size_t m = 0; m < l; ++m) {
                                    const Rational tripleCost = problem_.costOfTriple({k, l, m});
                                    if (tripleCost <= 0) {
                                        rhs += tripleCost / 2;
                                    }
                                }

                                // (k, m, l)
                                for (size_t m = l + 1; m < k; ++m) {
                                    const Rational tripleCost = problem_.costOfTriple({k, m, l});
                                    if (tripleCost <= 0) {
                                        rhs += tripleCost / 2;
                                    }
                                }

                                // (m, k, l)
                                for (size_t m = k + 1; m < problem_.numberOfElements(); ++m) {
                                    const Rational tripleCost = problem_.costOfTriple({m, k, l});
                                    if (tripleCost <= 0) {
                                        rhs += tripleCost / 2;
                                    }
                                }
                            }
                        }

                        // all edges with {l, k} with l \neq j
                        for (size_t l = k + 1; l < problem_.numberOfElements(); ++l) {
                            if (l != j) {
                                const Rational cost = problem_.costOfEdge({j, l});
                                if (cost <= 0) {
                                    rhs += cost;
                                }

                                // triple contributions
                                // (l, k, m)
                                for (size_t m = 0; m < k; ++m) {
                                    const Rational tripleCost = problem_.costOfTriple({l, k, m});
                                    if (tripleCost <= 0) {
                                        rhs += tripleCost / 2;
                                    }
                                }

                                // (l, m, k)
                                for (size_t m = k + 1; m < l; ++m) {
                                    const Rational tripleCost = problem_.costOfTriple({l, m, k});
                                    if (tripleCost <= 0) {
                                        rhs += tripleCost / 2;
                                    }
                                }

                                // (m, l, k)
                                for (size_t m = l + 1; m < problem_.numberOfElements(); ++m) {
                                    const Rational tripleCost = problem_.costOfTriple({m, l, k});
                                    if (tripleCost <= 0) {
                                        rhs += tripleCost / 2;
                                    }
                                }
                            }
                        }

                        if (edgeCost <= rhs) {
                            *mergeEdge = {j, k};
                            return true;
                        }


                    }
                }


                return false;
            }

            template<class T>
            void Persistency<T>::tripletSubgraphCriterion() {

                auto start = std::chrono::system_clock::now();

                if (subProblemPersistencies_.empty()) {
                    IndexEdge edge;

                    while (checkTripletSubgraph(&edge)) {
                        const size_t s = edge[0] > edge[1] ? edge[0] : edge[1];
                        const size_t t = edge[0] > edge[1] ? edge[1] : edge[0];
                        problem_.contract(s, t);
                    }
                } else {
                    for (auto &subPersistency: subProblemPersistencies_) {
                        subPersistency.tripletSubgraphCriterion();
                    }
                }

                auto end = std::chrono::system_clock::now();
                tripletSubgraphDuration_ += end - start;
            }


            template<class T>
            bool Persistency<T>::checkTripletSubgraph(Persistency::IndexTriple triple, Persistency::IndexEdge edge) {

                if (problem_.costOfEdge({triple[0], triple[1]}) + problem_.costOfEdge({triple[0], triple[2]}) > 0 ||
                    problem_.costOfEdge({triple[0], triple[1]}) + problem_.costOfEdge({triple[1], triple[2]}) > 0 ||
                    problem_.costOfEdge({triple[1], triple[2]}) + problem_.costOfEdge({triple[0], triple[2]}) > 0 ||
                    problem_.costOfTriple(triple) / 2 + problem_.costOfEdge({triple[0], triple[1]}) +
                    problem_.costOfEdge({triple[0], triple[2]}) + problem_.costOfEdge({triple[1], triple[2]}) > 0) {
                    return false;
                }

                const size_t u = edge[0];
                const size_t w = edge[1];
                size_t v;


                for (size_t vPrime: triple) {
                    if (vPrime != u && vPrime != w) {
                        v = vPrime;
                        break;
                    }
                }

                IndexEdge uv;
                if (u > v) {
                    uv = {u, v};
                } else {
                    uv = {v, u};
                }

                IndexEdge vw;
                if (v > w) {
                    vw = {v, w};
                } else {
                    vw = {w, v};
                }

                const Rational firstLhs =
                        problem_.costOfTriple(triple) + problem_.costOfEdge(uv) + problem_.costOfEdge(edge);
                const Rational secondLhs =
                        problem_.costOfTriple(triple) + problem_.costOfEdge(vw) + problem_.costOfEdge(edge);

                Rational rhs = 0;
                // calculate RHS
                // EDGE CONTRIBUTIONS
                // edges which contain u
                // (j, u)
                for (size_t j = u + 1; j < problem_.numberOfElements(); ++j) {
                    if (j != v && j != w) {
                        const Rational edgeCost = problem_.costOfEdge({j, u});
                        if (edgeCost < 0) {
                            rhs += edgeCost;
                        }

                        // contributions of triplets T_{\delta(j, u)}
                        // contributions (j, u, k)
                        for (size_t k = 0; k < u; ++k) {
                            const Rational tripleCost = problem_.costOfTriple({j, u, k});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }

                        // (j, k, u)
                        for (size_t k = u + 1; k < j; ++k) {
                            const Rational tripleCost = problem_.costOfTriple({j, k, u});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }

                        // (k, j, u)
                        for (size_t k = j + 1; k < problem_.numberOfElements(); ++k) {
                            const Rational tripleCost = problem_.costOfTriple({k, j, u});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }
                    }
                }

                // (u, j)
                for (size_t j = 0; j < u; ++j) {
                    if (j != v && j != w) {
                        const Rational edgeCost = problem_.costOfEdge({u, j});
                        if (edgeCost < 0) {
                            rhs += edgeCost;
                        }

                        // triplet contributions
                        // (u, j, k)
                        for (size_t k = 0; k < j; ++k) {
                            const Rational tripleCost = problem_.costOfTriple({u, j, k});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }

                        // (u, k, j)
                        for (size_t k = j + 1; k < u; ++k) {
                            const Rational tripleCost = problem_.costOfTriple({u, k, j});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }

                        // (k, u, j)
                        for (size_t k = u + 1; k < problem_.numberOfElements(); ++k) {
                            const Rational tripleCost = problem_.costOfTriple({k, u, j});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }
                    }
                }

                // (j, v)
                for (size_t j = v + 1; j < problem_.numberOfElements(); ++j) {
                    if (j != u && j != w) {
                        const Rational edgeCost = problem_.costOfEdge({j, v});
                        if (edgeCost < 0) {
                            rhs += edgeCost;
                        }

                        // triplet edgeCost contributions
                        // (j, v, k)
                        for (size_t k = 0; k < v; ++k) {
                            const Rational tripleCost = problem_.costOfTriple({j, v, k});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }

                        // (j, k, v)
                        for (size_t k = v + 1; k < j; ++k) {
                            const Rational tripleCost = problem_.costOfTriple({j, k, v});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }

                        // (k, j, v)
                        for (size_t k = j + 1; k < problem_.numberOfElements(); ++k) {
                            const Rational tripleCost = problem_.costOfTriple({k, j, v});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }
                    }
                }

                // (v, j)
                for (size_t j = 0; j < v; ++j) {
                    if (j != u && j != w) {
                        const Rational edgeCost = problem_.costOfEdge({v, j});
                        if (edgeCost < 0) {
                            rhs += edgeCost;
                        }

                        // triplet edgeCost contributions
                        // (v, j, k)
                        for (size_t k = 0; k < j; ++k) {
                            const Rational tripleCost = problem_.costOfTriple({v, j, k});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }

                        // (v, k, j)
                        for (size_t k = j + 1; k < v; ++k) {
                            const Rational tripleCost = problem_.costOfTriple({v, k, j});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }

                        // (k, v, j)
                        for (size_t k = v + 1; k < problem_.numberOfElements(); ++k) {
                            const Rational tripleCost = problem_.costOfTriple({k, v, j});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }
                    }
                }

                // (j, w)
                for (size_t j = w + 1; j < problem_.numberOfElements(); ++j) {
                    if (j != u && j != v) {
                        const Rational edgeCost = problem_.costOfEdge({j, w});
                        if (edgeCost < 0) {
                            rhs += edgeCost;
                        }

                        // triplet edgeCost contributions
                        // (j, w, k)
                        for (size_t k = 0; k < w; ++k) {
                            const Rational tripleCost = problem_.costOfTriple({j, w, k});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }

                        // (j, k, w)
                        for (size_t k = w + 1; k < j; ++k) {
                            const Rational tripleCost = problem_.costOfTriple({j, k, w});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }

                        // (k, j, w)
                        for (size_t k = j + 1; k < problem_.numberOfElements(); ++k) {
                            const Rational tripleCost = problem_.costOfTriple({k, j, w});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }
                    }
                }

                // (w, j)
                for (size_t j = 0; j < w; ++j) {
                    if (j != u && j != v) {
                        const Rational cost = problem_.costOfEdge({w, j});
                        if (cost < 0) {
                            rhs += cost;
                        }

                        // triplet edgeCost contributions
                        // (w, j, k)
                        for (size_t k = 0; k < j; ++k) {
                            const Rational tripleCost = problem_.costOfTriple({w, j, k});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }

                        // (w, k, j)
                        for (size_t k = j + 1; k < w; ++k) {
                            const Rational tripleCost = problem_.costOfTriple({w, k, j});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }

                        // (k, w, j)
                        for (size_t k = w + 1; k < problem_.numberOfElements(); ++k) {
                            const Rational tripleCost = problem_.costOfTriple({k, w, j});
                            if (tripleCost < 0) {
                                rhs += tripleCost / 2;
                            }
                        }
                    }
                }

                return firstLhs <= rhs && secondLhs <= rhs;
            }

            template<class T>
            bool Persistency<T>::checkTripletSubgraph(Persistency::IndexEdge *edge) {
                for (size_t i = 0; i < problem_.numberOfElements(); ++i) {
                    for (size_t j = 0; j < i; ++j) {
                        for (size_t k = 0; k < j; ++k) {
                            if (edges1_.find({i, j}) == edges1_.end() && checkTripletSubgraph({i, j, k}, {i, j})) {
                                *edge = {i, j};
                                return true;
                            }

                            if (edges1_.find({i, k}) == edges1_.end() && checkTripletSubgraph({i, j, k}, {i, k})) {
                                *edge = {i, k};
                                return true;
                            }

                            if (edges1_.find({j, k}) == edges1_.end() && checkTripletSubgraph({i, j, k}, {i, k})) {
                                *edge = {j, k};
                                return true;
                            }
                        }
                    }
                }
                return false;
            }


            template<class T>
            void Persistency<T>::tripletCutCriterion() {
                auto start = std::chrono::system_clock::now();

                if (subProblemPersistencies_.empty()) {

                    NegativePartCostAdaptor onlyNegativeCostProblem(problem_);

                    for (size_t i = 0; i < problem_.numberOfElements(); ++i) {
                        for (size_t j = 0; j < i; ++j) {
                            for (size_t k = 0; k < j; ++k) {


                                // test condition for i being separated from j and k
                                Rational firstRHS = getMinimalSeparationCost(
                                        onlyNegativeCostProblem,
                                        std::vector<size_t>({i}),
                                        std::vector<size_t>({j, k})
                                );

                                Rational firstLHS = positivePart(problem_.costOfTriple({i, j, k})) +
                                                    positivePart(problem_.costOfEdge({i, j})) +
                                                    positivePart(problem_.costOfEdge({i, k}));

                                if (firstLHS >= firstRHS) {
                                    triples0_.insert({i, j, k});
                                    continue;
                                }

                                // test condition for j being separated from i and j
                                Rational secondRHS = getMinimalSeparationCost(
                                        onlyNegativeCostProblem,
                                        std::vector<size_t>({j}),
                                        std::vector<size_t>({i, k})
                                );

                                Rational secondLHS = positivePart(problem_.costOfTriple({i, j, k}))
                                                     + positivePart(problem_.costOfEdge({i, j}))
                                                     + positivePart(problem_.costOfEdge({j, k}));

                                if (secondLHS >= secondRHS) {
                                    triples0_.insert({i, j, k});
                                    continue;
                                }

                                // test condition for k being separated from i and j
                                Rational thirdRHS = getMinimalSeparationCost(
                                        onlyNegativeCostProblem,
                                        std::vector<size_t>({k}),
                                        std::vector<size_t>({i, j})
                                );

                                Rational thirdLHS = positivePart(problem_.costOfTriple({i, j, k}))
                                                    + positivePart(problem_.costOfEdge({i, k}))
                                                    + positivePart(problem_.costOfEdge({j, k}));

                                if (thirdLHS >= thirdRHS) {
                                    triples0_.insert({i, j, k});
                                    continue;
                                }

                            }
                        }
                    }
                } else {
                    for (auto &subPersistency: subProblemPersistencies_) {
                        subPersistency.tripletCutCriterion();
                    }
                }

                auto end = std::chrono::system_clock::now();

                tripletCutDuration_ += end - start;
            }


            template<class T>
            template<class PROBLEM_VIEW>
            typename Persistency<T>::Rational Persistency<T>::getMinimalSeparationCost(
                    const PROBLEM_VIEW &problem,
                    const std::vector<size_t> &inSet,
                    const std::vector<size_t> &outOfSet
            ) {
                typedef QPBOCutProblemAdaptor<T> ProblemAdaptor;
                typedef qpbo::MaxFlowSolver MaxFlowSolver;

                ProblemAdaptor adaptor(problem, inSet, outOfSet);

                Rational minimalValue;

                andres::graph::qpbo::solveSubmodularQPBO<Rational, MaxFlowSolver::BOOST_PUSH_RELABEL>(
                        adaptor, minimalValue);


                return minimalValue;
            }

            template<class T>
            inline
            typename Persistency<T>::Rational
            Persistency<T>::negativePart(const Rational value) {
                if (value < 0) {
                    return -value;
                }
                return 0;
            }

            template<class T>
            inline
            typename Persistency<T>::Rational
            Persistency<T>::positivePart(const Rational value) {
                if (value > 0) {
                    return value;
                }
                return 0;
            }

            template<class T>
            bool Persistency<T>::tripletEdgeJoinCriterion(Persistency::IndexEdge *mergeEdge) {
                for (size_t j = 0; j < problem_.numberOfElements(); ++j) {
                    for (size_t k = 0; k < j; ++k) {
                        for (size_t l = 0; l < k; ++l) {
                            if (edges1_.find({j, k}) == edges1_.end() && tripletEdgeJoinCriterion({j, k, l}, {j, k})) {
                                *mergeEdge = {j, k};
                                return true;
                            }
                            if (edges1_.find({j, l}) == edges1_.end() && tripletEdgeJoinCriterion({j, k, l}, {j, l})) {
                                *mergeEdge = {j, l};
                                return true;
                            }
                            if (edges1_.find({k, l}) == edges1_.end() && tripletEdgeJoinCriterion({j, k, l}, {k, l})) {
                                *mergeEdge = {k, l};
                                return true;
                            }
                        }
                    }
                }
                return false;
            }

            template<class T>
            inline
            bool
            Persistency<T>::tripletEdgeJoinCriterion(IndexTriple triple, IndexEdge edge) {
                const size_t i = edge[0];
                const size_t k = edge[1];
                size_t j;

                for (size_t index: triple) {
                    if (index != i && index != k) {
                        j = index;
                    }
                }

                IndexEdge ijEdge;
                IndexEdge jkEdge;

                if (i > j) {
                    ijEdge = {i, j};
                } else {
                    ijEdge = {j, i};
                }

                if (j > k) {
                    jkEdge = {j, k};
                } else {
                    jkEdge = {k, j};
                }

                const AbsoluteCostProblemAdaptor absoluteCostProblem = AbsoluteCostProblemAdaptor(problem_);


                // third condition
                Rational lhsThird = problem_.costOfTriple(triple)
                                    + problem_.costOfEdge({triple[0], triple[1]})
                                    + problem_.costOfEdge({triple[0], triple[2]})
                                    + problem_.costOfEdge({triple[1], triple[2]});

                if (lhsThird > 0) return false;

                Rational rhsThird = calculateTripletEdgeJoinCriterionThirdRightHandSide(triple);


                if (lhsThird > rhsThird) {
                    return false;
                }

                // first condition
                // set U which separates edge[0]=i from edge[1]=k and j
                Rational lhsFirst = calculateTripletEdgeJoinCriterionFirstOrSecondLeftHandSide(ijEdge, edge);
                Rational rhsFirst = getMinimalSeparationCost(
                        absoluteCostProblem,
                        std::vector<size_t>({i}),
                        std::vector<size_t>({j, k})
                );

                if (lhsFirst < rhsFirst) {
                    return false;
                }

                // second condition
                // set W which separates edge[1]=k from edge[0]=i and j
                Rational lhsSecond = calculateTripletEdgeJoinCriterionFirstOrSecondLeftHandSide(jkEdge, edge);
                Rational rhsSecond = getMinimalSeparationCost(
                        absoluteCostProblem,
                        std::vector<size_t>({k}),
                        std::vector<size_t>({i, j})
                );

                if (lhsSecond < rhsSecond) {
                    return false;
                }


                return true;
            }

            template<class T>
            inline
            void
            Persistency<T>::tripletEdgeJoinCriterion() {

                auto start = std::chrono::system_clock::now();

                if (subProblemPersistencies_.empty()) {
                    IndexEdge edge;

                    while (tripletEdgeJoinCriterion(&edge)) {
                        const size_t s = edge[0] > edge[1] ? edge[0] : edge[1];
                        const size_t t = edge[0] > edge[1] ? edge[1] : edge[0];
                        problem_.contract(s, t);
                    }
                } else {
                    for (auto &subPersistency: subProblemPersistencies_) {
                        subPersistency.tripletEdgeJoinCriterion();
                    }
                }


                auto end = std::chrono::system_clock::now();
                tripletEdgeJoinDuration_ += end - start;
            }

            template<class T>
            inline
            typename Persistency<T>::Rational
            Persistency<T>::calculateTripletEdgeJoinCriterionThirdRightHandSide(IndexTriple triple) {

                Rational rhsThird = 0;

                // all with only triple[0]
                // all triples (triple[0], j, k)
                // all edges (triple[0], j)
                for (size_t j = 0; j < triple[0]; ++j) {
                    if (j == triple[1] || j == triple[2]) continue;
                    // edge cost
                    Rational edgeCost = problem_.costOfEdge({triple[0], j});
                    if (edgeCost < 0) {
                        rhsThird += edgeCost;
                    }

                    for (size_t k1 = 0; k1 < j; ++k1) {
                        if (k1 == triple[1] || k1 == triple[2]) continue;

                        Rational tripletCost = problem_.costOfTriple({triple[0], j, k1});
                        if (tripletCost < 0) {
                            rhsThird += tripletCost;
                        }
                    }
                }

                // all triples (j, triple[0], k)
                // all edges (j, triple[0])
                for (size_t j = triple[0] + 1; j < problem_.numberOfElements(); ++j) {
                    if (j == triple[1] || j == triple[2]) continue;

                    Rational edgeCost = problem_.costOfEdge({j, triple[0]});
                    if (edgeCost < 0) {
                        rhsThird += edgeCost;
                    }

                    for (size_t k1 = 0; k1 < triple[0]; ++k1) {
                        if (k1 == triple[1] || k1 == triple[2]) continue;

                        Rational tripleCost = problem_.costOfTriple({j, triple[0], k1});
                        if (tripleCost < 0) {
                            rhsThird += tripleCost;
                        }
                    }
                }

                // all triples (j, k, triple[0])
                for (size_t j = triple[0] + 2; j < problem_.numberOfElements(); ++j) {
                    if (j == triple[1] || j == triple[2]) continue;


                    for (size_t k1 = triple[0] + 1; k1 < j; ++k1) {
                        if (k1 == triple[1] || k1 == triple[2]) continue;

                        Rational tripleCost = problem_.costOfTriple({j, k1, triple[0]});
                        if (tripleCost < 0) {
                            rhsThird += tripleCost;
                        }
                    }
                }

                // all with only triple[1]
                // all triples (triple[1], j, k)
                // all edges (triple[1], j)
                for (size_t j = 0; j < triple[1]; ++j) {
                    if (j == triple[0] || j == triple[2]) continue;
                    // edge cost
                    Rational edgeCost = problem_.costOfEdge({triple[1], j});
                    if (edgeCost < 0) {
                        rhsThird += edgeCost;
                    }

                    for (size_t k1 = 0; k1 < j; ++k1) {
                        if (k1 == triple[0] || k1 == triple[2]) continue;

                        Rational tripletCost = problem_.costOfTriple({triple[1], j, k1});
                        if (tripletCost < 0) {
                            rhsThird += tripletCost;
                        }
                    }
                }

                // all triples (j, triple[1], k)
                // all edges (j, triple[1])
                for (size_t j = triple[1] + 1; j < problem_.numberOfElements(); ++j) {
                    if (j == triple[0] || j == triple[2]) continue;

                    Rational edgeCost = problem_.costOfEdge({j, triple[1]});
                    if (edgeCost < 0) {
                        rhsThird += edgeCost;
                    }

                    for (size_t k1 = 0; k1 < triple[1]; ++k1) {
                        if (k1 == triple[0] || k1 == triple[2]) continue;

                        Rational tripleCost = problem_.costOfTriple({j, triple[1], k1});
                        if (tripleCost < 0) {
                            rhsThird += tripleCost;
                        }
                    }
                }

                // all triples (j, k, triple[1])
                for (size_t j = triple[1] + 2; j < problem_.numberOfElements(); ++j) {
                    if (j == triple[0] || j == triple[2]) continue;

                    for (size_t k1 = triple[1] + 1; k1 < j; ++k1) {
                        if (k1 == triple[0] || k1 == triple[2]) continue;

                        Rational tripleCost = problem_.costOfTriple({j, k1, triple[1]});
                        if (tripleCost < 0) {
                            rhsThird += tripleCost;
                        }
                    }
                }

                // all with only triple[2]
                // all triples (triple[2], j, k)
                // all edges (triple[2], j)
                for (size_t j = 0; j < triple[2]; ++j) {
                    if (j == triple[0] || j == triple[1]) continue;
                    // edge cost
                    Rational edgeCost = problem_.costOfEdge({triple[2], j});
                    if (edgeCost < 0) {
                        rhsThird += edgeCost;
                    }

                    for (size_t k1 = 0; k1 < j; ++k1) {
                        if (k1 == triple[0] || k1 == triple[1]) continue;

                        Rational tripletCost = problem_.costOfTriple({triple[2], j, k1});
                        if (tripletCost < 0) {
                            rhsThird += tripletCost;
                        }
                    }
                }

                // all triples (j, triple[2], k)
                // all edges (j, triple[2])
                for (size_t j = triple[2] + 1; j < problem_.numberOfElements(); ++j) {
                    if (j == triple[1] || j == triple[2]) continue;

                    Rational edgeCost = problem_.costOfEdge({j, triple[2]});
                    if (edgeCost < 0) {
                        rhsThird += edgeCost;
                    }

                    for (size_t k1 = 0; k1 < triple[2]; ++k1) {
                        if (k1 == triple[0] || k1 == triple[1]) continue;

                        Rational tripleCost = problem_.costOfTriple({j, triple[2], k1});
                        if (tripleCost < 0) {
                            rhsThird += tripleCost;
                        }
                    }
                }

                // all triples (j, k, triple[2])
                for (size_t j = triple[2] + 2; j < problem_.numberOfElements(); ++j) {
                    if (j == triple[0] || j == triple[1]) continue;

                    for (size_t k1 = triple[2] + 1; k1 < j; ++k1) {
                        if (k1 == triple[0] || k1 == triple[1]) continue;

                        Rational tripleCost = problem_.costOfTriple({j, k1, triple[2]});
                        if (tripleCost < 0) {
                            rhsThird += tripleCost;
                        }
                    }
                }

                return rhsThird;
            }

            template<class T>
            inline
            typename Persistency<T>::Rational
            Persistency<T>::calculateTripletEdgeJoinCriterionFirstOrSecondLeftHandSide(const IndexEdge firstEdge,
                                                                                       const IndexEdge secondEdge) {
                Rational leftHandSide = 2 * negativePart(problem_.costOfEdge(firstEdge))
                                        + 2 * negativePart(problem_.costOfEdge(secondEdge));


                // all triples with (firstEdge[0], firstEdge[1], i)
                for (size_t i = 0; i < firstEdge[1]; ++i) {
                    leftHandSide += negativePart(problem_.costOfTriple({firstEdge[0], firstEdge[1], i}));
                }

                // all triples with (firstEdge[0], i, firstEdge[1])
                for (size_t i = firstEdge[1] + 1; i < firstEdge[0]; ++i) {
                    leftHandSide += negativePart(problem_.costOfTriple({firstEdge[0], i, firstEdge[1]}));
                }

                // all triples with (i, firstEdge[0], firstEdge[1])
                for (size_t i = firstEdge[0] + 1; i < problem_.numberOfElements(); ++i) {
                    leftHandSide += negativePart(problem_.costOfTriple({i, firstEdge[0], firstEdge[1]}));
                }

                // all triples with (secondEdge[0], secondEdge[1], i)
                for (size_t i = 0; i < secondEdge[1]; ++i) {
                    leftHandSide += negativePart(problem_.costOfTriple({secondEdge[0], secondEdge[1], i}));
                }

                // all triples with (secondEdge[0], i, secondEdge[1])
                for (size_t i = secondEdge[1] + 1; i < secondEdge[0]; ++i) {
                    leftHandSide += negativePart(problem_.costOfTriple({secondEdge[0], i, secondEdge[1]}));
                }

                // all triples with (i, secondEdge[0], secondEdge[1])
                for (size_t i = secondEdge[0] + 1; i < problem_.numberOfElements(); ++i) {
                    leftHandSide += negativePart(problem_.costOfTriple({i, secondEdge[0], secondEdge[1]}));
                }

                return leftHandSide;
            }

            template<class T>
            void Persistency<T>::edgeCutCriterion() {
                auto start = std::chrono::system_clock::now();

                if (subProblemPersistencies_.empty()) {
                    NegativePartCostAdaptor negativeCostProblem(problem_);


                    for (size_t i = 0; i < problem_.numberOfElements(); ++i) {
                        for (size_t j = 0; j < i; ++j) {
                            Rational lhs = positivePart(problem_.costOfEdge({i, j}));

                            Rational rhs = getMinimalSeparationCost(
                                    negativeCostProblem,
                                    std::vector<size_t>({i}),
                                    std::vector<size_t>({j})
                            );

                            if (lhs >= rhs) {
                                edges0_.insert({i, j});
                            }
                        }
                    }
                } else {
                    for (auto &subPersistency: subProblemPersistencies_) {
                        subPersistency.edgeCutCriterion();
                    }
                }


                auto end = std::chrono::system_clock::now();

                edgeCutDuration_ += end - start;
            }

            template<class T>
            void Persistency<T>::edgeJoinCriterion() {
                auto start = std::chrono::system_clock::now();

                if (subProblemPersistencies_.empty()) {
                    IndexEdge edge;

                    while (edgeJoinCriterion(&edge)) {
                        problem_.contract(edge[0], edge[1]);
                    }

                } else {
                    for (auto &subPersistency: subProblemPersistencies_) {
                        subPersistency.edgeJoinCriterion();
                    }
                }

                auto end = std::chrono::system_clock::now();
                edgeJoinDuration_ += end - start;
            }

            template<class T>
            bool Persistency<T>::edgeJoinCriterion(IndexEdge *edge) {
                AbsoluteCostAdaptor absoluteCostProblem(problem_);


                for (size_t i = 0; i < problem_.numberOfElements(); ++i) {
                    for (size_t j = 0; j < i; ++j) {
                        Rational lhs = 2 * negativePart(problem_.costOfEdge({i, j}));

                        // LHS: all (i, j, r)
                        for (size_t r = 0; r < j; ++r) {
                            lhs += negativePart(problem_.costOfTriple({i, j, r}));
                        }
                        // LHS: all (i, r, j)
                        for (size_t r = j + 1; r < i; ++r) {
                            lhs += negativePart(problem_.costOfTriple({i, r, j}));
                        }
                        // LHS: all (r, i, j)
                        for (size_t r = i + 1; r < problem_.numberOfElements(); ++r) {
                            lhs += negativePart(problem_.costOfTriple({r, i, j}));
                        }

                        Rational rhs = getMinimalSeparationCost(
                                absoluteCostProblem,
                                std::vector<size_t>({i}),
                                std::vector<size_t>({j})
                        );

                        if (lhs >= rhs) {
                            *edge = {i, j};
                            return true;
                        }

                    }
                }

                return false;
            }

            template<class T>
            bool Persistency<T>::tripletJoinCriterion(IndexTriple *triple) {
                NegativePartCostAdaptor negativeCostProblem(problem_);


                for (size_t i = 0; i < problem_.numberOfElements(); ++i) {
                    for (size_t j = 0; j < i; ++j) {
                        for (size_t k = 0; k < j; ++k) {
                            const Rational firstLHS = calculateTripletJoinCriterionLeftHandSide({i, j, k}, {j, k});

                            const Rational firstRHS = getMinimalSeparationCost(
                                    negativeCostProblem,
                                    std::vector<size_t>({i}),
                                    std::vector<size_t>({j, k})
                            );

                            if (firstLHS >= firstRHS) {
                                *triple = {i, j, k};
                                return true;
                            }

                            const Rational secondLHS = calculateTripletJoinCriterionLeftHandSide({i, j, k}, {i, k});

                            const Rational secondRHS = getMinimalSeparationCost(
                                    negativeCostProblem,
                                    std::vector<size_t>({j}),
                                    std::vector<size_t>({i, k})
                            );

                            if (secondLHS >= secondRHS) {
                                *triple = {i, j, k};
                                return true;
                            }

                            const Rational thirdLHS = calculateTripletJoinCriterionLeftHandSide({i, j, k}, {i, j});

                            const Rational thirdRHS = getMinimalSeparationCost(
                                    negativeCostProblem,
                                    std::vector<size_t>({k}),
                                    std::vector<size_t>({i, j})
                            );

                            if (thirdLHS >= thirdRHS) {
                                *triple = {i, j, k};
                                return true;
                            }
                        }
                    }
                }

                return false;
            }

            template<class T>
            void Persistency<T>::tripletJoinCriterion() {
                auto start = std::chrono::system_clock::now();

                if (subProblemPersistencies_.empty()) {
                    IndexTriple triple;

                    while (tripletJoinCriterion(&triple)) {
                        problem_.contract(triple[0], triple[1], triple[2]);
                    }
                } else {
                    // do subPersistency criterion
                    for (auto &subPersistency: subProblemPersistencies_) {
                        subPersistency.tripletJoinCriterion();
                    }
                }

                auto end = std::chrono::system_clock::now();
                tripletJoinDuration_ += end - start;
            }

            template<class T>
            typename Persistency<T>::Rational
            Persistency<T>::calculateTripletJoinCriterionLeftHandSide(const IndexTriple triple, const IndexEdge edge) {
                IndexEdge ijEdge;
                IndexEdge ikEdge;
                // jkEdge = passed edge

                const size_t j = edge[0];
                const size_t k = edge[1];

                // find i
                size_t i;
                for (auto index: triple) {
                    if (index != j && index != k) {
                        i = index;
                    }
                }

                if (i > j) {
                    ijEdge = {i, j};
                } else {
                    ijEdge = {j, i};
                }

                if (i > k) {
                    ikEdge = {i, k};
                } else {
                    ikEdge = {k, i};
                }

                Rational lhs = 2 * negativePart(problem_.costOfTriple(triple))
                               + 2 * negativePart(problem_.costOfEdge(ijEdge))
                               + 2 * negativePart(problem_.costOfEdge(ikEdge))
                               + negativePart(problem_.costOfEdge(edge));

                std::initializer_list<Rational> possibleMinTerms = {
                        0,
                        problem_.costOfEdge({triple[0], triple[1]}),
                        problem_.costOfEdge({triple[0], triple[2]}),
                        problem_.costOfEdge({triple[1], triple[2]})
                };

                lhs += min(possibleMinTerms);

                for (size_t i = 0; i < problem_.numberOfElements(); ++i) {
                    for (size_t j = 0; j < i; ++j) {
                        lhs -= positivePart(problem_.costOfEdge({i, j}));

                        for (size_t k = 0; k < j; ++k) {
                            lhs -= positivePart(problem_.costOfTriple({i, j, k}));
                        }
                    }
                }

                return lhs;
            }

            template<class T>
            typename Persistency<T>::Duration
            const &Persistency<T>::totalDurationEdgeCutCriterion() const {
                return edgeCutDuration_;
            }

            template<class T>
            typename Persistency<T>::Duration
            const &Persistency<T>::totalDurationEdgeJoinCriterion() const {
                return edgeJoinDuration_;
            }

            template<class T>
            typename Persistency<T>::Duration
            const &Persistency<T>::totalDurationEdgeSubgraphCriterion() const {
                return edgeSubgraphDuration_;
            }

            template<class T>
            typename Persistency<T>::Duration
            const &Persistency<T>::totalDurationFindIndependentProblems() const {
                return findIndependentProblemsDuration_;
            }

            template<class T>
            typename Persistency<T>::Duration
            const &Persistency<T>::totalDurationTripletCutCriterion() const {
                return tripletCutDuration_;
            }

            template<class T>
            typename Persistency<T>::Duration
            const &Persistency<T>::totalDurationTripletEdgeJoinCriterion() const {
                return tripletEdgeJoinDuration_;
            }

            template<class T>
            typename Persistency<T>::Duration
            const &Persistency<T>::totalDurationTripletJoinCriterion() const {
                return tripletJoinDuration_;
            }

            template<class T>
            typename Persistency<T>::Duration
            const &Persistency<T>::totalDurationTripletSubgraphCriterion() const {
                return tripletSubgraphDuration_;
            }


            template<class T>
            inline
            typename Persistency<T>::Duration const &
            Persistency<T>::totalDurationSubsetJoinCriterion() const {
                return subsetJoinDuration_;
            }

            template<class T>
            inline
            void
            Persistency<T>::subsetJoinCriterion() {
                auto start = std::chrono::system_clock::now();

                if (subProblemPersistencies_.empty()) {
                    std::set<size_t> subsetToMerge;

                    while (subsetJoinCriterion(&subsetToMerge)) {
                        std::vector<size_t> subsetArray(subsetToMerge.size());
                        std::copy(subsetToMerge.begin(), subsetToMerge.end(), subsetArray.begin());
                        std::sort(subsetArray.begin(), subsetArray.end(), std::greater<>());

                        problem_.contract(subsetArray);
                    }
                } else {
                    for (auto &subPersistency: subProblemPersistencies_) {
                        subPersistency.subsetJoinCriterion();
                    }
                }

                auto end = std::chrono::system_clock::now();
                subsetJoinDuration_ = end - start;
            }

            template<class T>
            inline
            bool
            Persistency<T>::subsetJoinCriterion(std::set<size_t> *subset) {
                for (size_t i = 0; i < problem_.numberOfElements(); ++i) {
                    for (size_t j = 0; j < i; ++j) {
                        if (problem_.costOfEdge({i, j}) > 0) continue;
                        std::set<size_t> subsetProposal = findMaximalComponentWithOnlyNegativeWeights({i, j});

                        Rational rhs = subsetJoinCalculateRightHandSide(subsetProposal);
                        Rational lhs = subsetJoinCalculateLeftHandSide(subsetProposal);

                        // do merge
                        if (lhs <= rhs) {
                            *subset = subsetProposal;
                            return true;

                        }
                    }
                }

                return false;
            }


            template<class T>
            inline
            std::pair<bool, typename Persistency<T>::Rational>
            Persistency<T>::canMergeSubsetJoin(std::set<size_t> const &subset, size_t const &nextNode) {

                Rational offset = 0;

                for (auto &s1: subset) {
                    // check edge cost
                    IndexEdge edge = {s1, nextNode};
                    std::sort(edge.begin(), edge.end(), std::greater_equal<>());

                    Rational const edgeCost = problem_.costOfEdge(edge);

                    if (edgeCost > 0) return std::pair<bool, Rational>(false, -1.0);

                    offset += edgeCost;

                    for (auto &s2: subset) {
                        if (s2 != s1) {
                            IndexTriple triple = {s1, s2, nextNode};
                            std::sort(triple.begin(), triple.end(), std::greater_equal<>());

                            Rational const tripleCost = problem_.costOfTriple(triple);

                            if (tripleCost > 0) return std::pair<bool, Rational>(false, -1.0);

                            offset += tripleCost;
                        }
                    }
                }


                return std::pair<bool, Rational>(true, offset);
            }

            template<class T>
            inline
            std::set<size_t>
            Persistency<T>::findMaximalComponentWithOnlyNegativeWeights(size_t const seed) {

                std::vector<size_t> queue;
                for (size_t index = 0; index < problem_.numberOfElements(); ++index) {
                    if (index != seed) {
                        queue.push_back(index);
                    }
                }

                std::set<size_t> subset({seed});

                while (!queue.empty()) {
                    size_t nextNode = queue[0];

                    std::pair<bool, Rational> canMerge = canMergeSubsetJoin(subset, nextNode);

                    if (canMerge.first) {
                        subset.insert(nextNode);
                    }
                    queue.erase(std::remove(queue.begin(), queue.end(), nextNode), queue.end());
                }

                return subset;
            }

            template<class T>
            inline
            std::set<size_t>
            Persistency<T>::findMaximalComponentWithOnlyNegativeWeights(std::set<size_t> const &seedSet) {
                // assert seed Set is already negative component
                for (auto &i: seedSet) {
                    for (auto &j: seedSet) {
                        if (i == j) continue;
                        assert(problem_.costOfEdge({i, j}) <= 0);
                        for (auto &k: seedSet) {
                            if (k == i || k == j) continue;
                            assert(problem_.costOfTriple({i, j, k}) <= 0);
                        }
                    }
                }

                std::vector<size_t> queue;
                for (size_t index = 0; index < problem_.numberOfElements(); ++index) {
                    if (seedSet.find(index) == seedSet.end()) {
                        queue.push_back(index);
                    }
                }

                std::set<size_t> subset(seedSet);

                bool didMerge = true;

                while (didMerge) {
                    didMerge = false;
                    bool foundNodeToAdd = false;

                    std::vector<Rational> negativeOffsets(queue.size(), std::numeric_limits<Rational>().max());

                    // calc all negative offsets to merge possible nodes
                    for (size_t queueIndex = 0; queueIndex < queue.size(); ++queueIndex) {
                        size_t element = queue[queueIndex];

                        std::pair<bool, Rational> canMerge = canMergeSubsetJoin(subset, element);

                        if (canMerge.first) {
                            foundNodeToAdd = true;
                            negativeOffsets[queueIndex] = canMerge.second;
                        }
                    }

                    // find minimum of offsets if foundNodeToAdd = true
                    if (foundNodeToAdd) {
                        size_t minIndex;
                        Rational minOffset = std::numeric_limits<Rational>().max();

                        for (size_t queueIndex = 0; queueIndex < queue.size(); ++queueIndex) {
                            if (negativeOffsets[queueIndex] < minOffset) {
                                // replace minIndex
                                minIndex = queueIndex;
                                minOffset = negativeOffsets[queueIndex];
                            }
                        }

                        // determined most negative node
                        subset.insert(queue[minIndex]);

                        // remove from queue
                        queue.erase(queue.begin() + static_cast<long>(minIndex));
                        didMerge = true;
                    }
                }

                return subset;
            }

            template<class T>
            inline
            typename Persistency<T>::Rational
            Persistency<T>::subsetJoinCalculateLeftHandSide(const std::set<size_t> &subset) {
                ProblemType subsetProblem = problem_.unorderedProjectTo(subset);


                struct edge_t {
                    size_t first;
                    size_t second;
                };

                typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                        boost::no_property, boost::property<boost::edge_weight_t, Rational> > undirected_graph;
                typedef typename boost::property_map<undirected_graph, boost::edge_weight_t>::type weight_map_type;
                typedef typename boost::property_traits<weight_map_type>::value_type weight_type;

                size_t const numberOfNodes = subsetProblem.numberOfElements();
                size_t const numberOfEdges = numberOfNodes * (numberOfNodes - 1) / 2;

                edge_t edges[numberOfEdges];
                weight_type ws[numberOfEdges];

                size_t edgeIndex = 0;
                for (size_t i = 0; i < numberOfNodes; ++i) {
                    for (size_t j = 0; j < i; ++j) {
                        // impose e1 > e2
                        Rational edgeCost = -subsetProblem.costOfEdge({i, j});

                        // (i, j, k)
                        for (size_t k = 0; k < j; ++k) {
                            edgeCost += -subsetProblem.costOfTriple({i, j, k}) / 2.0;
                        }

                        // (i, k, j)
                        for (size_t k = j + 1; k < i; ++k) {
                            edgeCost += -subsetProblem.costOfTriple({i, k, j}) / 2.0;
                        }

                        // (k, i, j)
                        for (size_t k = i + 1; k < numberOfNodes; ++k) {
                            edgeCost += -subsetProblem.costOfTriple({k, i, j}) / 2.0;
                        }

                        edges[edgeIndex] = {i, j};
                        ws[edgeIndex] = edgeCost;

                        edgeIndex++;
                    }
                }

                undirected_graph g(edges, edges + numberOfEdges, ws, numberOfNodes, numberOfEdges);

                Rational cutWeight = boost::stoer_wagner_min_cut(g, get(boost::edge_weight, g));

                return -cutWeight;
            }

            template<class T>
            inline
            typename Persistency<T>::Rational
            Persistency<T>::subsetJoinCalculateRightHandSide(const std::set<size_t> &subset) {
                Rational lhs = 0;

                // vertex in U
                for (auto inVertex: subset) {
                    // vertex not in U
                    for (size_t outVertex = 0; outVertex < inVertex; ++outVertex) {
                        // edge fully contained in subset
                        if (subset.find(outVertex) != subset.end()) continue;

                        Rational edgeCost = problem_.costOfEdge({inVertex, outVertex});
                        if (edgeCost <= 0) {
                            lhs += edgeCost;
                        }

                        // (inVertex, outVertex, thirdVertex)
                        for (size_t thirdVertex = 0; thirdVertex < outVertex; ++thirdVertex) {
                            // all distinct
                            Rational tripleCost = problem_.costOfTriple({inVertex, outVertex, thirdVertex});
                            if (tripleCost <= 0) {
                                lhs += tripleCost / 2;
                            }
                        }

                        // (inVertex, thirdVertex, outVertex)
                        for (size_t thirdVertex = outVertex + 1; thirdVertex < inVertex; ++thirdVertex) {
                            Rational tripleCost = problem_.costOfTriple({inVertex, thirdVertex, outVertex});
                            if (tripleCost <= 0) {
                                lhs += tripleCost / 2;
                            }
                        }

                        // (thirdVertex, inVertex, outVertex)
                        for (size_t thirdVertex = inVertex + 1;
                             thirdVertex < problem_.numberOfElements(); ++thirdVertex) {
                            Rational tripleCost = problem_.costOfTriple({thirdVertex, inVertex, outVertex});
                            if (tripleCost <= 0) {
                                lhs += tripleCost / 2;
                            }
                        }
                    }

                    // (outVertex, inVertex)
                    for (size_t outVertex = inVertex; outVertex < problem_.numberOfElements(); ++outVertex) {
                        // edge fully contained in subset
                        if (subset.find(outVertex) != subset.end()) continue;

                        Rational edgeCost = problem_.costOfEdge({inVertex, outVertex});
                        if (edgeCost < 0) {
                            lhs += edgeCost;
                        }

                        // (outVertex, inVertex, thirdVertex)
                        for (size_t thirdVertex = 0; thirdVertex < inVertex; ++thirdVertex) {
                            // all distinct
                            Rational tripleCost = problem_.costOfTriple({outVertex, inVertex, thirdVertex});
                            if (tripleCost <= 0) {
                                lhs += tripleCost / 2;
                            }
                        }

                        // (outVertex, thirdVertex, inVertex)
                        for (size_t thirdVertex = inVertex + 1; thirdVertex < outVertex; ++thirdVertex) {
                            Rational tripleCost = problem_.costOfTriple({outVertex, thirdVertex, inVertex});
                            if (tripleCost <= 0) {
                                lhs += tripleCost / 2;
                            }
                        }

                        // (thirdVertex, outVertex, inVertex)
                        for (size_t thirdVertex = outVertex + 1;
                             thirdVertex < problem_.numberOfElements(); ++thirdVertex) {
                            Rational tripleCost = problem_.costOfTriple({thirdVertex, outVertex, inVertex});
                            if (tripleCost <= 0) {
                                lhs += tripleCost / 2;
                            }
                        }

                    }


                }

                return lhs;
            }

            template<class T>
            inline
            size_t
            Persistency<T>::remainingNodes() const {
                if (subProblemPersistencies_.empty()) {
                    return problem_.numberOfElements();
                } else {
                    size_t numberOfNodes = 0;
                    for (auto &subPersistency: subProblemPersistencies_) {
                        numberOfNodes += subPersistency.remainingNodes();
                    }
                    return numberOfNodes;
                }
            }

            template<class T>
            inline
            size_t
            Persistency<T>::remainingVariables() const {
                if (subProblemPersistencies_.empty()) {
                    size_t nodes = remainingNodes();

                    return nodes * (nodes - 1) / 2 - edges0_.size();
                } else {
                    size_t numberOfVariables = 0;
                    for (auto &subPersistency: subProblemPersistencies_) {
                        numberOfVariables += subPersistency.remainingVariables();
                    }
                    return numberOfVariables;
                }
            }

            template<class T>
            inline
            size_t
            Persistency<T>::remainingTriples() const {
                if (subProblemPersistencies_.empty()) {
                    size_t nodes = remainingNodes();
                    return nodes * (nodes - 1) * (nodes - 2) / 6 - triples0_.size();
                } else {
                    size_t numberOfTriplets = 0;
                    for (auto &subPersistency: subProblemPersistencies_) {
                        numberOfTriplets += subPersistency.remainingTriples();
                    }
                    return numberOfTriplets;
                }
            }


        }

    }
}


#endif
