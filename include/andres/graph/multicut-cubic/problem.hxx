#pragma once
#ifndef ANDRES_GRAPH_MULTICUT_CUBIC_PROBLEM_HXX
#define ANDRES_GRAPH_MULTICUT_CUBIC_PROBLEM_HXX

#include <cassert>
#include <array>
#include <vector>
#include <functional>
#include <algorithm>
#include <map>
#include <set>

#include "andres/partition.hxx" // todo: move this out to a different class

namespace andres {
namespace graph {
namespace multicut_cubic {

template<class T = double, bool ASSUME_STRICTLY_DECREASING_TRIPLES = false>
class Problem {
public:
    typedef T Rational;
    typedef std::array<size_t, 2> IndexPair;
    typedef std::array<size_t, 3> IndexTriple;

    Problem(size_t const = 0);
    size_t numberOfElements() const;
    Rational & costOfEdge(IndexPair const &);
    Rational costOfEdge(IndexPair const &) const;
    Rational & costOfTriple(IndexTriple const &);
    Rational costOfTriple(IndexTriple const &) const;
    Rational & costConstant();
    Rational costConstant() const;

    Problem onlyNegativePartCosts();
    Problem absoluteCosts();
    Problem unorderedProjectTo(std::set<size_t> const &);
    Problem orderedProjectTo(std::vector<size_t> const &);

    // todo: move this out to a different class
    typedef Partition<size_t> OnePartition;
    void contract(size_t const, size_t const);
    void contract(size_t const, size_t const, size_t const);
    // neither const nor ref, since we need to sort eventually and the client should not be confused by an altered array
    void contract(std::vector<size_t>);
    void swap(size_t const, size_t const);
    OnePartition partition() const;

private:    
    size_t indexOfTriple(IndexTriple const &) const;
    size_t indexOfTriple(size_t const, size_t const, size_t const) const;
    size_t indexOfEdge(IndexPair const &) const;
    size_t indexOfEdge(size_t const, size_t const) const;

    size_t numberOfElements_;
    Rational constant_;
    std::vector<Rational> costsOfEdges_;
    std::vector<Rational> costsOfTriples_;

    // todo: move this out to a different class
    // maps indices to initial node indices
    std::map<size_t, size_t> labels_;
    OnePartition partition_;
};

template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
inline
Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::Problem(
    size_t const numberOfElements
)
:   numberOfElements_(numberOfElements),
    constant_(),
    costsOfEdges_(numberOfElements * (numberOfElements - 1) / 2),
    costsOfTriples_(numberOfElements * (numberOfElements - 1) * (numberOfElements - 2) / 6),
    labels_(),
    partition_(numberOfElements)
{
        for (size_t i = 0; i < numberOfElements; i++){
            labels_[i] = i;
        }
}

template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
inline
typename Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::Rational
Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::costConstant() const {
    return constant_;
}

template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
inline
typename Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::Rational &
Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::costConstant() {
    return constant_;
}

template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
inline
typename Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::Rational &
Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::costOfTriple(
    IndexTriple const & triple
) {
    return costsOfTriples_[indexOfTriple(triple)];
}

template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
inline
typename Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::Rational
Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::costOfTriple(
        IndexTriple const & indexTriple
) const {
    return costsOfTriples_[indexOfTriple(indexTriple)];
}

template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
inline
typename Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::Rational &
Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::costOfEdge(
        IndexPair const & indexPair
        ) {
    return costsOfEdges_[indexOfEdge(indexPair)];
}

template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
inline
typename Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::Rational
Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::costOfEdge(
        IndexPair const & indexPair
) const {
    return costsOfEdges_[indexOfEdge(indexPair)];
}

template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
inline
size_t
Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::numberOfElements() const {
    return numberOfElements_;
}

template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
inline
void
Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::contract(
    size_t const s,
    size_t const t
) {
    assert(s < numberOfElements());
    assert(t < numberOfElements());
    assert (s > t);

    if(s == t) {
        return;
    }

    // adjust constant
    if (ASSUME_STRICTLY_DECREASING_TRIPLES){
        IndexPair stEdge = {s, t};
        std::sort(stEdge.begin(), stEdge.end(), std::greater<size_t>());

        constant_ += costOfEdge(stEdge);
    } else {
        constant_ += costOfEdge({s, t});
    }

    // adjust edge weights
    for (size_t j = 0; j < numberOfElements(); ++j){
        if (j == s || j == t) continue;
        if (ASSUME_STRICTLY_DECREASING_TRIPLES){
            IndexPair sEdge = {j, s};
            IndexPair tEdge = {j, t};
            IndexTriple stTriple = {j, s, t};
            std::sort(sEdge.begin(), sEdge.end(), std::greater<size_t>());
            std::sort(tEdge.begin(), tEdge.end(), std::greater<size_t>());
            std::sort(stTriple.begin(), stTriple.end(), std::greater<size_t>());

            costOfEdge(tEdge) += costOfEdge(sEdge) + costOfTriple(stTriple);
        } else {
            costOfEdge({j, t}) += costOfEdge({j, s}) + costOfTriple({j, s, t});
        }
    }

    // add cost of triples (j, k, s) to costs of triples (j, k, t)
    for(size_t j = 0; j < numberOfElements(); ++j) {
        if(j == s || j == t) continue;
        for(size_t k = 0; k < j; ++k) {
            if(k == s || k == t) continue;
            if(ASSUME_STRICTLY_DECREASING_TRIPLES) {
                IndexTriple tTriple = {j, k, t};
                IndexTriple sTriple = {j, k, s};
                IndexTriple jTriple = {j, s, t};
                IndexTriple kTriple = {k, s, t};
                std::sort(tTriple.begin(), tTriple.end(), std::greater<size_t>());
                std::sort(sTriple.begin(), sTriple.end(), std::greater<size_t>());
                std::sort(jTriple.begin(), jTriple.end(), std::greater<size_t>());
                std::sort(kTriple.begin(), kTriple.end(), std::greater<size_t>());
//                costOfTriple(tTriple) += costOfTriple(sTriple) + costOfTriple(jTriple) + costOfTriple(kTriple);

                costOfTriple(tTriple) += costOfTriple(sTriple);
            }
            else {
//                costOfTriple({j, k, t}) += costOfTriple({j, k, s}) + costOfTriple({j, s, t}) + costOfTriple({k, s, t});

                costOfTriple({j, k, t}) += costOfTriple({j, k, s});
            }
        }
    }

    size_t const last = numberOfElements() - 1;
    if(s != last) {
        // swap s with last element
        for(size_t j = 0; j < numberOfElements(); ++j) {
            if(j == s || j == last) continue;
            for(size_t k = 0; k < j; ++k) {
                if(k == s || k == last) continue;
                if(ASSUME_STRICTLY_DECREASING_TRIPLES) {
                    IndexTriple sTriple = {j, k, s};
                    IndexTriple lastTriple = {j, k, last};
                    std::sort(sTriple.begin(), sTriple.end(), std::greater<size_t>());
                    std::sort(lastTriple.begin(), lastTriple.end(), std::greater<size_t>());
                    costOfTriple(sTriple) = costOfTriple(lastTriple);
                }
                else {
                    costOfTriple({j, k, s}) = costOfTriple({j, k, last});
                }
            }
        }

        for (size_t j = 0; j < numberOfElements(); ++j){
            if (j == s || j == last) continue;
            if (ASSUME_STRICTLY_DECREASING_TRIPLES){
                IndexPair sEdge = {j, s};
                IndexPair lastEdge = {j, last};
                std::sort(sEdge.begin(), sEdge.end(), std::greater<size_t>());
                std::sort(lastEdge.begin(), lastEdge.end(), std::greater<size_t>());
                costOfEdge(sEdge) = costOfEdge(lastEdge);
            } else {
                costOfEdge({j, s}) = costOfEdge({j, last});
            }
        }


    }
    --numberOfElements_;
    costsOfTriples_.resize(numberOfElements() * (numberOfElements() - 1) * (numberOfElements() - 2) / 6);
    costsOfEdges_.resize(numberOfElements() * (numberOfElements() - 1) / 2);

    // labels_ is a mapping of index (in next iteration) -> actual initial element
    // last is labeled s in subsequent iteration

    partition_.merge(labels_[s], labels_[t]);
    if (s != last){
        labels_[s] = labels_[last];
    }

}


template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
inline
size_t
Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::indexOfTriple(
    IndexTriple const & indexTriple
) const {
    if(ASSUME_STRICTLY_DECREASING_TRIPLES) {
        assert(indexTriple[0] > indexTriple[1] && indexTriple[1] > indexTriple[2]);
        return indexOfTriple(indexTriple[0], indexTriple[1], indexTriple[2]);
    }
    else {
        IndexTriple t = indexTriple; // copy
        std::sort(t.begin(), t.end(), std::greater<size_t>());
        assert(t[0] > t[1] && t[1] > t[2]);
        return indexOfTriple(t[0], t[1], t[2]);
    }
}

/// Scalar index of strictly decreasing sequence of three indices.
///
/// This function always assumes j > k > l, regardless of the parameter
/// ASSUME_STRICTLY_DECREASING_TRIPLES.
///
template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
inline
size_t
Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::indexOfTriple(
    size_t const j,
    size_t const k,
    size_t const l
) const {
    assert(j > k && k > l);
    return j * (j - 1) * (j - 2) / 6 + k * (k - 1) / 2 + l;
}

template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
inline
size_t
Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::indexOfEdge(
        IndexPair const & indexPair
) const {
    if(ASSUME_STRICTLY_DECREASING_TRIPLES) {
        assert(indexPair[0] > indexPair[1]);
        return indexOfEdge(indexPair[0], indexPair[1]);
    }
    else {
        IndexPair e = indexPair; // copy
        std::sort(e.begin(), e.end(), std::greater<>());
        assert(e[0] > e[1]);
        return indexOfEdge(e[0], e[1]);
    }
}

/// Scalar index of strictly decreasing sequence of three indices.
///
/// This function always assumes j > k > l, regardless of the parameter
/// ASSUME_STRICTLY_DECREASING_TRIPLES.
///
template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
inline
size_t
Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::indexOfEdge(
        size_t const j,
        size_t const k
) const {
    assert(j > k);
    return j * (j - 1) / 2 + k;
}


template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
typename Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::OnePartition
Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::partition() const {
    return partition_;
}

template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
inline
void
Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::contract(const size_t u, const size_t v, const size_t w) {
    assert(u != v && v != w && u != w);


    size_t oldNumberOfElements = numberOfElements();

    IndexTriple triple = {u, v, w};
    std::sort(triple.begin(), triple.end(), std::greater<>());

//    contract(triple);

    // triple is sorted decreasing
    contract(triple[1], triple[2]);

    // index triple[1] corresponds to node last / triple[2] corresponds to merge node
    // if triple[0] was last, then we need to merge triple[2], triple[1] again
    if (triple[0] == oldNumberOfElements - 1){
        // triple[0] was last node, and has been renamed to triple[2]
        contract(triple[1], triple[2]);
    } else {
        // triple[0] has the correct index
        contract(triple[0], triple[2]);
    }

    // merging (u, v) means that the min(u, v) corresponds to merge node {u, v}
    // max(u, v) corresponds to last or is deleted if max(u, v) = numElements - 1
//    // case 1: u and v < numElements - 1
//    if (u < v){
//        // u corresponds to merge node
//        if (w != oldNumberOfElements - 1){
//            // w was not last node and not renamed
//            contract(u, w);
//        } else {
//            // w was last node and renamed to v
//            contract(u, v);
//        }
//    } else {
//        // v corresponds to merge node
//        if (w != oldNumberOfElements - 1){
//            // w was not last node and not renamed
//            contract(v, w);
//        } else {
//            // w was last node and renamed to u
//            contract(u, v);
//        }
//    }

}

template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
inline
void
Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::contract(
        std::vector<size_t> indicesToContract) {
    // assert all unequal
    const size_t m = indicesToContract.size();

    for (size_t index = 0; index < m - 1; ++index){
        assert(indicesToContract[index] != indicesToContract[index + 1]);
    }

    std::sort(indicesToContract.begin(), indicesToContract.end(), std::greater<>());

    for (size_t i = 0; i < m; ++i){
        size_t index = indicesToContract[i];
        size_t last = numberOfElements_ - 1 - i;

        swap(index, last);
    }


    // at this point, the indices to be merged should be
    // [n-1], [n-2], [n-indicesToContract.size()]

    // we define [n-indicesToContract.size()] to be the new merge node
    const size_t mergeNodeIndex = numberOfElements_ - m;

    // update constant
    for (size_t i = mergeNodeIndex; i < numberOfElements_; ++i){
        for (size_t j = mergeNodeIndex; j < i; ++j){
            constant_ += costOfEdge({i, j});

            for (size_t k = mergeNodeIndex; k < j; ++k){
                constant_ += costOfTriple({i, j, k});
            }
        }
    }

    // update edges (mergeNodeIndex, i)
    // i not in merge set
    for (size_t i = 0; i < mergeNodeIndex; ++i){

        Rational newCost = 0;
        // j in merge set
        // add all edgeCost(j, i)
        for (size_t j = mergeNodeIndex; j < numberOfElements_; ++j){
            newCost += costOfEdge({j, i});

            // all (j, k, i)
            for (size_t k = mergeNodeIndex; k < j; ++k){
                newCost += costOfTriple({j, k, i});
            }
        }

        costOfEdge({mergeNodeIndex, i}) = newCost;
    }

    // update triplet costs
    // update triples (mergeNodeIndex, i, j)
    for (size_t i = 1; i < mergeNodeIndex; ++i){
        for (size_t j = 0; j < i; ++j){
            Rational newCost = 0;

            for (size_t k = mergeNodeIndex; k < numberOfElements_; ++k){
                // add (k, i, j)
                newCost += costOfTriple({k, i, j});
            }

            costOfTriple({mergeNodeIndex, i, j}) = newCost;
        }
    }

    // contract in partition
    for (size_t i = mergeNodeIndex; i < numberOfElements_ - 1; ++i){
        partition_.merge(labels_[i], labels_[i+1]);
    }

    // then deduct indicesToContract.size() - 1 from numberOfElements_
    numberOfElements_ = mergeNodeIndex + 1;
    costsOfEdges_.resize(numberOfElements_ * (numberOfElements_ - 1) / 2);
    costsOfTriples_.resize(numberOfElements_ * (numberOfElements_ - 1) * (numberOfElements_ - 2) / 6);

}

template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
inline
void
Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::swap(const size_t s, const size_t t) {
    // already correct indices, return
    if (s == t) return;

    // edge costs are not affected

    // triplet costs
    // c(s, t, i) is not affected
    // swap all c(s, i, j) with c(t, i, j)
    for (size_t i = 0; i < numberOfElements_; ++i){
        if (i == s || i == t) continue;

        // swap edgeCost(s, i) with edgeCost(t, i)
        IndexPair sEdge = {s, i};
        IndexPair tEdge = {t, i};

        if (ASSUME_STRICTLY_DECREASING_TRIPLES){
            std::sort(sEdge.begin(), sEdge.end(), std::greater<>());
            std::sort(tEdge.begin(), tEdge.end(), std::greater<>());

            Rational sEdgeCost = costOfEdge(sEdge);
            costOfEdge(sEdge) = costOfEdge(tEdge);
            costOfEdge(tEdge) = sEdgeCost;
        } else {
            Rational sEdgeCost = costOfEdge(sEdge);
            costOfEdge(sEdge) = costOfEdge(tEdge);
            costOfEdge(tEdge) = sEdgeCost;
        }


        for (size_t j = 0; j < i; ++j){
            if (j == s || j == t) continue;

            IndexTriple sTriple = {i, j, s};
            IndexTriple tTriple = {i, j, t};

            if (ASSUME_STRICTLY_DECREASING_TRIPLES){
                std::sort(sTriple.begin(), sTriple.end(), std::greater<>());
                std::sort(tTriple.begin(), tTriple.end(), std::greater<>());

                Rational sTripleCost = costOfTriple(sTriple);
                costOfTriple(sTriple) = costOfTriple(tTriple);
                costOfTriple(tTriple) = sTripleCost;
            } else {
                Rational sTripleCost = costOfTriple(sTriple);
                costOfTriple(sTriple) = costOfTriple(tTriple);
                costOfTriple(tTriple) = sTripleCost;
            }
        }
    }

    // swap index labels
    size_t sLabel = labels_[s];
    labels_[s] = labels_[t];
    labels_[t] = sLabel;
}


    // This function returns a new problem instance with all positive edge and triplet costs set to zero
    template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
    Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES> Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::onlyNegativePartCosts() {
        Problem problem(numberOfElements());
        for (size_t i = 0; i < numberOfElements(); ++i){
            for (size_t j = 0; j < i; ++j){
                for (size_t k = 0; k < j; ++k){
                    if (costOfTriple({i, j, k}) > 0){
                        problem.costOfTriple({i, j, k}) = 0;
                    } else {
                        problem.costOfTriple({i, j, k}) = -costOfTriple({i, j, k});
                    }
                }
            }
        }

        for (size_t i = 0; i < numberOfElements(); ++i){
            for (size_t j = 0; j < i; ++j){
                if (costOfEdge({i, j}) > 0){
                    problem.costOfEdge({i, j}) = 0;
                } else {
                    problem.costOfEdge({i, j}) = -costOfEdge({i, j});
                }
            }
        }

        return problem;
    }

    template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
    Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES> Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::absoluteCosts() {
        Problem problem(numberOfElements());
        for (size_t i = 0; i < numberOfElements(); ++i){
            for (size_t j = 0; j < i; ++j){
                for (size_t k = 0; k < j; ++k){
                    problem.costOfTriple({i, j, k}) = costOfTriple({i, j, k}) < 0 ? - costOfTriple({i, j, k}) : costOfTriple({i, j, k});
                }
            }
        }

        for (size_t i = 0; i < numberOfElements(); ++i){
            for (size_t j = 0; j < i; ++j){
                problem.costOfEdge({i, j}) = costOfEdge({i, j}) < 0 ? - costOfEdge({i, j}) : costOfEdge({i, j});
                }
            }

        return problem;
    }

    template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
    Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES> Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::unorderedProjectTo(
            const std::set<size_t> & subset) {

        std::vector<size_t> subsetArray(subset.size());
        std::copy(subset.begin(), subset.end(), subsetArray.begin());
        std::sort(subsetArray.begin(), subsetArray.end(), std::less<>());

        return orderedProjectTo(subsetArray);
    }

    template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
    Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES> Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::orderedProjectTo(
            const std::vector<size_t> & subset
            ) {
        Problem problem(subset.size());

        for (size_t i = 0; i < subset.size(); ++i){
            for (size_t j = 0; j < i; ++j){
                if (ASSUME_STRICTLY_DECREASING_TRIPLES){
                    IndexPair fromProjectedEdge = {subset[i], subset[j]};
                    std::sort(fromProjectedEdge.begin(), fromProjectedEdge.end(), std::greater<>());
                    problem.costOfEdge({i, j}) = costOfEdge(fromProjectedEdge);
                } else {
                    problem.costOfEdge({i, j}) = costOfEdge({subset[i], subset[j]});
                }


                for (size_t k = 0; k < j; ++k){
                    if (ASSUME_STRICTLY_DECREASING_TRIPLES){
                        IndexTriple fromProjectedTriple = {subset[i], subset[j], subset[k]};
                        std::sort(fromProjectedTriple.begin(), fromProjectedTriple.end(), std::greater<>());
                        problem.costOfTriple({i, j, k}) = costOfTriple(fromProjectedTriple);
                    } else {
                        problem.costOfTriple({i, j, k}) = costOfTriple({subset[i], subset[j], subset[k]});
                    }
                }
            }
        }
        return problem;
    }
}
}
}

#endif
