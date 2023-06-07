#pragma once
#ifndef DSSC_PROBLEM_RESCALING_ADAPTOR_HXX
#define DSSC_PROBLEM_RESCALING_ADAPTOR_HXX

#include <cassert>

#include "dssc/math.hxx"

namespace dssc {

    enum struct RescalingCostType {
        PAIRS, TRIPLES, PAIRS_AND_TRIPLES
    };

/// Adaptor any problem obeying the contract look like a multicut cubic problem
/// This adapter scales the negative costs to [cmin, 0] and the positive costs to [0, cmax]
/// This is done for triplets or pairs or triplets and pairs
    template<class R, RescalingCostType CostType = RescalingCostType::TRIPLES>
    class CostRescalingAdaptor {
    public:
        typedef R Rational;
        typedef std::array<size_t, 3> IndexTriple;
        typedef std::array<size_t, 2> IndexEdge;

        // ProblemView needs to have the same public contract as below
        template<class PROBLEM_VIEW>
        CostRescalingAdaptor(
                PROBLEM_VIEW const &,
                Rational const &,
                Rational const &
        );

        template<class PROBLEM_VIEW>
        CostRescalingAdaptor(
                PROBLEM_VIEW const &,
                Rational const &,
                Rational const &,
                Rational const &,
                Rational const &
        );


        size_t numberOfElements() const;

        Rational costOfTriple(IndexTriple const &) const;

        Rational costOfEdge(IndexEdge const &) const;

    private:

        size_t indexOfTriple(IndexTriple const &) const;

        size_t indexOfTriple(size_t const, size_t const, size_t const) const;

        size_t indexOfEdge(IndexEdge const &) const;

        size_t indexOfEdge(size_t const, size_t const) const;

        template<class PROBLEM_VIEW>
        std::pair<Rational, Rational> calculateTripletBaseCosts(
                PROBLEM_VIEW const &
        ) const;

        template<class PROBLEM_VIEW>
        std::pair<Rational, Rational> calculateEdgeBaseCosts(
                PROBLEM_VIEW const &
        ) const;

        template<class PROBLEM_VIEW>
        void rescaleTripletCosts(
                PROBLEM_VIEW const &,
                Rational const &,
                Rational const &,
                Rational const &,
                Rational const &
        );

        template<class PROBLEM_VIEW>
        void rescalePairCosts(
                PROBLEM_VIEW const &,
                Rational const &,
                Rational const &,
                Rational const &,
                Rational const &
        );


        size_t numberOfElements_;
        std::vector<Rational> costsOfEdges_;
        std::vector<Rational> costsOfTriples_;
    };

    template<class R, RescalingCostType CostType>
    template<class PROBLEM_VIEW>
    void
    CostRescalingAdaptor<R, CostType>::rescalePairCosts(
            const PROBLEM_VIEW &problemView,
            const Rational &maxPositiveCost,
            const Rational &currentPositiveMaximum,
            const Rational &minNegativeCost,
            const Rational &currentNegativeMinimum
    ) {
        for (size_t j = 0; j < numberOfElements_; ++j) {
            for (size_t k = 0; k < j; ++k) {
                Rational edgeCost = problemView.costOfEdge({j, k});
                if (edgeCost > 0 && currentPositiveMaximum > 0) {
                    costsOfEdges_[indexOfEdge(j, k)] = edgeCost * maxPositiveCost / currentPositiveMaximum;
                } else if (edgeCost < 0 && currentNegativeMinimum < 0) {
                    costsOfEdges_[indexOfEdge(j, k)] = edgeCost * minNegativeCost / currentNegativeMinimum;
                }
            }
        }
    }

    template<class R, RescalingCostType CostType>
    template<class PROBLEM_VIEW>
    void
    CostRescalingAdaptor<R, CostType>::rescaleTripletCosts(
            const PROBLEM_VIEW &problemView,
            const Rational &maxPositiveCost,
            const Rational &currentPositiveMaximum,
            const Rational &minNegativeCost,
            const Rational &currentNegativeMinimum
    ) {
        for (size_t j = 0; j < numberOfElements_; ++j) {
            for (size_t k = 0; k < j; ++k) {
                for (size_t l = 0; l < k; ++l) {
                    Rational tripletCost = problemView.costOfTriple({j, k, l});
                    if (tripletCost > 0 && currentPositiveMaximum > 0) {
                        costsOfTriples_[indexOfTriple(j, k, l)] =
                                tripletCost * maxPositiveCost / currentPositiveMaximum;
                    } else if (tripletCost < 0 && currentNegativeMinimum < 0) {
                        costsOfTriples_[indexOfTriple(j, k, l)] =
                                tripletCost * minNegativeCost / currentNegativeMinimum;
                    }
                }
            }
        }

    }

    template<class R, RescalingCostType CostType>
    template<class PROBLEM_VIEW>
    std::pair<R, R>
    CostRescalingAdaptor<R, CostType>::calculateEdgeBaseCosts(const PROBLEM_VIEW &problemView) const {
        Rational edgePositiveMaximum = 0;
        Rational edgeNegativeMinimum = 0;

        for (size_t j = 0; j < numberOfElements_; ++j) {
            for (size_t k = 0; k < j; ++k) {
                Rational edgeCost = problemView.costOfEdge({j, k});

                if (edgeCost > 0) {
                    edgePositiveMaximum = std::max(edgeCost, edgePositiveMaximum);
                } else {
                    edgeNegativeMinimum = std::min(edgeCost, edgeNegativeMinimum);
                }
            }
        }

        return std::pair<Rational, Rational>(edgeNegativeMinimum, edgePositiveMaximum);
    }

    template<class R, RescalingCostType M>
    template<class PROBLEM_VIEW>
    std::pair<R, R>
    CostRescalingAdaptor<R, M>::calculateTripletBaseCosts(
            const PROBLEM_VIEW &problemView
    ) const {
        Rational tripletPositiveMaximum = 0;
        Rational tripletNegativeMinimum = 0;

        for (size_t j = 0; j < numberOfElements_; ++j) {
            for (size_t k = 0; k < j; ++k) {
                for (size_t l = 0; l < k; ++l) {
                    Rational tripletCost = problemView.costOfTriple({j, k, l});

                    if (tripletCost > 0) {
                        tripletPositiveMaximum = std::max(tripletCost, tripletPositiveMaximum);
                    } else {
                        tripletNegativeMinimum = std::min(tripletCost, tripletNegativeMinimum);
                    }
                }
            }
        }

        return std::pair<Rational, Rational>(tripletNegativeMinimum, tripletPositiveMaximum);
    }

    // constructor which calculates the min-max of the current problem and rescales such that these values are rescaled to minNegativeCost and maxPositiveCost
    template<class R, RescalingCostType M>
    template<class PROBLEM_VIEW>
    inline
    CostRescalingAdaptor<R, M>::CostRescalingAdaptor(
            PROBLEM_VIEW const &problemView,
            Rational const &minNegativeCost,
            Rational const &maxPositiveCost
    )
            :
            numberOfElements_(problemView.numberOfElements()),
            costsOfTriples_(numberOfElements_ * (numberOfElements_ - 1) * (numberOfElements_ - 2) / 6),
            costsOfEdges_(numberOfElements_ * (numberOfElements_ - 1) / 2) {

        if (M == RescalingCostType::TRIPLES || M == RescalingCostType::PAIRS_AND_TRIPLES) {
            std::pair<Rational, Rational> const &tripletRange = calculateTripletBaseCosts(problemView);
            rescaleTripletCosts(
                    problemView,
                    maxPositiveCost,
                    tripletRange.second,
                    minNegativeCost,
                    tripletRange.first
            );
        }

        if (M == RescalingCostType::PAIRS || M == RescalingCostType::PAIRS_AND_TRIPLES) {
            std::pair<Rational, Rational> const &edgeRange = calculateEdgeBaseCosts(problemView);
            rescalePairCosts(
                    problemView,
                    maxPositiveCost,
                    edgeRange.second,
                    minNegativeCost,
                    edgeRange.first
            );
        }
    }


    template<class R, RescalingCostType M>
    template<class PROBLEM_VIEW>
    inline
    CostRescalingAdaptor<R, M>::CostRescalingAdaptor(
            const PROBLEM_VIEW &problemView,
            const Rational &maxPositiveCost,
            const Rational &currentPositiveMaximum,
            const Rational &minNegativeCost,
            const Rational &currentNegativeMinimum
    ):
            numberOfElements_(problemView.numberOfElements()),
            costsOfTriples_(numberOfElements_ * (numberOfElements_ - 1) * (numberOfElements_ - 2) / 6),
            costsOfEdges_(numberOfElements_ * (numberOfElements_ - 1) / 2) {
        if (M == RescalingCostType::TRIPLES || M == RescalingCostType::PAIRS_AND_TRIPLES) {
            rescaleTripletCosts(
                    problemView,
                    maxPositiveCost,
                    currentPositiveMaximum,
                    minNegativeCost,
                    currentNegativeMinimum
            );
        }

        if (M == RescalingCostType::PAIRS || M == RescalingCostType::PAIRS_AND_TRIPLES) {
            rescalePairCosts(
                    problemView,
                    maxPositiveCost,
                    currentPositiveMaximum,
                    minNegativeCost,
                    currentNegativeMinimum
            );
        }
    }


    template<class R, RescalingCostType M>
    inline
    size_t
    CostRescalingAdaptor<R, M>::numberOfElements() const {
        return numberOfElements_;
    }

    template<class R, RescalingCostType M>
    inline
    R
    CostRescalingAdaptor<R, M>::costOfTriple(
            IndexTriple const &indexTriple
    ) const {

        return costsOfTriples_[indexOfTriple(indexTriple)];
    }

    template<class R, RescalingCostType M>
    inline
    R
    CostRescalingAdaptor<R, M>::costOfEdge(const IndexEdge &indexEdge) const {
        return costsOfEdges_[indexOfEdge(indexEdge)];
    }


    template<class R, RescalingCostType M>
    inline
    size_t
    CostRescalingAdaptor<R, M>::indexOfTriple(
            IndexTriple const &indexTriple
    ) const {
        assert(indexTriple[0] > indexTriple[1] && indexTriple[1] > indexTriple[2]);
        return indexOfTriple(indexTriple[0], indexTriple[1], indexTriple[2]);
    }

    template<class R, RescalingCostType M>
    inline
    size_t
    CostRescalingAdaptor<R, M>::indexOfEdge(
            IndexEdge const &indexEdge
    ) const {
        assert(indexEdge[0] > indexEdge[1]);
        return indexOfEdge(indexEdge[0], indexEdge[1]);

    }

/// Scalar index of strictly decreasing sequence of three indices.
///
/// This function always assumes j > k > l, regardless of the parameter
/// ASSUME_STRICTLY_DECREASING_TRIPLES.
///
    template<class R, RescalingCostType M>
    inline
    size_t
    CostRescalingAdaptor<R, M>::indexOfTriple(
            size_t const j,
            size_t const k,
            size_t const l
    ) const {
        assert(j > k && k > l);
        size_t const index = j * (j - 1) * (j - 2) / 6 + k * (k - 1) / 2 + l;
        return index;
    }

    /// Scalar index of strictly decreasing sequence of three indices.
    ///
    /// This function always assumes j > k > l, regardless of the parameter
    /// ASSUME_STRICTLY_DECREASING_TRIPLES.
    ///
    template<class R, RescalingCostType M>
    inline
    size_t
    CostRescalingAdaptor<R, M>::indexOfEdge(
            size_t const j,
            size_t const k
    ) const {
        assert(j > k);
        return j * (j - 1) / 2 + k;
    }


}

#endif
