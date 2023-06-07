#pragma once
#ifndef DSSC_COST_MULTIPLICATION_ADAPTOR_HXX
#define DSSC_COST_MULTIPLICATION_ADAPTOR_HXX

#include <cassert>

#include "dssc/math.hxx"

namespace dssc {

/// Adaptor any problem obeying the contract look like a multicut cubic problem
/// This adapter scales the negative costs to [cmin, 0] and the positive costs to [0, cmax]
/// This is done independently for edges and triplets
    template<class PROBLEM_VIEW, class R>
    class CostMultiplicationAdaptor {
    public:
        typedef R Rational;
        typedef PROBLEM_VIEW ProblemView;
        typedef std::array<size_t, 3> IndexTriple;
        typedef std::array<size_t, 2> IndexEdge;

        // ProblemView needs to have the same public contract as below
        explicit CostMultiplicationAdaptor(
                ProblemView const &,
                Rational const & = 1.0,
                Rational const & = 1.0
        );
        size_t numberOfElements() const;
        Rational costOfTriple(IndexTriple const &) const;
        Rational costOfEdge(IndexEdge const &) const;

    private:
        Rational const & edgeScale_;
        Rational const & tripletScale_;
        ProblemView  const & problemView_;
    };

    template<class PROBLEM_VIEW, class R>
    inline
    CostMultiplicationAdaptor<PROBLEM_VIEW, R>::CostMultiplicationAdaptor(
            ProblemView const & problemView,
            Rational const & edgeScale,
            Rational const & tripletScale
    )
            :
            problemView_(problemView),
            edgeScale_(edgeScale),
            tripletScale_(tripletScale)
    {
    }

    template<class PROBLEM_VIEW, class T>
    inline
    size_t
    CostMultiplicationAdaptor<PROBLEM_VIEW, T>::numberOfElements() const {
        return problemView_.numberOfElements();
    }

    template<class PROBLEM_VIEW, class T>
    inline
    typename CostMultiplicationAdaptor<PROBLEM_VIEW, T>::Rational
    CostMultiplicationAdaptor<PROBLEM_VIEW, T>::costOfEdge(IndexEdge const & indexEdge) const {
        return edgeScale_*problemView_.costOfEdge(indexEdge);
    }

    template<class PROBLEM_VIEW, class T>
    inline
    typename CostMultiplicationAdaptor<PROBLEM_VIEW, T>::Rational
    CostMultiplicationAdaptor<PROBLEM_VIEW, T>::costOfTriple(IndexTriple const & indexTriple) const {
        return tripletScale_*problemView_.costOfTriple(indexTriple);
    }
}

#endif
