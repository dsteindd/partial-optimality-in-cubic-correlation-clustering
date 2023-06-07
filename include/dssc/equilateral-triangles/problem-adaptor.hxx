#pragma once
#ifndef DSSC_EQUILATERAL_ADAPTOR_HXX
#define DSSC_EQUILATERAL_ADAPTOR_HXX

#include <cassert>

#include "math.h"
#include "dssc/math.hxx"
#include "problem.hxx"

namespace dssc {

    enum struct AngleCostMode {
        SUM_ANGLE_DIFFERENCES_ABS,
    };


    /// Adaptor that makes AffineProblem look like a cubic multicut problem.
    template<class R, AngleCostMode = AngleCostMode::SUM_ANGLE_DIFFERENCES_ABS>
    class ProblemEquilateralTrianglesAdaptor {
    public:
        typedef R Rational;
        typedef ProblemEquilateralTriangles<Rational> ProblemEquilateral;
        typedef typename ProblemEquilateral::Point2 Point2;
        typedef std::array<size_t, 3> IndexTriple;
        typedef std::array<size_t, 2> IndexEdge;

        // constructor for triplet and edge costs
        ProblemEquilateralTrianglesAdaptor(
                ProblemEquilateral const &,
                Rational const &,
                Rational const &,
                Rational const &,
                Rational const &,
                Rational const &
        );

        ProblemEquilateralTrianglesAdaptor(
                ProblemEquilateral const &,
                Rational const &,
                Rational const &,
                Rational const &,
                Rational const &,
                Rational const &,
                Rational const &
        );

        size_t numberOfElements() const;

        Rational costOfTriple(IndexTriple const &) const;

        Rational costOfEdge(IndexEdge const &) const;

    private:
        ProblemEquilateral const &problemEquilateral_;
        Rational const tripletCostAtAllPointsClose_;
        Rational const edgeCostAtDistanceZero_;
        Rational const edgeLengthThreshold_;
        Rational const tripletMaxEdgeLengthThreshold_;

        // constructor which we always used should stay intact
        Rational const angleDifferenceWeight_;
        Rational const costAtDifferenceZero_;

        Rational maxPairwiseDistance(Point2 const &p1, Point2 const &p2, Point2 const &p3) const;

        Rational minPairwiseDistance(Point2 const &p1, Point2 const &p2, Point2 const &p3) const;

        // calculates the metric for the angle difference
        Rational calculateAngleDifferenceToEquilateral(Point2 const &p1, Point2 const &p2, Point2 const &p3) const;
    };

    template<class R, AngleCostMode M>
    inline
    ProblemEquilateralTrianglesAdaptor<R, M>::ProblemEquilateralTrianglesAdaptor(
            ProblemEquilateral const &problemEquilateral,
            Rational const &angleDistanceOfCostZero,
            Rational const &tripletCostAtAngleDistanceZero,
            Rational const &edgeCostAtDistanceZero,
            Rational const &tripletMaxEdgeLengthThreshold,
            Rational const &edgeLengthThreshold
    )
            :   problemEquilateral_(problemEquilateral),
                tripletCostAtAllPointsClose_(tripletCostAtAngleDistanceZero),
                edgeCostAtDistanceZero_(edgeCostAtDistanceZero),
                tripletMaxEdgeLengthThreshold_(tripletMaxEdgeLengthThreshold),
                edgeLengthThreshold_(edgeLengthThreshold),
                angleDifferenceWeight_(-tripletCostAtAngleDistanceZero / angleDistanceOfCostZero),
                costAtDifferenceZero_(tripletCostAtAngleDistanceZero) {
        assert(tripletCostAtAngleDistanceZero < 0);
        assert(edgeCostAtDistanceZero < 0);
        assert(tripletMaxEdgeLengthThreshold >= 0);
        assert(edgeLengthThreshold >= 0);
    }

    template<class R, AngleCostMode M>
    inline
    ProblemEquilateralTrianglesAdaptor<R, M>::ProblemEquilateralTrianglesAdaptor(
            ProblemEquilateral const &problemEquilateral,
            Rational const &tripletAngleCostWeight,
            Rational const &costAtDistanceZero,
            Rational const &edgeCostAtDistanceZero,
            Rational const &tripletMaxEdgeLengthThreshold,
            Rational const &tripletCostAtAllPointsClose,
            Rational const &edgeLengthThreshold
    )
            :   problemEquilateral_(problemEquilateral),
                tripletCostAtAllPointsClose_(tripletCostAtAllPointsClose),
                edgeCostAtDistanceZero_(edgeCostAtDistanceZero),
                tripletMaxEdgeLengthThreshold_(tripletMaxEdgeLengthThreshold),
                edgeLengthThreshold_(edgeLengthThreshold),
                angleDifferenceWeight_(tripletAngleCostWeight),
                costAtDifferenceZero_(costAtDistanceZero) {
        assert(tripletCostAtAllPointsClose < 0);
        assert(edgeCostAtDistanceZero < 0);
        assert(tripletMaxEdgeLengthThreshold >= 0);
        assert(edgeLengthThreshold >= 0);
    }

    template<class R, AngleCostMode M>
    inline
    size_t
    ProblemEquilateralTrianglesAdaptor<R, M>::numberOfElements() const {
        return problemEquilateral_.points().size();
    }

    template<class R, AngleCostMode M>
    inline
    R
    ProblemEquilateralTrianglesAdaptor<R, M>::costOfTriple(
            IndexTriple const &indexTriple
    ) const {
        Point2 const &p0 = problemEquilateral_.points()[indexTriple[0]];
        Point2 const &p1 = problemEquilateral_.points()[indexTriple[1]];
        Point2 const &p2 = problemEquilateral_.points()[indexTriple[2]];

        Rational const maxDistance = maxPairwiseDistance(p0, p1, p2);

        if (maxDistance < tripletMaxEdgeLengthThreshold_) {
            // all three points are close
            return tripletCostAtAllPointsClose_ * (1 - maxDistance / tripletMaxEdgeLengthThreshold_);
        }
        Rational const minDistance = minPairwiseDistance(p0, p1, p2);

        if (minDistance < edgeLengthThreshold_) {
            return 0;
        }
        Rational const angleDifference = calculateAngleDifferenceToEquilateral(p0, p1, p2);

        Rational const cost = angleDifferenceWeight_ * angleDifference + costAtDifferenceZero_;

        return cost;
    }

    template<class R, AngleCostMode M>
    inline
    R
    ProblemEquilateralTrianglesAdaptor<R, M>::calculateAngleDifferenceToEquilateral(Point2 const &p1,
                                                                                    Point2 const &p2,
                                                                                    Point2 const &p3) const {

        Rational const dphi0 = angleBetween(p2 - p1, p3 - p1) - M_PI / 3;
        Rational const dphi1 = angleBetween(p3 - p2, p1 - p2) - M_PI / 3;
        Rational const dphi2 = angleBetween(p2 - p3, p1 - p3) - M_PI / 3;

        if (M == AngleCostMode::SUM_ANGLE_DIFFERENCES_ABS) {
            return abs(dphi0) + abs(dphi1) + abs(dphi2);
        } else {
            throw std::runtime_error("unknown cost mode to calculate an angle difference.");
        }

    }


    // now returns negative cost if distance is sufficiently small and positive if not
    template<class R, AngleCostMode M>
    inline
    R
    ProblemEquilateralTrianglesAdaptor<R, M>::costOfEdge(const IndexEdge &indexEdge) const {
        Point2 const &p0 = problemEquilateral_.points()[indexEdge[0]];
        Point2 const &p1 = problemEquilateral_.points()[indexEdge[1]];

        double const d = distance(p0, p1);

        // if distance is small return something negative
        if (d < edgeLengthThreshold_) {
            // all d <= 2sigma -> want to have -3 for d = 0 and d = 0 for d = 2sigma
            // c = -3 + d* (--3)/ 2sigma = -3 - (-3)*d/(2sigma) = c0 - c0 * d/(d0) = c0 (1- d/d0)
            return edgeCostAtDistanceZero_ * (1 - d / edgeLengthThreshold_);
        } else {
            return 0;
        }
    }

    template<class R, AngleCostMode M>
    inline
    R
    ProblemEquilateralTrianglesAdaptor<R, M>::maxPairwiseDistance(
            Point2 const &p1,
            Point2 const &p2,
            Point2 const &p3
    ) const {
        Rational const d12 = distance(p1, p2);
        Rational const d13 = distance(p1, p3);
        Rational const d23 = distance(p2, p3);

        return std::max({d12, d13, d23});
    }

    template<class R, AngleCostMode M>
    inline
    R
    ProblemEquilateralTrianglesAdaptor<R, M>::minPairwiseDistance(
            const Point2 &p1,
            const Point2 &p2,
            const Point2 &p3
    ) const {
        Rational const d12 = distance(p1, p2);
        Rational const d13 = distance(p1, p3);
        Rational const d23 = distance(p2, p3);

        return std::min({d12, d13, d23});
    }

}

#endif
