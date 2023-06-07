#pragma once

#include <cmath>
#include <random>
#include <array>
#include <vector>
#include <fstream>

#include "dssc/math.hxx"

namespace dssc {


/// Problem of clustering 1-dimensional affine subspaces in 2D space.
    template<class R = double>
    class ProblemEquilateralTriangles {
    public:
        typedef R Rational;
        typedef Point<Rational, 2> Point2;
        typedef std::vector<Point2> PointContainer;

        struct Triangle {
            Triangle(){};
            Triangle(const Point2 midPoint, const Rational radius = 0.5, const Rational angleOffset = 0.2) {
                point1_ = toPoint(midPoint, radius, angleOffset);
                point2_ = toPoint(midPoint, radius, angleOffset + 2*M_PI / 3);
                point3_ = toPoint(midPoint, radius, angleOffset + 2*2*M_PI / 3);
            }

            Triangle(const Point2 p1, const Point2 p2, const Point2 p3){
                point1_ = p1;
                point2_ = p2;
                point3_ = p3;
            }

            std::vector<Point2> const points() const{
                return {point1_, point2_, point3_};
            };

            Point2 point1_;
            Point2 point2_;
            Point2 point3_;
        };

        typedef std::vector<Triangle> TriangleContainer;
        typedef std::array<size_t, 3> IndexTriple;
        typedef size_t Label;
        typedef std::vector<Label> LabelContainer;

        // generate problem, given mid points of triangles, radii and angle offsets
        ProblemEquilateralTriangles(
                std::vector<Point2> const &,
                std::vector<Rational> const &,
                std::vector<Rational> const &,
                std::vector<std::vector<size_t>> const &,
                Rational const,
                size_t const & =4242);

        ProblemEquilateralTriangles(size_t const, size_t const, Rational const, Rational const = 0.2, Rational const = 0.2, size_t const & = 4242);

        ProblemEquilateralTriangles(
                    std::vector<Point2> const &
                );

        PointContainer const &points() const;

        TriangleContainer const &triangles() const;

        LabelContainer const &labels() const;

    private:
        Rational sigma_;
        TriangleContainer triangles_;
        PointContainer points_;
        LabelContainer labels_;

        Point2 const drawRandomDisplacementVector(std::default_random_engine &);
        static Point2 toPoint(Point2 const, Rational const, Rational const );
        bool isTriangleInUnitSquare(Triangle const &);
        bool isPointInUnitSquare(Point2 const &);

    };

/// Define problem with fixed midpoints, radii and angleOffsets
    template<class C>
    inline
    ProblemEquilateralTriangles<C>::ProblemEquilateralTriangles(
            std::vector<Point2> const &midpoints,
            std::vector<Rational> const &radii,
            std::vector<Rational> const &angleOffsets,
            std::vector<std::vector<size_t>> const & pointsAtTriangleNodes,
            Rational const sigma,
            size_t const& seed
    )
            :   sigma_(sigma),
                triangles_(midpoints.size()) {
        // assert size correctness
        assert(midpoints.size() == radii.size() && radii.size() == angleOffsets.size());

        for (const Rational radius: radii) assert(radius > 0);

        assert(pointsAtTriangleNodes.size() == midpoints.size());
        // calculate number of points
        size_t numOfPoints = 0;
        for (auto const & trianglePoints : pointsAtTriangleNodes){
            assert(trianglePoints.size() == 3);
            for (auto const & pointsAtNode : trianglePoints){
                numOfPoints += pointsAtNode;
            }
        }

        points_ = std::vector<Point2>(numOfPoints);
        labels_ = std::vector<size_t >(numOfPoints);

        std::default_random_engine randomEngine(seed);

        size_t totalDrawnPoints = 0;

        for (size_t tIndex = 0; tIndex < midpoints.size(); ++tIndex) {
            Triangle triangle(midpoints[tIndex], radii[tIndex], angleOffsets[tIndex]);
//            assert(isTriangleInUnitSquare(triangle));
            triangles_[tIndex] = triangle;

            for (size_t nodeIndex = 0; nodeIndex < 3; ++nodeIndex){
                size_t const pointsToBeDrawn = pointsAtTriangleNodes[tIndex][nodeIndex];


                for (size_t pIndex = totalDrawnPoints; pIndex < totalDrawnPoints + pointsToBeDrawn; ++pIndex){

                    Point2 p = triangle.points()[nodeIndex] + drawRandomDisplacementVector(randomEngine);

//                    do {
//                        p = triangle.points()[nodeIndex] + drawRandomDisplacementVector(randomEngine);
//                    } while (!isPointInUnitSquare(p));

                    points_[pIndex] = p;
                    labels_[pIndex] = tIndex;
                }

                totalDrawnPoints += pointsAtTriangleNodes[tIndex][nodeIndex];
            }
        }

    }

/// Define problem with random triangles
    template<class C>
    inline
    ProblemEquilateralTriangles<C>::ProblemEquilateralTriangles(
            size_t const numberOfTriangles,
            size_t const numberOfNoisePoints,
            Rational const sigma,
            Rational const minRadius,
            Rational const maxRadius,
            size_t const & seed
    )
            :   sigma_(sigma),
                triangles_(numberOfTriangles),
                points_(3 * numberOfTriangles + numberOfNoisePoints),
                labels_(3 * numberOfTriangles + numberOfNoisePoints)
                {
        // distribute midpoints in unit square
        std::default_random_engine randomEngine(seed);
        std::uniform_real_distribution<Rational> coordinateDistribution(0, 1);
        std::uniform_real_distribution<Rational> angleDistribution(0, 2 * M_PI);
        std::uniform_real_distribution<Rational> radiusDistribution(minRadius, maxRadius);

        for (size_t index = 0; index < numberOfTriangles; ++index) {
            Triangle distortedTriangle;
            Triangle trueTriangle;

            do {

                const Rational x = coordinateDistribution(randomEngine);
                const Rational y = coordinateDistribution(randomEngine);
                const Rational r = radiusDistribution(randomEngine);
                const Rational phi = angleDistribution(randomEngine);
                trueTriangle = Triangle(Point2({x, y}), r, phi);

                const Point2 p1 = trueTriangle.point1_ + drawRandomDisplacementVector(randomEngine);
                const Point2 p2 = trueTriangle.point2_ + drawRandomDisplacementVector(randomEngine);
                const Point2 p3 = trueTriangle.point3_ + drawRandomDisplacementVector(randomEngine);

                distortedTriangle = Triangle(p1, p2, p3);
            } while (!isTriangleInUnitSquare(distortedTriangle) || !isTriangleInUnitSquare(trueTriangle));


            triangles_[index] = trueTriangle;

            points_[3 * index] = distortedTriangle.point1_;
            labels_[3 * index] = index;

            points_[3 * index + 1] = distortedTriangle.point2_;
            labels_[3 * index + 1] = index;

            points_[3 * index + 2] = distortedTriangle.point3_;
            labels_[3 * index + 2] = index;

        }

        for (size_t index = 3*numberOfTriangles; index < 3*numberOfTriangles + numberOfNoisePoints; ++index){
            const Rational x = coordinateDistribution(randomEngine);
            const Rational y = coordinateDistribution(randomEngine);

            Point2 p({x, y});
            points_[index] = p;
            labels_[index] = index;
        }
    }

    template<class C>
    inline
    ProblemEquilateralTriangles<C>::ProblemEquilateralTriangles(
            const std::vector<Point2> & points
            ):
            points_(points)
            {}

    template<class C>
    inline
    typename ProblemEquilateralTriangles<C>::PointContainer const &
    ProblemEquilateralTriangles<C>::points() const {
        return points_;
    }

    template<class C>
    inline
    typename ProblemEquilateralTriangles<C>::LabelContainer const &
    ProblemEquilateralTriangles<C>::labels() const {
        return labels_;
    }

    template<class C>
    inline
    typename ProblemEquilateralTriangles<C>::Point2 const
    ProblemEquilateralTriangles<C>::drawRandomDisplacementVector(std::default_random_engine & randomEngine) {

        std::normal_distribution<Rational > displacementDistribution(0.0, sigma_);

        const Rational dx = displacementDistribution(randomEngine);
        const Rational dy = displacementDistribution(randomEngine);

        return Point2({dx, dy});
    }

    template<class C>
    inline
    typename ProblemEquilateralTriangles<C>::TriangleContainer const &
    ProblemEquilateralTriangles<C>::triangles() const {
        return triangles_;
    }

    template<class C>
    inline
    typename ProblemEquilateralTriangles<C>::Point2
    ProblemEquilateralTriangles<C>::toPoint(const Point2 midPoint, const Rational radius, const Rational angleOffset) {
        const Point2 d({cos(angleOffset), sin(angleOffset)});
        const Point2 p = midPoint + radius * d;

        return p;
    }

    template<class C>
    inline
    bool
    ProblemEquilateralTriangles<C>::isTriangleInUnitSquare(const Triangle & triangle) {
        return (0 <= triangle.point1_[0] && triangle.point1_[0] <= 1) &&
                (0 <= triangle.point1_[1] && triangle.point1_[1] <= 1) &&
                (0 <= triangle.point2_[0] && triangle.point2_[0] <= 1) &&
                (0 <= triangle.point2_[1] && triangle.point2_[1] <= 1) &&
                (0 <= triangle.point3_[0] && triangle.point3_[0] <= 1) &&
                (0 <= triangle.point3_[1] && triangle.point3_[1] <= 1);
    }

    template<class C>
    inline
    bool
    ProblemEquilateralTriangles<C>::isPointInUnitSquare(const Point2 & p) {
        return (0 <= p[0] && p[0] <= 1) && (0 <= p[1] && p[1] <= 1);
    }

} // namespace dssc
