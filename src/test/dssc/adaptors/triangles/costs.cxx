
#include "dssc/equilateral-triangles/problem.hxx"
#include "dssc/equilateral-triangles/problem-adaptor.hxx"
#include "dssc/cost-adjustment-adaptors/problem-rescaling-adaptor.hxx"

typedef dssc::ProblemEquilateralTriangles<double> Problem;
typedef dssc::ProblemEquilateralTriangles<double>::Point2 Point;
typedef dssc::ProblemEquilateralTriangles<double>::Triangle Triangle;

template<dssc::AngleCostMode ANGLE_COST_MODE>
using ProblemAdaptor = dssc::ProblemEquilateralTrianglesAdaptor<double, ANGLE_COST_MODE>;

typedef dssc::CostRescalingAdaptor<double, dssc::RescalingCostType::TRIPLES> Rescaler;


inline void test(bool const condition) {
    if (!condition) throw std::logic_error("test failed.");
}

inline bool isClose(double value, double compareValue, double tolerance = 0.001) {
    return (abs(value - compareValue) <= tolerance);
}

inline
Triangle
getTriangle(double alpha, double beta, double gamma){
    Point const p0({0, 0});
    Point const p1({1, 0});

    double x = 1/(1 + tan(alpha)/tan(beta));
    double y = 1 / (1 / tan(alpha) + 1 / tan(beta));

    Point const p2({x, y});
    Triangle triangle(p0, p1, p2);

    return triangle;
}

inline
Triangle
getDegenerateTriangle(double x){
    Triangle triangle({0, 0}, {1, 0}, {x, 0});
    return triangle;
}

inline
double
degreesToRads(double degrees){
    return degrees*M_PI / 180;
}

void testAngleDifferenceSumAbs() {
    typedef ProblemAdaptor<dssc::AngleCostMode::SUM_ANGLE_DIFFERENCES_ABS> ProblemAdaptor;

    double const maximumAngleDifference = 4*M_PI / 3;

    {
        // right angle triangle

        // 30°
        double const &angleDifferenceAtCostZero = M_PI / 6;
        double const &tripletCostAtAngleDistanceZero = -1;
        double const &edgeCostAtDistanceZero = -1;
        double const &tripletMaxEdgeLengthThreshold = 0;
        double const &edgeLengthThreshold = 0;

        double const maximumPossibleCost = -tripletCostAtAngleDistanceZero / angleDifferenceAtCostZero*(maximumAngleDifference - angleDifferenceAtCostZero);

        Triangle const triangle = getTriangle(M_PI / 2, M_PI / 4, M_PI / 4);

        Problem problem(triangle.points());

        ProblemAdaptor adaptor(
                problem,
                angleDifferenceAtCostZero,
                tripletCostAtAngleDistanceZero,
                edgeCostAtDistanceZero,
                tripletMaxEdgeLengthThreshold,
                edgeLengthThreshold
        );

        Rescaler rescaled(
                adaptor,
                -tripletCostAtAngleDistanceZero,
                maximumPossibleCost,
                tripletCostAtAngleDistanceZero,
                tripletCostAtAngleDistanceZero
        );

        // angle difference is here 30 + 2*15 = 60
        // expected positive: 1*(240 - 60)/ (240 - 30)
        double const expected = (60. - 30) / (240 - 30);
        double const actual = rescaled.costOfTriple({2, 1, 0});

        test(isClose(actual, expected));
        test(rescaled.costOfEdge({1, 0}) == 0);
        test(rescaled.costOfEdge({2, 1}) == 0);
        test(rescaled.costOfEdge({2, 0}) == 0);
    }

    {
        // 30°
        double const &angleDifferenceAtCostZero = M_PI / 6;
        double const &tripletCostAtAngleDistanceZero = -1;
        double const &edgeCostAtDistanceZero = -1;
        double const &tripletMaxEdgeLengthThreshold = 0;
        double const &edgeLengthThreshold = 0;

        double const maximumPossibleCost = -tripletCostAtAngleDistanceZero / angleDifferenceAtCostZero*(maximumAngleDifference - angleDifferenceAtCostZero);

        Triangle const triangle = getTriangle(M_PI / 3, M_PI / 3, M_PI / 3);


        Problem problem(triangle.points());

        ProblemAdaptor adaptor(
                problem,
                angleDifferenceAtCostZero,
                tripletCostAtAngleDistanceZero,
                edgeCostAtDistanceZero,
                tripletMaxEdgeLengthThreshold,
                edgeLengthThreshold
        );

        Rescaler rescaler(
                adaptor,
                -tripletCostAtAngleDistanceZero,
                maximumPossibleCost,
                tripletCostAtAngleDistanceZero,
                tripletCostAtAngleDistanceZero
                );

        double const expected = -1;
        double const actual = rescaler.costOfTriple({2, 1, 0});

        test(isClose(actual, expected));
        test(adaptor.costOfEdge({1, 0}) == 0);
        test(adaptor.costOfEdge({2, 1}) == 0);
        test(adaptor.costOfEdge({2, 0}) == 0);
    }

    {
        // 30°
        double const &angleDifferenceAtCostZero = M_PI / 6;
        double const &tripletCostAtAngleDistanceZero = -1;
        double const &edgeCostAtDistanceZero = -1;
        double const &tripletMaxEdgeLengthThreshold = 0;
        double const &edgeLengthThreshold = 0;

        double const maximumPossibleCost = -tripletCostAtAngleDistanceZero / angleDifferenceAtCostZero*(maximumAngleDifference - angleDifferenceAtCostZero);

        Triangle const triangle = getDegenerateTriangle(0.2);


        Problem problem(triangle.points());

        ProblemAdaptor adaptor(
                problem,
                angleDifferenceAtCostZero,
                tripletCostAtAngleDistanceZero,
                edgeCostAtDistanceZero,
                tripletMaxEdgeLengthThreshold,
                edgeLengthThreshold
        );

        Rescaler rescaler(
                adaptor,
                -tripletCostAtAngleDistanceZero,
                maximumPossibleCost,
                tripletCostAtAngleDistanceZero,
                tripletCostAtAngleDistanceZero
        );

        double const expected = 1;
        double const actual = rescaler.costOfTriple({2, 1, 0});

        test(isClose(actual, expected));
        test(adaptor.costOfEdge({1, 0}) == 0);
        test(adaptor.costOfEdge({2, 1}) == 0);
        test(adaptor.costOfEdge({2, 0}) == 0);
    }

    {
        // 30°
        double const &angleDifferenceAtCostZero = M_PI / 6;
        double const &tripletCostAtAngleDistanceZero = -1;
        double const &edgeCostAtDistanceZero = -1;
        double const &tripletMaxEdgeLengthThreshold = 0;
        double const &edgeLengthThreshold = 0;

        double const maximumPossibleCost = -tripletCostAtAngleDistanceZero / angleDifferenceAtCostZero*(maximumAngleDifference - angleDifferenceAtCostZero);

        Triangle const triangle = getTriangle(
                degreesToRads(55),
                degreesToRads(55),
                degreesToRads(70)
                );

        //
        Problem problem(triangle.points());

        ProblemAdaptor adaptor(
                problem,
                angleDifferenceAtCostZero,
                tripletCostAtAngleDistanceZero,
                edgeCostAtDistanceZero,
                tripletMaxEdgeLengthThreshold,
                edgeLengthThreshold
        );

        Rescaler rescaler(
                adaptor,
                -tripletCostAtAngleDistanceZero,
                maximumPossibleCost,
                tripletCostAtAngleDistanceZero,
                tripletCostAtAngleDistanceZero
        );

        // expected delta = 20
        // expected cost = 1 /(30)*(20 - 30) = - 1/3
        double const expected = - 1.0 / 3;
        double const actual = rescaler.costOfTriple({2, 1, 0});

        test(isClose(actual, expected));
        test(adaptor.costOfEdge({1, 0}) == 0);
        test(adaptor.costOfEdge({2, 1}) == 0);
        test(adaptor.costOfEdge({2, 0}) == 0);
    }
}


int main() {
    testAngleDifferenceSumAbs();

    return 0;
}
