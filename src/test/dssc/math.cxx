#include <stdexcept>

#include "dssc/math.hxx"

inline void test(const bool condition) {
    if(!condition) throw std::logic_error("test failed.");
}

void testFitting() {
    typedef dssc::Point<double, 3> Point;

    {
        Point p0 = {0.0, 0.0, 0.0};
        Point p1 = {1.0, 0.0, 0.0};
        Point p2 = {0.0, 1.0, 0.0};
        Point unitNormalVector = dssc::fitLinearTLS({p0, p1, p2});
        test(unitNormalVector[0] == 0.0 && unitNormalVector[1] == 0 && unitNormalVector[2] == 1.0);
        test(dssc::distanceFromHyperplane(p0, unitNormalVector) == 0.0);
        test(dssc::distanceFromHyperplane(p1, unitNormalVector) == 0.0);
        test(dssc::distanceFromHyperplane(p2, unitNormalVector) == 0.0);
    }
    {
        Point p0 = {0.0, 0.0, 0.5};
        Point p1 = {1.0, 0.0, 0.0};
        Point p2 = {0.0, 1.0, 0.0};
        Point unitNormalVector = dssc::fitLinearTLS({p0, p1, p2});

        test(unitNormalVector[0] == 0.0 && unitNormalVector[1] == 0.0 && unitNormalVector[2] == 1.0);
        test(dssc::distanceFromHyperplane(p0, unitNormalVector) == 0.5);
        test(dssc::distanceFromHyperplane(p1, unitNormalVector) == 0.0);
        test(dssc::distanceFromHyperplane(p2, unitNormalVector) == 0.0);
    }

    {
        Point p0 = {0.0, 0.0, 0.0};
        Point p1 = {0.0, 1.0, 0.0};
        Point p2 = {0.0, 0.0, 1.0};
        Point unitNormalVector = dssc::fitLinearTLS({p0, p1, p2});
        test(unitNormalVector[0] == 1.0 && unitNormalVector[1] == 0 && unitNormalVector[2] == 0.0);
        test(dssc::distanceFromHyperplane(p0, unitNormalVector) == 0.0);
        test(dssc::distanceFromHyperplane(p1, unitNormalVector) == 0.0);
        test(dssc::distanceFromHyperplane(p2, unitNormalVector) == 0.0);
    }
    {
        Point p0 = {0.5, 0.0, 0.0};
        Point p1 = {0.0, 1.0, 0.0};
        Point p2 = {0.0, 0.0, 1.0};
        Point unitNormalVector = dssc::fitLinearTLS({p0, p1, p2});

        test(unitNormalVector[0] == 1.0 && unitNormalVector[1] == 0.0 && unitNormalVector[2] == 0.0);
        test(dssc::distanceFromHyperplane(p0, unitNormalVector) == 0.5);
        test(dssc::distanceFromHyperplane(p1, unitNormalVector) == 0.0);
        test(dssc::distanceFromHyperplane(p2, unitNormalVector) == 0.0);
    }

    {
        Point p0 = {0.0, 0.0, 0.0};
        Point p1 = {1.0, 0.0, 0.0};
        Point p2 = {0.0, 0.0, 1.0};
        Point unitNormalVector = dssc::fitLinearTLS({p0, p1, p2});
        test(unitNormalVector[0] == 0.0 && unitNormalVector[1] == 1.0 && unitNormalVector[2] == 0.0);
        test(dssc::distanceFromHyperplane(p0, unitNormalVector) == 0.0);
        test(dssc::distanceFromHyperplane(p1, unitNormalVector) == 0.0);
        test(dssc::distanceFromHyperplane(p2, unitNormalVector) == 0.0);
    }
    {
        Point p0 = {0.0, 0.5, 0.0};
        Point p1 = {1.0, 0.0, 0.0};
        Point p2 = {0.0, 0.0, 1.0};
        Point unitNormalVector = dssc::fitLinearTLS({p0, p1, p2});

        test(unitNormalVector[0] == 0.0 && unitNormalVector[1] == 1.0 && unitNormalVector[2] == 0.0);
        test(dssc::distanceFromHyperplane(p0, unitNormalVector) == 0.5);
        test(dssc::distanceFromHyperplane(p1, unitNormalVector) == 0.0);
        test(dssc::distanceFromHyperplane(p2, unitNormalVector) == 0.0);
    }
}

int main() {
    testFitting();

    return 0;
}
