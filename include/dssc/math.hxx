#pragma once
#ifndef ANDRES_DSSC_FIT_HXX
#define ANDRES_DSSC_FIT_HXX

#include <cmath>
#include <array>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Jacobi>

namespace dssc {

template<class C = double, size_t D = 2>
using Point = std::array<C, D>;

template<class C, size_t D>
inline
Point<C, D>
operator+(
    Point<C, D> const & a,
    Point<C, D> const & b
) {
    Point<C, D> c;
    for(size_t j = 0; j < D; ++j) {
        c[j] = a[j] + b[j];
    }
    return c;
}

template<class C, size_t D>
inline
Point<C, D>
operator-(
    Point<C, D> const & a,
    Point<C, D> const & b
) {
    Point<C, D> c;
    for(size_t j = 0; j < D; ++j) {
        c[j] = a[j] - b[j];
    }
    return c;
}

template<class C, size_t D>
inline
Point<C, D>
operator/(
    Point<C, D> const & a,
    C const alpha
) {
    Point<C, D> b;
    for(size_t j = 0; j < D; ++j) {
        b[j] = a[j] / alpha;
    }
    return b;
}

template<class C, size_t D>
inline
Point<C, D>
operator*(
    Point<C, D> const & a,
    C const alpha
) {
    Point<C, D> b;
    for(size_t j = 0; j < D; ++j) {
        b[j] = a[j] * alpha;
    }
    return b;
}

template<class C, size_t D>
inline
Point<C, D>
operator*(
    C const alpha,
    Point<C, D> const & a
) {
    return operator*(a, alpha);
}

template<class C, size_t D>
inline
C
scalarProduct(
    Point<C, D> const & a,
    Point<C, D> const & b
){
    C result = C();
    for(size_t j = 0; j < D; ++j) {
        result += a[j] * b[j];
    }
    return result;
}

template<class C>
inline
Point<C, 3>
crossProduct(
    Point<C, 3> const & a,
    Point<C, 3> const & b
){
    return {
        a[1]*b[2] - a[2]*b[1],
        a[2]*b[0] - a[0]*b[2],
        a[0]*b[1] - a[1]*b[0]
    };
}

template<class C, size_t D>
inline
C
norm(
    Point<C, D> const & p
) {
    return std::sqrt(scalarProduct(p, p));
}

template<class C, size_t D>
inline
Point<C, D>
normalize(
    Point<C, D> const & p
){
    return p / norm(p);
}

template<class C, size_t D>
inline
C
distance(
    Point<C, D> const & p,
    Point<C, D> const & q
) {
    return norm(p - q);
}

template<class C, size_t D>
inline
C
angleBetween(
        Point<C, D> const & p,
        Point<C, D> const & q
             ){
    return acos(scalarProduct(p, q) / (norm(p)* norm(q)));
}

template<class C, size_t D>
inline
C
distanceFromHyperplane(
    Point<C, D> const & point,
    Point<C, D> const & unitNormalVectorHyperplane,
    C affineOffset = C()
) {
    C const q = scalarProduct(unitNormalVectorHyperplane, point) - affineOffset;
    return q < 0 ? -q : q;
}

template<class C, size_t D>
inline
Point<C, D>
footOfPerpendicularOnHyperplane(
    Point<C, D> const & point,
    Point<C, D> const & unitNormalVectorHyperplane,
    C affineOffset = C()
) {
    C const alpha = affineOffset - scalarProduct(unitNormalVectorHyperplane, point);
    return point + alpha * unitNormalVectorHyperplane;
}

template<class C>
inline
std::vector<Point<C, 2>>
intersectionCircleLine(
    C radius,
    Point<C, 2> const & unitNormalVectorLine,
    C affineOffsetLine = C()
) {
    C const a = unitNormalVectorLine[0];
    C const b = unitNormalVectorLine[1];
    C const c = -affineOffsetLine;

    if(affineOffsetLine > radius) {
        return std::vector<Point<C, 2>>(); // no intersection
    }
    else {
        C const alpha = std::sqrt(radius*radius - c*c);
        C const x0 = -a * c;
        C const y0 = -b * c;
        return {
            {x0 + b * alpha, y0 - a * alpha},
            {x0 - b * alpha, y0 + a * alpha}
        };
    }
}

template<class POINT_ITERATOR>
inline
Point<double, 3>
fitLinearTLS(
    POINT_ITERATOR it,
    POINT_ITERATOR end
) {
    typedef Eigen::Matrix<double, 3, 1> Vector;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 3> Matrix;
    typedef Eigen::JacobiSVD<Matrix> SVD;
    // typedef Eigen::BDCSVD<Matrix> SVD;

    size_t const numberOfPoints = std::distance(it, end);
    Matrix pointMatrix(numberOfPoints, 3);
    for(size_t j = 0; j < numberOfPoints; ++j, ++it) {
        for(size_t d = 0; d < 3; ++d) {
            pointMatrix(j, d) = (*it)[d];
        }
    }

    SVD svd(pointMatrix, Eigen::ComputeFullU | Eigen::ComputeFullV); // perform SVD
    Vector const unitNormalVector = svd.matrixV().col(2); // rigth singular vector wrt. smallest singular value
    return {unitNormalVector(0), unitNormalVector(1), unitNormalVector(2)};
}

template<class C>
inline
Point<C, 3>
fitLinearTLS(
    std::initializer_list<Point<C, 3>> const & points
) {
    return fitLinearTLS(points.begin(), points.end());
}

template<class POINT_ITERATOR>
inline
void
fitAffineTLS(
    POINT_ITERATOR it,
    POINT_ITERATOR end,
    Point<double, 2> & unitNormalVector,
    double & offset
) {
    typedef Eigen::Matrix<double, 3, 1> Vector;
    typedef Eigen::Matrix<double, Eigen::Dynamic, 3> Matrix;
    typedef Eigen::JacobiSVD<Matrix> SVD;
    // typedef Eigen::BDCSVD<Matrix> SVD;

    size_t const numberOfPoints = std::distance(it, end);
    Matrix pointMatrix(numberOfPoints, 3);
    for(size_t j = 0; j < numberOfPoints; ++j, ++it) {
        for(size_t d = 0; d < 2; ++d) {
            pointMatrix(j, d) = (*it)[d];
        }
        pointMatrix(j, 2) = 1.0; // affine
    }

    SVD svd(pointMatrix, Eigen::ComputeFullU | Eigen::ComputeFullV); // perform SVD
    Vector const v = svd.matrixV().col(2); // rigth singular vector wrt. smallest singular value
    double const factor = std::sqrt(v(0)*v(0) + v(1)*v(1));
    unitNormalVector[0] = v(0) / factor;
    unitNormalVector[1] = v(1) / factor;
    offset = -v(2) / factor;
}

template<class C>
inline
void
fitAffineTLS(
    std::initializer_list<Point<C, 2>> const & points,
    Point<double, 2> & unitNormalVector,
    double & offset
) {
    fitAffineTLS(points.begin(), points.end(), unitNormalVector, offset);
}

} // namespace dssc

#endif
