#pragma once
#ifndef ANDRES_DSSC_FIT_HXX
#define ANDRES_DSSC_FIT_HXX

#include <cmath>
#include <array>
#include <vector>

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

} // namespace dssc

#endif
