#pragma once
#ifndef ANDRES_GRAPH_QPBO_PROBLEM_HXX
#define ANDRES_GRAPH_QPBO_PROBLEM_HXX

#include "cassert"

namespace andres {
    namespace graph {

        template<class T = double>
        class QPBOProblem {
        public:
            typedef T Rational;

            QPBOProblem(size_t const);

            Rational &vertexCost(size_t const &);

            Rational vertexCost(size_t const &) const;

            Rational &edgeCost(size_t const &, size_t const &);

            Rational edgeCost(size_t const &, size_t const &) const;

            Rational &constant();

            Rational constant() const;

            size_t numberOfVertices() const;

        protected:
            std::vector<Rational> edgeCosts_;
            std::vector<Rational> vertexCosts_;
            size_t numberOfVertices_;
            Rational constant_;

            size_t edgeIndex(size_t const &, size_t const &);
        };

        template<class T>
        inline
        QPBOProblem<T>::QPBOProblem(const size_t numberOfVertices)
                :
                numberOfVertices_(numberOfVertices),
                edgeCosts_(numberOfVertices * (numberOfVertices - 1) / 2),
                vertexCosts_(numberOfVertices),
                constant_() {

        }

        template<class T>
        inline
        typename QPBOProblem<T>::Rational &
        QPBOProblem<T>::vertexCost(const size_t &u) {
            return vertexCosts_[u];
        }

        template<class T>
        inline
        typename QPBOProblem<T>::Rational
        QPBOProblem<T>::vertexCost(const size_t &u) const {
            return vertexCosts_[u];
        }

        template<class T>
        inline
        typename QPBOProblem<T>::Rational &
        QPBOProblem<T>::edgeCost(const size_t &u, const size_t &v) {
            return edgeCosts_[edgeIndex(u, v)];
        }

        template<class T>
        inline
        typename QPBOProblem<T>::Rational
        QPBOProblem<T>::edgeCost(const size_t &u, const size_t &v) const {
            return edgeCosts_[edgeIndex(u, v)];
        }

        template<class T>
        inline
        typename QPBOProblem<T>::Rational &
        QPBOProblem<T>::constant() {
            return constant_;
        }

        template<class T>
        inline
        typename QPBOProblem<T>::Rational
        QPBOProblem<T>::constant() const {
            return constant_;
        }

        template<class T>
        size_t QPBOProblem<T>::edgeIndex(const size_t &u, const size_t &v) {
            assert(u > v);
            return u * (u - 1) / 2 + v;
        }

        template<class T>
        size_t QPBOProblem<T>::numberOfVertices() const {
            return numberOfVertices_;
        }
    }
}


#endif
