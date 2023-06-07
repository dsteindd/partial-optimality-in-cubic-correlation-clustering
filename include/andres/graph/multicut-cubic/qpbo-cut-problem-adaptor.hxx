#pragma once
#ifndef ANDRES_GRAPH_QPBO_ADAPTOR_HXX
#define ANDRES_GRAPH_QPBO_ADAPTOR_HXX

#include <cstddef>
#include <vector>
#include <set>
#include <array>

#include "problem.hxx"
#include "andres/graph/qpbo.hxx"
#include "andres/graph/qpbo-problem.hxx"


namespace andres {
    namespace graph {
        namespace multicut_cubic {

            template<class T = double>
            class QPBOCutProblemAdaptor : public QPBOProblem<T> {
            public:
                typedef T Rational;

                template<class PROBLEM_VIEW>
                QPBOCutProblemAdaptor(PROBLEM_VIEW const &, std::vector<size_t> const &, std::vector<size_t> const &);
                template<class PROBLEM_VIEW>
                QPBOCutProblemAdaptor(PROBLEM_VIEW const &);

            private:
                void persistOneVertex(size_t const &);
                void persistZeroVertex(size_t const &);
                void moveToEnd(size_t const &);
                std::vector<size_t> oneVertexConstraints_;
                std::vector<size_t> zeroVertexConstraints_;
                std::map<size_t, size_t> vertexLabels_;

            };


            template<class T>
            template<class PROBLEM_VIEW>
            inline
            QPBOCutProblemAdaptor<T>::QPBOCutProblemAdaptor(
                    PROBLEM_VIEW const & problem,
                    std::vector<size_t> const & oneVertexConstraints,
                    std::vector<size_t> const & zeroVertexConstraints
            ) : QPBOProblem<T>(
                    problem.numberOfElements()
                ),
                oneVertexConstraints_(oneVertexConstraints),
                zeroVertexConstraints_(zeroVertexConstraints),
                vertexLabels_()
            {
                    for (size_t i = 0; i < problem.numberOfElements(); ++i){
                        vertexLabels_[i] = i;
                    }


                // first do adjustments without taking into consideration the constraints

                (*this).constant_ = problem.costConstant();
                for (size_t i = 0; i < problem.numberOfElements(); ++i){
                    for (size_t j = 0; j < i; ++j){
                        for (size_t k = 0; k < j; ++k){
                            Rational tripleCost = problem.costOfTriple({i, j, k});
                            (*this).edgeCosts_[(*this).edgeIndex(i, j)] -= tripleCost;
                            (*this).edgeCosts_[(*this).edgeIndex(i, k)] -= tripleCost;
                            (*this).edgeCosts_[(*this).edgeIndex(j, k)] -= tripleCost;

                            (*this).vertexCosts_[i] += tripleCost;
                            (*this).vertexCosts_[j] += tripleCost;
                            (*this).vertexCosts_[k] += tripleCost;
                        }
                    }
                }

                for (size_t i = 0; i < problem.numberOfElements(); ++i){
                    for (size_t j = 0; j < i; ++j){
                        Rational edgeCost = problem.costOfEdge({i, j});
                        (*this).edgeCosts_[(*this).edgeIndex(i, j)] += -2*edgeCost;

                        (*this).vertexCosts_[i] += edgeCost;
                        (*this).vertexCosts_[j] += edgeCost;
                    }
                }

                // add constraints
                while (!oneVertexConstraints_.empty()){
                    const size_t nextOneVertex = oneVertexConstraints_[0];
                    persistOneVertex(nextOneVertex);
                    oneVertexConstraints_.erase(oneVertexConstraints_.begin());
                }

                while (!zeroVertexConstraints_.empty()){
                    const size_t nextZeroVertex = zeroVertexConstraints_[0];
                    persistZeroVertex(nextZeroVertex);
                    zeroVertexConstraints_.erase(zeroVertexConstraints_.begin());
                }

            }

            template<class T>
            template<class PROBLEM_VIEW>
            inline
            QPBOCutProblemAdaptor<T>::QPBOCutProblemAdaptor(
                    PROBLEM_VIEW const & problem
            ) : QPBOCutProblemAdaptor(problem, std::vector<size_t>(), std::vector<size_t>()){}


            template<class T>
            void QPBOCutProblemAdaptor<T>::persistOneVertex(const size_t & u) {
                moveToEnd(u);

                size_t last = (*this).numberOfVertices_ - 1;

                // adjust costs
                (*this).constant_ += (*this).vertexCosts_[last];

                for (size_t i = 0; i < last; ++i){
                    (*this).vertexCosts_[i] += (*this).edgeCosts_[(*this).edgeIndex(last, i)];
                }

                // resize costs
                (*this).numberOfVertices_--;
                (*this).vertexCosts_.resize((*this).numberOfVertices_);
                (*this).edgeCosts_.resize((*this).numberOfVertices_ * ((*this).numberOfVertices_ - 1)/2);
            }



            template<class T>
            void QPBOCutProblemAdaptor<T>::persistZeroVertex(const size_t & v) {
                moveToEnd(v);

                // resize costs
                (*this).numberOfVertices_--;
                (*this).vertexCosts_.resize((*this).numberOfVertices_);
                (*this).edgeCosts_.resize((*this).numberOfVertices_ * ((*this).numberOfVertices_ - 1)/2);
            }

            template<class T>
            void QPBOCutProblemAdaptor<T>::moveToEnd(const size_t & u) {
                // moved to end by simply switching last and u
                const size_t last = (*this).numberOfVertices_ - 1;

                vertexLabels_[u] = vertexLabels_[last];

                // adjust costs
                // c[u] = c[last]
                Rational oneVertexCost = (*this).vertexCosts_[u];
                (*this).vertexCosts_[u] = (*this).vertexCosts_[last];
                (*this).vertexCosts_[last] = oneVertexCost;

                // c[u, v] = c[last, v] for all v = 0... n - 2
                for (size_t v = 0; v < u; ++v){
                    if (u != v){
                        Rational uvEdgeCost = (*this).edgeCosts_[(*this).edgeIndex(u, v)];
                        (*this).edgeCosts_[(*this).edgeIndex(u, v)] = (*this).edgeCosts_[(*this).edgeIndex(last, v)];
                        (*this).edgeCosts_[(*this).edgeIndex(last, v)] = uvEdgeCost;
                    }
                }
                for (size_t v = u+1; v < (*this).numberOfVertices_ - 1; ++v){
                    if (u != v){
                        Rational vuEdgeCost = (*this).edgeCosts_[(*this).edgeIndex(v, u)];
                        (*this).edgeCosts_[(*this).edgeIndex(v, u)] = (*this).edgeCosts_[(*this).edgeIndex(last, v)];
                        (*this).edgeCosts_[(*this).edgeIndex(last, v)] = vuEdgeCost;
                    }
                }

                // adjust constraint vertices, each vertex which is higher than u needs to be decremented
                // last is in any of the constraint vertices, then relabel to u
                for (size_t i = 0; i < (*this).oneVertexConstraints_.size(); ++i){
                    if ((*this).oneVertexConstraints_[i] == last){
                        (*this).oneVertexConstraints_[i] = u;
                    }
                }
                for (size_t i = 0; i < (*this).zeroVertexConstraints_.size(); ++i){
                    if ((*this).zeroVertexConstraints_[i] == last){
                        (*this).zeroVertexConstraints_[i] = u;
                    }
                }
            }
        }
    }
}


#endif
