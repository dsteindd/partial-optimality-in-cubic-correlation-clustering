//
// Created by dstein on 13.01.23.
//

#ifndef CUBIC_MULTICUT_NEGATIVE_PART_COST_PROBLEM_ADAPTOR_HXX
#define CUBIC_MULTICUT_NEGATIVE_PART_COST_PROBLEM_ADAPTOR_HXX

#include <cassert>
#include <array>
#include <vector>
#include <functional>
#include <algorithm>
#include <map>
#include <set>
#include "problem.hxx"

namespace andres{
    namespace graph{
        namespace multicut_cubic{

            template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES = false>
            class NegativePartCostProblemAdaptor{
            public:
                typedef T Rational;
                typedef std::array<size_t, 2> IndexPair;
                typedef std::array<size_t, 3> IndexTriple;

                NegativePartCostProblemAdaptor(Problem<T> const &);

                size_t numberOfElements() const;
                Rational costOfEdge(IndexPair const &) const;
                Rational costOfTriple(IndexTriple const &) const;
                Rational costConstant() const;

            private:
                Problem<T, ASSUME_STRICTLY_DECREASING_TRIPLES> const & problemAdaptee_;
            };

            template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
            NegativePartCostProblemAdaptor<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::NegativePartCostProblemAdaptor(
                    const Problem<T> & problem):
                    problemAdaptee_(problem){}

            template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
            inline
            size_t
            NegativePartCostProblemAdaptor<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::numberOfElements() const {
                return problemAdaptee_.numberOfElements();
            }

            template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
            inline
            typename NegativePartCostProblemAdaptor<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::Rational
            NegativePartCostProblemAdaptor<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::costOfEdge(IndexPair const & indexPair) const {
                Rational const cost = problemAdaptee_.costOfEdge(indexPair);
                if (cost < 0){
                    return -cost;
                }
                return 0;
            }

            template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
            inline
            typename NegativePartCostProblemAdaptor<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::Rational
            NegativePartCostProblemAdaptor<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::costOfTriple(const IndexTriple & indexTriple) const {
                Rational const cost = problemAdaptee_.costOfTriple(indexTriple);
                if (cost < 0){
                    return -cost;
                }
                return 0;
            }

            template<class T, bool ASSUME_STRICTLY_DECREASING_TRIPLES>
            inline
            typename NegativePartCostProblemAdaptor<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::Rational
            NegativePartCostProblemAdaptor<T, ASSUME_STRICTLY_DECREASING_TRIPLES>::costConstant() const {
                return 0;
            }


        }
    }
}

#endif //CUBIC_MULTICUT_ABSOLUTE_COST_PROBLEM_ADAPTOR_HXX
