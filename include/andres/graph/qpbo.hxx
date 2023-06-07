#pragma once
#ifndef ANDRES_GRAPH_QPBO_HXX
#define ANDRES_GRAPH_QPBO_HXX
#define BOOST_DISABLE_ASSERTS

#include <functional>
#include <algorithm>
#include <stdexcept>
#include <set>
#include <vector>

#include "qpbo-problem.hxx"
#include <boost/graph/push_relabel_max_flow.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/read_dimacs.hpp>
#include <boost/assert.hpp>


namespace andres {
    namespace graph {
        namespace qpbo {

            using namespace boost;

            enum struct MaxFlowSolver {
                BOOST_PUSH_RELABEL
            };

            typedef adjacency_list_traits<vecS, vecS, directedS> Traits;

            typedef adjacency_list<vecS, vecS, directedS, no_property,
                    property<edge_capacity_t, double,
                            property<edge_residual_capacity_t, double,
                                    property<edge_reverse_t, Traits::edge_descriptor> > >
            > DoubleGraph;
            typedef property_map<DoubleGraph, edge_capacity_t>::type DoubleCapacityMap;
            typedef property_map<DoubleGraph, edge_reverse_t>::type DoubleReverseMap;


            typedef Traits::vertex_descriptor VertexDescriptor;
            typedef Traits::edge_descriptor EdgeDescriptor;


            template<class T, MaxFlowSolver Solver = MaxFlowSolver::BOOST_PUSH_RELABEL>
            inline
            void solveSubmodularQPBO(
                    QPBOProblem<T> &problem,
                    T &optimalValue
            ) {
                if (Solver == MaxFlowSolver::BOOST_PUSH_RELABEL) {
                    typedef T Rational;

                    size_t const numberOfVertices = problem.numberOfVertices() + 2;

                    DoubleGraph g(numberOfVertices);

                    DoubleCapacityMap capacity = get(edge_capacity, g);
                    DoubleReverseMap rev = get(edge_reverse, g);

                    VertexDescriptor s, t;
                    s = vertex(0, g);
                    t = vertex(numberOfVertices - 1, g);

                    Rational constant = problem.constant();

                    // prevVertex to prevVertex costs
                    for (size_t p = 1; p < numberOfVertices - 1; ++p) {
                        for (size_t q = p + 1; q < numberOfVertices - 1; ++q) {
                            Rational edgeCost = -problem.edgeCost(q - 1, p - 1);
                            assert(edgeCost >= 0);

                            if (edgeCost > 0) {
                                VertexDescriptor p_descriptor = vertex(p, g);
                                VertexDescriptor q_descriptor = vertex(q, g);

                                EdgeDescriptor e, er;
                                bool in1, in2;

                                tie(e, in1) = add_edge(p_descriptor, q_descriptor, g);
                                tie(er, in2) = add_edge(q_descriptor, p_descriptor, g);

                                capacity[e] = edgeCost;
                                capacity[er] = 0;

                                rev[e] = er;
                                rev[er] = e;
                            }
                        }
                    }

                    for (size_t p = 1; p < numberOfVertices - 1; ++p) {
                        Rational adjustedCost = problem.vertexCost(p - 1);
                        for (size_t q = p + 1; q < numberOfVertices - 1; ++q) {
                            adjustedCost += problem.edgeCost(q - 1, p - 1);
                        }

                        if (adjustedCost < 0) {
                            // source connection
                            VertexDescriptor p_descriptor = vertex(p, g);

                            EdgeDescriptor e, er;
                            bool in1, in2;

                            tie(e, in1) = add_edge(s, p_descriptor, g);
                            tie(er, in2) = add_edge(p_descriptor, s, g);
                            capacity[e] = -adjustedCost;
                            capacity[er] = 0;

                            rev[e] = er;
                            rev[er] = e;

                            constant += adjustedCost;

                        } else if (adjustedCost > 0) {
                            // sink connection
                            VertexDescriptor p_descriptor = vertex(p, g);

                            EdgeDescriptor e, er;
                            bool in1, in2;

                            tie(e, in1) = add_edge(p_descriptor, t, g);
                            tie(er, in2) = add_edge(t, p_descriptor, g);

                            capacity[e] = adjustedCost;
                            capacity[er] = 0;

                            rev[e] = er;
                            rev[er] = e;
                        }
                    }

                    Rational const flow = push_relabel_max_flow(g, s, t);

                    optimalValue = constant + flow;
                }
                else {
                    throw std::runtime_error("unknown solver.");
                }

            }
        }
    }
}

#endif
