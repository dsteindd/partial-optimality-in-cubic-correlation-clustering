#pragma once
#ifndef DSSC_PERSISTENCY_RESULTS_HXX
#define DSSC_PERSISTENCY_RESULTS_HXX

#include <algorithm>
#include <array>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "dssc/math.hxx"
#include "dssc/latex.hxx"
#include "andres/graph/multicut-cubic/problem.hxx"
#include "andres/graph/multicut-cubic/persistency.hxx"
namespace dssc {
namespace persistency {

template<class PROBLEM_VIEW, class R>
void
printPersistencyResults(
    std::ostream & out,
    PROBLEM_VIEW const & problem,
    andres::graph::multicut_cubic::Persistency<R> const & persistency,
    std::map<std::string, std::string> const & configurationValues = std::map<std::string, std::string>()
) {
    out << "Problem Instance" << std::endl;
    out << std::string(20, '=') << std::endl;
    out << "numberOfNodes=" << problem.numberOfElements() << std::endl;

    for (auto & [key, value]: configurationValues){
        out << key << "=" << value << std::endl;
    }

    out << std::endl;

    out << "Persistency Results" << std::endl;
    out << std::string(20, '=') << std::endl;

    out << "independentSubProblems=" << persistency.subProblems().size() << std::endl;
    for (auto subProblem : persistency.subProblems()){
        out << "{";
        for (size_t n: subProblem){
            out << n << ",";
        }
        out << "}";
        out << std::endl;
    }


    out << "Zero Edges: " << persistency.edges0().size() << std::endl;
    out << "One Edges: " << persistency.edges1().size() << std::endl;
    out << "Zero Triplets: " << persistency.triples0().size() << std::endl;

    out << std::endl;

    out << "Durations" << std::endl;
    out << std::string(20, '=') << std::endl;

    out << "findIndependentProblems=" << persistency.totalDurationFindIndependentProblems().count() * 1000 << "ms" << std::endl;
    out << "edgeCut=" << persistency.totalDurationEdgeCutCriterion().count() * 1000 << "ms" << std::endl;
    out << "tripletCut=" << persistency.totalDurationTripletCutCriterion().count() * 1000 << "ms" << std::endl;
    out << "edgeJoin=" << persistency.totalDurationEdgeJoinCriterion().count() * 1000 << "ms" << std::endl;
    out << "tripletEdgeJoin=" << persistency.totalDurationTripletEdgeJoinCriterion().count() * 1000 << "ms" << std::endl;
    out << "tripletJoin=" << persistency.totalDurationTripletJoinCriterion().count() * 1000 << "ms" << std::endl;
    out << "edgeSubgraph=" << persistency.totalDurationEdgeSubgraphCriterion().count() * 1000 << "ms" << std::endl;
    out << "tripletSubgraph=" << persistency.totalDurationTripletSubgraphCriterion().count() * 1000 << "ms" << std::endl;

    double totalRuntime = persistency.totalDurationFindIndependentProblems().count()
            + persistency.totalDurationEdgeCutCriterion().count()
            + persistency.totalDurationTripletCutCriterion().count()
            + persistency.totalDurationEdgeJoinCriterion().count()
            + persistency.totalDurationTripletEdgeJoinCriterion().count()
            + persistency.totalDurationTripletJoinCriterion().count()
            + persistency.totalDurationEdgeSubgraphCriterion().count()
            + persistency.totalDurationTripletSubgraphCriterion().count();

    out << "total=" << totalRuntime << "s" << std::endl;

}

template<class PROBLEM_VIEW, class R>
void
printPersistencyResults(
    std::string const & fileName,
    PROBLEM_VIEW const & problem,
    andres::graph::multicut_cubic::Persistency<R> const & persistency,
    std::map<std::string, std::string> const & configurationValues = std::map<std::string, std::string>()
) {
    std::ofstream fileStream;
    fileStream.open(fileName, std::ios::out);
    if(fileStream) {
        printPersistencyResults(fileStream, problem, persistency, configurationValues);
    }
    fileStream.close();
}

}
}

#endif
