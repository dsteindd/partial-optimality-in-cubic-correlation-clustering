#pragma once
#ifndef DSSC_RESULT_WRITER_HXX
#define DSSC_RESULT_WRITER_HXX

#include <cassert>
#include <iostream>
#include <fstream>


#include "dssc/math.hxx"
#include "random"
#include "andres/graph/multicut-cubic/persistency.hxx"

namespace dssc {

/// Gaussian Mixture Model
    template<class R>
    class PersistencyResultWriter {
    public:
        typedef R Rational;
        typedef andres::graph::multicut_cubic::Persistency<Rational> Persistency;

        // constructor for only triplet cost
        PersistencyResultWriter(
                std::string const &,
                Persistency const &,
                std::map<std::string, std::string> const &
        );

        void writeResults();
        void collectResults(std::string const &);
        void writeRuntimes();
        void close();

    private:
        Persistency const & persistency_;
        std::ofstream* out_;
        std::vector<std::string> subsections_;
        std::vector<std::string> oneEdges_;
        std::vector<std::string> zeroEdges_;
        std::vector<std::string> zeroTriples_;
        std::vector<std::string> remainingNodes_;
        std::vector<std::string> remainingVariables_;
    };

    template<class R>
    inline
    PersistencyResultWriter<R>::PersistencyResultWriter(
            std::string const & fileName,
            Persistency const & persistency,
            std::map<std::string, std::string> const & configurationValues
    )
            :
           persistency_(persistency),
           subsections_(),
           oneEdges_(),
           zeroEdges_(),
           zeroTriples_(),
           remainingNodes_(),
           remainingVariables_()
    {
        std::ofstream* fileStream;
        fileStream = new std::ofstream(fileName, std::ios::out);

        if (fileStream->is_open()){
            out_ = fileStream;

            *out_ << "Instance Parameters" << std::endl << std::string (20, '=') << std::endl;

            for (auto & [key, value]: configurationValues){
                *out_ << key << "=" << value << std::endl;
            }

            *out_ << std::endl;
        }
    }

    template<class R>
    inline
    void PersistencyResultWriter<R>::collectResults(const std::string & subsection) {
        subsections_.push_back(subsection);
        oneEdges_.push_back(std::to_string(persistency_.edges1().size()));
        zeroEdges_.push_back(std::to_string(persistency_.edges0().size()));
        zeroTriples_.push_back(std::to_string(persistency_.triples0().size()));
        remainingNodes_.push_back(std::to_string(persistency_.remainingNodes()));
        remainingVariables_.push_back(std::to_string(persistency_.remainingVariables()));
    }

    template<class R>
    inline
    void PersistencyResultWriter<R>::writeResults() {
        *out_ << "Results" << std::endl << std::string (20, '=') << std::endl;
        *out_ << "| criterion | oneEdges | zeroEdges | zeroTriples | remainingNodes | remainingVariables" << std::endl;

        for (size_t index = 0; index < this->subsections_.size(); ++index){
            *out_ <<  "| " << this->subsections_[index];
            *out_ <<  "| " << this->oneEdges_[index];
            *out_ <<  "| " << this->zeroEdges_[index];
            *out_ <<  "| " << this->zeroTriples_[index];
            *out_ << "| " << this->remainingNodes_[index];
            *out_ << "| " << this->remainingVariables_[index] << "|";
            *out_ << std::endl;
        }

        *out_ << std::endl;
    }

    template<class R>
    inline
    void PersistencyResultWriter<R>::writeRuntimes() {
        *out_ << "Durations" << std::endl;
        *out_ << std::string(20, '=') << std::endl;

        *out_ << "findIndependentProblems=" << persistency_.totalDurationFindIndependentProblems().count() * 1000 << "ms" << std::endl;
        *out_ << "edgeCut=" << persistency_.totalDurationEdgeCutCriterion().count() * 1000 << "ms" << std::endl;
        *out_ << "tripletCut=" << persistency_.totalDurationTripletCutCriterion().count() * 1000 << "ms" << std::endl;
        *out_ << "edgeJoin=" << persistency_.totalDurationEdgeJoinCriterion().count() * 1000 << "ms" << std::endl;
        *out_ << "tripletEdgeJoin=" << persistency_.totalDurationTripletEdgeJoinCriterion().count() * 1000 << "ms" << std::endl;
        *out_ << "tripletJoin=" << persistency_.totalDurationTripletJoinCriterion().count() * 1000 << "ms" << std::endl;
        *out_ << "edgeSubgraph=" << persistency_.totalDurationEdgeSubgraphCriterion().count() * 1000 << "ms" << std::endl;
        *out_ << "tripletSubgraph=" << persistency_.totalDurationTripletSubgraphCriterion().count() * 1000 << "ms" << std::endl;
        *out_ << "subsetJoin=" << persistency_.totalDurationSubsetJoinCriterion().count() * 1000 << "ms" << std::endl;

        double totalRuntime = persistency_.totalDurationFindIndependentProblems().count()
                              + persistency_.totalDurationEdgeCutCriterion().count()
                              + persistency_.totalDurationTripletCutCriterion().count()
                              + persistency_.totalDurationEdgeJoinCriterion().count()
                              + persistency_.totalDurationTripletEdgeJoinCriterion().count()
                              + persistency_.totalDurationTripletJoinCriterion().count()
                              + persistency_.totalDurationEdgeSubgraphCriterion().count()
                              + persistency_.totalDurationTripletSubgraphCriterion().count()
                              + persistency_.totalDurationSubsetJoinCriterion().count();

        *out_ << "total=" << totalRuntime << "s" << std::endl;
    }

    template<class R>
    inline
    void PersistencyResultWriter<R>::close() {
        out_->close();
    }

}

#endif
