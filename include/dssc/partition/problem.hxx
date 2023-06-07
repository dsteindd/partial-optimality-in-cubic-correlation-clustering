#pragma once
#ifndef DSSC_PROBLEM_PARTITION_HXX
#define DSSC_PROBLEM_PARTITION_HXX

#include <cassert>

#include "dssc/math.hxx"
#include "random"

namespace dssc {

    enum struct CostMode {ONE_PLANE_LINEAR, ONE_PLANE_AFFINE, THREE_PLANES};

/// Gaussian Mixture Model
    template<class R, class D>
    class ProblemPartition {
    public:
        typedef R Rational;
        // D must implement the call operator with a random engine as argument
        typedef D Distribution;
        typedef std::array<size_t, 3> IndexTriple;
        typedef std::array<size_t, 2> IndexEdge;
        typedef std::vector<size_t > LabelContainer;

        // constructor for only triplet cost
        ProblemPartition(
                size_t const &,
                std::vector<Rational > const &,
                Distribution & = std::normal_distribution<Rational >(-1, 0.1),
                Distribution & = std::normal_distribution<Rational >(1, 0.1),
                size_t const & seed = 4242
        );
        // constructor for only triplet cost
        ProblemPartition(
                size_t const &,
                size_t const &,
                Distribution & = std::normal_distribution<Rational>(-1, 0.1),
                Distribution & = std::normal_distribution<Rational >(1, 0.1),
                size_t const & seed = 4242
        );
        ProblemPartition(
                std::vector<size_t> const &,
                Distribution & = std::normal_distribution<Rational>(-1, 0.1),
                Distribution & = std::normal_distribution<Rational >(1, 0.1),
                size_t const & seed = 4242
                );
        size_t numberOfElements() const;
        Rational costOfTriple(IndexTriple const &) const;
        Rational costOfEdge(IndexEdge const &) const;
        LabelContainer const & labels() const;

        // getter
        std::vector<size_t > const & edgeDrawStatistics() const;
        std::vector<size_t > const & tripletDrawStatistics() const;

    private:
        size_t indexOfTriple(IndexTriple const &) const;
        size_t indexOfTriple(size_t const, size_t const, size_t const) const;
        size_t indexOfEdge(IndexEdge const &) const;
        size_t indexOfEdge(size_t const, size_t const) const;
        Rational generateRandomNumber(size_t const);
        size_t computeDistributionIndex();

        void initCosts();

        size_t numberOfElements_;
        std::vector<Rational> costsOfEdges_;
        std::vector<Rational> costsOfTriples_;

        std::vector<Rational> probabilities_;

        std::default_random_engine engine_;
        std::uniform_real_distribution<Rational> uniform_;
        Distribution & intraDistribution_;
        Distribution & interDistribution_;

        std::vector<size_t > labels_;

    };

    // constructor for uniformly distributed cluster sizes with 'numberOfClusters' clusters
    template<class R, class D>
    inline
    ProblemPartition<R, D>::ProblemPartition(
            size_t const & numberOfElements,
            size_t const & numberOfClusters,
            Distribution & intraDistribution,
            Distribution & interDistribution,
            size_t const & seed
            ): ProblemPartition(
                    numberOfElements,
                    std::vector<Rational>(numberOfClusters, 1),
                    intraDistribution,
                    interDistribution,
                    seed
                    )
               {
                // init is done by called constructor
                }

            // constructor for specified cluster sizes
    template<class R, class D>
    inline
    ProblemPartition<R, D>::ProblemPartition(
            std::vector<size_t > const & clusterSizes,
            Distribution & intraDistribution,
            Distribution & interDistribution,
            size_t const & seed
    ):
       interDistribution_(interDistribution),
       intraDistribution_(intraDistribution),
       uniform_(0, 1),
       engine_(seed)
       {
           numberOfElements_ = 0;
           for (auto & clusterSize : clusterSizes){
               numberOfElements_ += clusterSize;
           }
           labels_ = std::vector<size_t >(numberOfElements_);
           costsOfEdges_ = std::vector<Rational >(numberOfElements_ * (numberOfElements_ - 1)/ 2);
           costsOfTriples_ = std::vector<Rational >(numberOfElements_ * (numberOfElements_ - 1) * (numberOfElements_ - 2)/ 6);

        size_t offset = 0;
        for (size_t clusterIndex = 0; clusterIndex < clusterSizes.size(); ++clusterIndex){
            const size_t clusterSize = clusterSizes[clusterIndex];

            for (size_t index = offset; index < offset + clusterSize; ++index){
                labels_[index] = clusterIndex;
            }

            offset += clusterSizes[clusterIndex];
        }

        initCosts();
    }

    // constructor for probabilistic cluster sizes
    template<class R, class D>
    inline
    ProblemPartition<R, D>::ProblemPartition(
            size_t const & numberOfElements,
            std::vector<Rational > const & probabilities,
            Distribution & intraDistribution,
            Distribution & interDistribution,
            size_t const & seed
    )
            :
            numberOfElements_(numberOfElements),
            costsOfEdges_(numberOfElements*(numberOfElements - 1)/2),
            costsOfTriples_(numberOfElements*(numberOfElements - 1)*(numberOfElements - 2)/6),
            probabilities_(probabilities),
            interDistribution_(interDistribution),
            intraDistribution_(intraDistribution),
            uniform_(0, 1),
            engine_(seed),
            labels_(numberOfElements)
    {
        // assert that all sigmas are positive,
        // assert that all probablities are positive
        for (auto prob: probabilities_) assert(prob > 0);

        // normalize probablities
        // calculate sum
        Rational probSum = 0;
        for (auto prob: probabilities_){
            probSum += prob;
        }

        for (size_t i = 0; i < probabilities_.size(); ++i){
            probabilities_[i] = probabilities_[i]/probSum;
        }

        // draw random cluster indices
        for (size_t index = 0; index < numberOfElements; ++index){
            labels_[index] = computeDistributionIndex();
        }

        initCosts();
    }

    template<class R, class D>
    inline
    size_t
    ProblemPartition<R, D>::computeDistributionIndex() {
        const Rational p = uniform_(engine_);

        // 2. find index of gaussian to use
        Rational currentSum = 0;
        for (size_t index = 0; index < probabilities_.size(); ++index){

            if (currentSum < p && p < currentSum + probabilities_[index]){
                return index;
            }

            currentSum += probabilities_[index];
        }
    }


    template<class R, class D>
    inline
    size_t
    ProblemPartition<R, D>::numberOfElements() const {
        return numberOfElements_;
    }

    template<class R, class D>
    inline
    R
    ProblemPartition<R, D>::costOfTriple(
            IndexTriple const & indexTriple
    ) const {

        return costsOfTriples_[indexOfTriple(indexTriple)];
    }

    template<class R, class D>
    inline
    R
    ProblemPartition<R, D>::costOfEdge(const IndexEdge & indexEdge) const {
        return costsOfEdges_[indexOfEdge(indexEdge)];
    }


    template<class R, class D>
    inline
    size_t
    ProblemPartition<R, D>::indexOfTriple(
            IndexTriple const & indexTriple
    ) const {
        assert(indexTriple[0] > indexTriple[1] && indexTriple[1] > indexTriple[2]);
        return indexOfTriple(indexTriple[0], indexTriple[1], indexTriple[2]);
    }

    template<class R, class D>
    inline
    size_t
    ProblemPartition<R, D>::indexOfEdge(
            IndexEdge const & indexEdge
    ) const {
        assert(indexEdge[0] > indexEdge[1]);
        return indexOfEdge(indexEdge[0], indexEdge[1]);

    }

    template<class R, class D>
    inline
    size_t
    ProblemPartition<R, D>::indexOfTriple(
            size_t const j,
            size_t const k,
            size_t const l
    ) const {
        assert(j > k && k > l);
        return j * (j - 1) * (j - 2) / 6 + k * (k - 1) / 2 + l;
    }

    template<class R, class D>
    inline
    size_t
    ProblemPartition<R, D>::indexOfEdge(
            size_t const j,
            size_t const k
    ) const {
        assert(j > k);
        return j * (j - 1) / 2 + k;
    }

    template<class R, class D>
    inline
    void
    ProblemPartition<R, D>::initCosts() {
        // at this point, the labels should be set

        for (size_t i = 0; i < numberOfElements_; ++i){
            for (size_t j = 0; j < i; ++j){


                if (labels_[i] == labels_[j]){
                    // joined edge, draw from intraDistribution
                    costsOfEdges_[indexOfEdge(i, j)] = intraDistribution_(engine_);
                } else {
                    // cut edge, draw from interDistribution
                    costsOfEdges_[indexOfEdge(i, j)] = interDistribution_(engine_);
                }


                for (size_t k = 0; k < j; ++k){
                    if (labels_[i] == labels_[j] &&
                        labels_[j] == labels_[k] &&
                        labels_[i] == labels_[k]
                    ){
                        // same cluster
                        costsOfTriples_[indexOfTriple(i, j, k)] = intraDistribution_(engine_);
                    } else {
                        costsOfTriples_[indexOfTriple(i, j, k)] = interDistribution_(engine_);
                    }
                }
            }
        }
    }

}

#endif
