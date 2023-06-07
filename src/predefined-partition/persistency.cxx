#include <andres/graph/multicut-cubic/persistency.hxx>
#include <random>
#include <iostream>
#include "dssc/partition/problem.hxx"
#include "dssc/csv-writer.hxx"
#include "dssc/cost-adjustment-adaptors/problem-cost-multiplication-adaptor.hxx"
#include "chrono"
#include "filesystem"


typedef dssc::ProblemPartition<double, std::normal_distribution<double>> Problem;
typedef andres::graph::multicut_cubic::Persistency<double> Persistency;
typedef dssc::CSVWriter CSVWriter;
typedef dssc::CostMultiplicationAdaptor<Problem, double> Scaler;


template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 3)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}


inline
CSVWriter
initializeCSVWriter(
        double const & alpha,
        double const & sigma0,
        double const & sigma1,
        size_t const & n,
        double const & edgeScale = 1.0,
        double const & tripletScale = 1.0,
        size_t const & seed = 4242,
        std::string const & fileName = "test.csv"
){
    CSVWriter csvWriter(fileName, {
            "seed",
            "sigma0",
            "sigma1",
            "alpha",
            "beta",
            "edgeScale",
            "tripletScale",
            "n",
            "orderIndex",
            "criterionName",
            "oneEdges",
            "zeroEdges",
            "zeroTriples",
            "remainingNodes",
            "remainingVariables",
            "remainingTriples",
            "runTime"
    }, {
                                std::to_string(seed),
                                to_string_with_precision(sigma0),
                                to_string_with_precision(sigma1),
                                to_string_with_precision(alpha),
                                to_string_with_precision(tripletScale),
                                to_string_with_precision(edgeScale),
                                to_string_with_precision(tripletScale),
                                std::to_string(n)
                        });

    return csvWriter;
}


inline
std::string
currentDateTime() {
    std::time_t t = std::time(nullptr);
    std::tm* now = std::localtime(&t);

    char buffer[128];
    strftime(buffer, sizeof(buffer), "%m-%d-%Y-%H-%M-%S", now);
    return buffer;
}

inline
void
doInitialSnapshot(
        Persistency & persistency,
        CSVWriter & csvWriter,
        int const & index = 0
        ){
    csvWriter.writeRow(
            {
                    std::to_string(index),
                    "initial",
                    std::to_string(persistency.edges1().size()),
                    std::to_string(persistency.edges0().size()),
                    std::to_string(persistency.triples0().size()),
                    std::to_string(persistency.remainingNodes()),
                    std::to_string(persistency.remainingVariables()),
                    std::to_string(persistency.remainingTriples())
            });
}

inline
void
doSubsetJoinCriterion(
        Persistency & persistency,
        CSVWriter & csvWriter,
        int const & index = 1
        ){

    persistency.subsetJoinCriterion();

    // orderIndex, criterionName, remainingNodes, remainingVariables
    csvWriter.writeRow({
                               std::to_string(index),
                               "subsetJoin",
                               std::to_string(persistency.edges1().size()),
                               std::to_string(persistency.edges0().size()),
                               std::to_string(persistency.triples0().size()),
                               std::to_string(persistency.remainingNodes()),
                               std::to_string(persistency.remainingVariables()),
                                std::to_string(persistency.remainingTriples()),
                               std::to_string(persistency.totalDurationSubsetJoinCriterion().count() * 1000)
                       });
}

inline
void
doEdgeSubgraphCriterion(
        Persistency & persistency,
        CSVWriter & csvWriter,
        int const & index = 2
        ){
    // join criteria
    persistency.edgeSubgraphCriterion();

    csvWriter.writeRow({
                               std::to_string(index),
                               "edgeSubgraph",
                               std::to_string(persistency.edges1().size()),
                               std::to_string(persistency.edges0().size()),
                               std::to_string(persistency.triples0().size()),
                               std::to_string(persistency.remainingNodes()),
                               std::to_string(persistency.remainingVariables()),
                               std::to_string(persistency.remainingTriples()),
                               std::to_string(persistency.totalDurationEdgeSubgraphCriterion().count() * 1000)
                       });
}

inline
void
doTripletSubgraphCriterion(
        Persistency & persistency,
        CSVWriter & csvWriter,
        int const & index = 3
        ){
    persistency.tripletSubgraphCriterion();

    csvWriter.writeRow({
                               std::to_string(index),
                               "tripletSubgraph",
                               std::to_string(persistency.edges1().size()),
                               std::to_string(persistency.edges0().size()),
                               std::to_string(persistency.triples0().size()),
                               std::to_string(persistency.remainingNodes()),
                               std::to_string(persistency.remainingVariables()),
                               std::to_string(persistency.remainingTriples()),
                               std::to_string(persistency.totalDurationTripletSubgraphCriterion().count() * 1000)
                       });
}

inline
void
doTripletEdgeJoinCriterion(
        Persistency & persistency,
        CSVWriter & csvWriter,
        int const & index = 4
){
    persistency.tripletEdgeJoinCriterion();

    csvWriter.writeRow({
                               std::to_string(index),
                               "tripletEdgeJoin",
                               std::to_string(persistency.edges1().size()),
                               std::to_string(persistency.edges0().size()),
                               std::to_string(persistency.triples0().size()),
                               std::to_string(persistency.remainingNodes()),
                               std::to_string(persistency.remainingVariables()),
                               std::to_string(persistency.remainingTriples()),
                               std::to_string(persistency.totalDurationTripletEdgeJoinCriterion().count() * 1000)
                       });
}


inline
void
doEdgeJoinCriterion(
        Persistency & persistency,
        CSVWriter & csvWriter,
        int const & index = 5
        ){
    persistency.edgeJoinCriterion();

    csvWriter.writeRow({
                               std::to_string(index),
                               "edgeJoin",
                               std::to_string(persistency.edges1().size()),
                               std::to_string(persistency.edges0().size()),
                               std::to_string(persistency.triples0().size()),
                               std::to_string(persistency.remainingNodes()),
                               std::to_string(persistency.remainingVariables()),
                               std::to_string(persistency.remainingTriples()),
                               std::to_string(persistency.totalDurationEdgeJoinCriterion().count() * 1000)
                       });
}

inline
void
doTripletJoinCriterion(
        Persistency & persistency,
        CSVWriter & csvWriter,
        int const & index = 6
        ){
    persistency.tripletJoinCriterion();

    csvWriter.writeRow({
                               std::to_string(index),
                               "tripletJoin",
                               std::to_string(persistency.edges1().size()),
                               std::to_string(persistency.edges0().size()),
                               std::to_string(persistency.triples0().size()),
                               std::to_string(persistency.remainingNodes()),
                               std::to_string(persistency.remainingVariables()),
                               std::to_string(persistency.remainingTriples()),
                               std::to_string(persistency.totalDurationTripletJoinCriterion().count() * 1000)
                       });
}

inline
void
doEdgeCutCriterion(
        Persistency & persistency,
        CSVWriter & csvWriter,
        int const & index = 7
        ){
    // cut criteria
    persistency.edgeCutCriterion();

    csvWriter.writeRow({
                               std::to_string(index),
                               "edgeCut",
                               std::to_string(persistency.edges1().size()),
                               std::to_string(persistency.edges0().size()),
                               std::to_string(persistency.triples0().size()),
                               std::to_string(persistency.remainingNodes()),
                               std::to_string(persistency.remainingVariables()),
                               std::to_string(persistency.remainingTriples()),
                               std::to_string(persistency.totalDurationEdgeCutCriterion().count() * 1000)
                       });
}

inline
void
doTripletCutCriterion(
        Persistency & persistency,
        CSVWriter & csvWriter,
        int const & index = 8
        ){
    persistency.tripletCutCriterion();

    csvWriter.writeRow({
                               std::to_string(index),
                               "tripletCut",
                               std::to_string(persistency.edges1().size()),
                               std::to_string(persistency.edges0().size()),
                               std::to_string(persistency.triples0().size()),
                               std::to_string(persistency.remainingNodes()),
                               std::to_string(persistency.remainingVariables()),
                               std::to_string(persistency.remainingTriples()),
                               std::to_string(persistency.totalDurationTripletCutCriterion().count() * 1000)
                       });
}

inline
void
doFindIndependentSubsets(
        Persistency & persistency,
        CSVWriter & csvWriter,
        int const & index = 9
){
    persistency.findIndependentSubProblems();

    csvWriter.writeRow({
                               std::to_string(index),
                               "findIndependentSubProblems",
                               std::to_string(persistency.edges1().size()),
                               std::to_string(persistency.edges0().size()),
                               std::to_string(persistency.triples0().size()),
                               std::to_string(persistency.remainingNodes()),
                               std::to_string(persistency.remainingVariables()),
                               std::to_string(persistency.remainingTriples()),
                               std::to_string(persistency.totalDurationFindIndependentProblems().count() * 1000)
                       });
}


inline
void
runPredefinedPartitionInstanceAllCriteria(
        double const & alpha,
        double const & sigma0,
        double const & sigma1,
        size_t const & n,
        double const & edgeScale = 1.0,
        double const & tripletScale = 1.0,
        size_t const & seed = 4242,
        std::string const & csvFileName = "test.csv"
){
    assert(alpha >= 0 && alpha <= 1);
    std::vector<size_t > const & clusterSizes = {n, 2*n, 2*n, 3*n};

    double intraMean = -1 + alpha;
    double interMean = 1 - alpha;
    double sigma = sigma0 + alpha*(sigma1 - sigma0);

    std::normal_distribution<double> interDistribution(interMean, sigma);
    std::normal_distribution<double> intraDistribution(intraMean, sigma);

    Problem problem = Problem (
            clusterSizes,
            intraDistribution,
            interDistribution,
            seed
    );

    Scaler scaler(problem, edgeScale, tripletScale);

    Persistency persistency(scaler);

    CSVWriter csvWriter = initializeCSVWriter(
            alpha,
            sigma0,
            sigma1,
            n,
            edgeScale,
            tripletScale,
            seed,
            csvFileName
    );

    doInitialSnapshot(persistency, csvWriter);
    doSubsetJoinCriterion(persistency, csvWriter);
    doEdgeSubgraphCriterion(persistency, csvWriter);
    doTripletSubgraphCriterion(persistency, csvWriter);
    doTripletEdgeJoinCriterion(persistency, csvWriter);
    doEdgeJoinCriterion(persistency, csvWriter);
    doTripletJoinCriterion(persistency, csvWriter);
    doEdgeCutCriterion(persistency, csvWriter);
    doTripletCutCriterion(persistency, csvWriter);

    csvWriter.close();
}


inline
void
runPredefinedPartitionSubsetJoin(
        double const & alpha,
        double const & sigma0,
        double const & sigma1,
        size_t const & n,
        double const & edgeScale = 1.0,
        double const & tripletScale = 1.0,
        size_t const & seed = 4242,
        std::string const & csvBaseFileName = "test",
        std::filesystem::path const & directory = "./test"
        ){
    if (!std::filesystem::exists(directory)){
        std::filesystem::create_directories(directory);
    }

    assert(alpha >= 0 && alpha <= 1);
    std::vector<size_t > const & clusterSizes = {n, 2*n, 2*n, 3*n};

    double intraMean = -1 + alpha;
    double interMean = 1 - alpha;
    double sigma = sigma0 + alpha*(sigma1 - sigma0);

    std::normal_distribution<double> interDistribution(interMean, sigma);
    std::normal_distribution<double> intraDistribution(intraMean, sigma);

    Problem problem = Problem (
            clusterSizes,
            intraDistribution,
            interDistribution,
            seed
    );

    Scaler scaler(problem, edgeScale, tripletScale);

    // subset join - criterion 1
    Persistency persistency(scaler);

    std::filesystem::path subsetJoinFilePath = directory / std::filesystem::path(csvBaseFileName + "_subsetJoin.csv");

    CSVWriter csvWriter = initializeCSVWriter(
            alpha,
            sigma0,
            sigma1,
            n,
            edgeScale,
            tripletScale,
            seed,
            subsetJoinFilePath
    );

    doInitialSnapshot(
            persistency,
            csvWriter,
            0
    );

    doSubsetJoinCriterion(
            persistency,
            csvWriter,
            1
    );
    csvWriter.close();
}

inline
void
runPredefinedPartitionEdgeSubgraph(
        double const & alpha,
        double const & sigma0,
        double const & sigma1,
        size_t const & n,
        double const & edgeScale = 1.0,
        double const & tripletScale = 1.0,
        size_t const & seed = 4242,
        std::string const & csvBaseFileName = "test",
        std::filesystem::path const & directory = "./test"
){
    if (!std::filesystem::exists(directory)){
        std::filesystem::create_directories(directory);
    }

    assert(alpha >= 0 && alpha <= 1);
    std::vector<size_t > const & clusterSizes = {n, 2*n, 2*n, 3*n};

    double intraMean = -1 + alpha;
    double interMean = 1 - alpha;
    double sigma = sigma0 + alpha*(sigma1 - sigma0);

    std::normal_distribution<double> interDistribution(interMean, sigma);
    std::normal_distribution<double> intraDistribution(intraMean, sigma);

    Problem problem = Problem (
            clusterSizes,
            intraDistribution,
            interDistribution,
            seed
    );

    Scaler scaler(problem, edgeScale, tripletScale);

    Persistency persistency(scaler);

    // edge subgraph join - criterion 1
    std::filesystem::path edgeSubgraphCriterionFilePath = directory / std::filesystem::path(csvBaseFileName + "_edgeSubgraph.csv");

    CSVWriter csvWriter = initializeCSVWriter(
            alpha,
            sigma0,
            sigma1,
            n,
            edgeScale,
            tripletScale,
            seed,
            edgeSubgraphCriterionFilePath
    );

    doInitialSnapshot(
            persistency,
            csvWriter,
            0
    );

    doEdgeSubgraphCriterion(
            persistency,
            csvWriter,
            1
    );
    csvWriter.close();
}

inline
void
runPredefinedPartitionTripletSubgraph(
        double const & alpha,
        double const & sigma0,
        double const & sigma1,
        size_t const & n,
        double const & edgeScale = 1.0,
        double const & tripletScale = 1.0,
        size_t const & seed = 4242,
        std::string const & csvBaseFileName = "test",
        std::filesystem::path const & directory = "./test"
){
    if (!std::filesystem::exists(directory)){
        std::filesystem::create_directories(directory);
    }

    assert(alpha >= 0 && alpha <= 1);
    std::vector<size_t > const & clusterSizes = {n, 2*n, 2*n, 3*n};

    double intraMean = -1 + alpha;
    double interMean = 1 - alpha;
    double sigma = sigma0 + alpha*(sigma1 - sigma0);

    std::normal_distribution<double> interDistribution(interMean, sigma);
    std::normal_distribution<double> intraDistribution(intraMean, sigma);

    Problem problem = Problem (
            clusterSizes,
            intraDistribution,
            interDistribution,
            seed
    );

    Scaler scaler(problem, edgeScale, tripletScale);

    Persistency persistency(scaler);

    // edge subgraph join - criterion 1
    std::filesystem::path tripletSubgraphFileName = directory / std::filesystem::path(csvBaseFileName + "_tripletSubgraph.csv");

    CSVWriter csvWriter = initializeCSVWriter(
            alpha,
            sigma0,
            sigma1,
            n,
            edgeScale,
            tripletScale,
            seed,
            tripletSubgraphFileName
    );

    doInitialSnapshot(
            persistency,
            csvWriter,
            0
    );

    doTripletSubgraphCriterion(
            persistency,
            csvWriter,
            1
    );
    csvWriter.close();
}


inline
void
runPredefinedPartitionTripletEdgeJoin(
        double const & alpha,
        double const & sigma0,
        double const & sigma1,
        size_t const & n,
        double const & edgeScale = 1.0,
        double const & tripletScale = 1.0,
        size_t const & seed = 4242,
        std::string const & csvBaseFileName = "test",
        std::filesystem::path const & directory = "./test"
){
    if (!std::filesystem::exists(directory)){
        std::filesystem::create_directories(directory);
    }

    assert(alpha >= 0 && alpha <= 1);
    std::vector<size_t > const & clusterSizes = {n, 2*n, 2*n, 3*n};

    double intraMean = -1 + alpha;
    double interMean = 1 - alpha;
    double sigma = sigma0 + alpha*(sigma1 - sigma0);

    std::normal_distribution<double> interDistribution(interMean, sigma);
    std::normal_distribution<double> intraDistribution(intraMean, sigma);

    Problem problem = Problem (
            clusterSizes,
            intraDistribution,
            interDistribution,
            seed
    );

    Scaler scaler(problem, edgeScale, tripletScale);

    Persistency persistency(scaler);

    // edge subgraph join - criterion 1
    std::filesystem::path tripletEdgeJoinFileName = directory / std::filesystem::path(csvBaseFileName + "_tripletEdgeJoin.csv");

    CSVWriter csvWriter = initializeCSVWriter(
            alpha,
            sigma0,
            sigma1,
            n,
            edgeScale,
            tripletScale,
            seed,
            tripletEdgeJoinFileName
    );

    doInitialSnapshot(
            persistency,
            csvWriter,
            0
    );

    doTripletEdgeJoinCriterion(
            persistency,
            csvWriter,
            1
    );
    csvWriter.close();
}

inline
void
runPredefinedPartitionEdgeJoin(
        double const & alpha,
        double const & sigma0,
        double const & sigma1,
        size_t const & n,
        double const & edgeScale = 1.0,
        double const & tripletScale = 1.0,
        size_t const & seed = 4242,
        std::string const & csvBaseFileName = "test",
        std::filesystem::path const & directory = "./test"
){
    if (!std::filesystem::exists(directory)){
        std::filesystem::create_directories(directory);
    }

    assert(alpha >= 0 && alpha <= 1);
    std::vector<size_t > const & clusterSizes = {n, 2*n, 2*n, 3*n};

    double intraMean = -1 + alpha;
    double interMean = 1 - alpha;
    double sigma = sigma0 + alpha*(sigma1 - sigma0);

    std::normal_distribution<double> interDistribution(interMean, sigma);
    std::normal_distribution<double> intraDistribution(intraMean, sigma);

    Problem problem = Problem (
            clusterSizes,
            intraDistribution,
            interDistribution,
            seed
    );

    Scaler scaler(problem, edgeScale, tripletScale);

    Persistency persistency(scaler);

    // edge subgraph join - criterion 1
    std::filesystem::path edgeJoinFileName = directory / std::filesystem::path(csvBaseFileName + "_edgeJoin.csv");

    CSVWriter csvWriter = initializeCSVWriter(
            alpha,
            sigma0,
            sigma1,
            n,
            edgeScale,
            tripletScale,
            seed,
            edgeJoinFileName
    );

    doInitialSnapshot(
            persistency,
            csvWriter,
            0
    );

    doEdgeJoinCriterion(
            persistency,
            csvWriter,
            1
    );
    csvWriter.close();
}

inline
void
runPredefinedPartitionTripletJoin(
        double const & alpha,
        double const & sigma0,
        double const & sigma1,
        size_t const & n,
        double const & edgeScale = 1.0,
        double const & tripletScale = 1.0,
        size_t const & seed = 4242,
        std::string const & csvBaseFileName = "test",
        std::filesystem::path const & directory = "./test"
){
    if (!std::filesystem::exists(directory)){
        std::filesystem::create_directories(directory);
    }

    assert(alpha >= 0 && alpha <= 1);
    std::vector<size_t > const & clusterSizes = {n, 2*n, 2*n, 3*n};

    double intraMean = -1 + alpha;
    double interMean = 1 - alpha;
    double sigma = sigma0 + alpha*(sigma1 - sigma0);

    std::normal_distribution<double> interDistribution(interMean, sigma);
    std::normal_distribution<double> intraDistribution(intraMean, sigma);

    Problem problem = Problem (
            clusterSizes,
            intraDistribution,
            interDistribution,
            seed
    );

    Scaler scaler(problem, edgeScale, tripletScale);

    Persistency persistency(scaler);

    // edge subgraph join - criterion 1
    std::filesystem::path tripletJoinFileName = directory / std::filesystem::path(csvBaseFileName + "_tripletJoin.csv");

    CSVWriter csvWriter = initializeCSVWriter(
            alpha,
            sigma0,
            sigma1,
            n,
            edgeScale,
            tripletScale,
            seed,
            tripletJoinFileName
    );

    doInitialSnapshot(
            persistency,
            csvWriter,
            0
    );

    doTripletJoinCriterion(
            persistency,
            csvWriter,
            1
    );
    csvWriter.close();
}



inline
void
runPredefinedPartitionEdgeCut(
        double const & alpha,
        double const & sigma0,
        double const & sigma1,
        size_t const & n,
        double const & edgeScale = 1.0,
        double const & tripletScale = 1.0,
        size_t const & seed = 4242,
        std::string const & csvBaseFileName = "test",
        std::filesystem::path const & directory = "./test"
){
    if (!std::filesystem::exists(directory)){
        std::filesystem::create_directories(directory);
    }

    assert(alpha >= 0 && alpha <= 1);
    std::vector<size_t > const & clusterSizes = {n, 2*n, 2*n, 3*n};

    double intraMean = -1 + alpha;
    double interMean = 1 - alpha;
    double sigma = sigma0 + alpha*(sigma1 - sigma0);

    std::normal_distribution<double> interDistribution(interMean, sigma);
    std::normal_distribution<double> intraDistribution(intraMean, sigma);

    Problem problem = Problem (
            clusterSizes,
            intraDistribution,
            interDistribution,
            seed
    );

    Scaler scaler(problem, edgeScale, tripletScale);

    Persistency persistency(scaler);

    // edge subgraph join - criterion 1
    std::filesystem::path edgeCutFileName = directory / std::filesystem::path(csvBaseFileName + "_edgeCut.csv");

    CSVWriter csvWriter = initializeCSVWriter(
            alpha,
            sigma0,
            sigma1,
            n,
            edgeScale,
            tripletScale,
            seed,
            edgeCutFileName
    );

    doInitialSnapshot(
            persistency,
            csvWriter,
            0
    );

    doEdgeCutCriterion(
            persistency,
            csvWriter,
            1
    );
    csvWriter.close();
}

inline
void
runPredefinedPartitionTripletCut(
        double const & alpha,
        double const & sigma0,
        double const & sigma1,
        size_t const & n,
        double const & edgeScale = 1.0,
        double const & tripletScale = 1.0,
        size_t const & seed = 4242,
        std::string const & csvBaseFileName = "test",
        std::filesystem::path const & directory = "./test"
){
    if (!std::filesystem::exists(directory)){
        std::filesystem::create_directories(directory);
    }

    assert(alpha >= 0 && alpha <= 1);
    std::vector<size_t > const & clusterSizes = {n, 2*n, 2*n, 3*n};

    double intraMean = -1 + alpha;
    double interMean = 1 - alpha;
    double sigma = sigma0 + alpha*(sigma1 - sigma0);

    std::normal_distribution<double> interDistribution(interMean, sigma);
    std::normal_distribution<double> intraDistribution(intraMean, sigma);

    Problem problem = Problem (
            clusterSizes,
            intraDistribution,
            interDistribution,
            seed
    );

    Scaler scaler(problem, edgeScale, tripletScale);

    Persistency persistency(scaler);

    // edge subgraph join - criterion 1
    std::filesystem::path tripletCutCriterion = directory / std::filesystem::path(csvBaseFileName + "_tripletCut.csv");

    CSVWriter csvWriter = initializeCSVWriter(
            alpha,
            sigma0,
            sigma1,
            n,
            edgeScale,
            tripletScale,
            seed,
            tripletCutCriterion
    );

    doInitialSnapshot(
            persistency,
            csvWriter,
            0
    );

    doTripletCutCriterion(
            persistency,
            csvWriter,
            1
    );
    csvWriter.close();
}



inline
void
runPredefinedPartitionFindIndependentSubsets(
        double const & alpha,
        double const & sigma0,
        double const & sigma1,
        size_t const & n,
        double const & edgeScale = 1.0,
        double const & tripletScale = 1.0,
        size_t const & seed = 4242,
        std::string const & csvBaseFileName = "test",
        std::filesystem::path const & directory = "./test"
){
    if (!std::filesystem::exists(directory)){
        std::filesystem::create_directories(directory);
    }

    assert(alpha >= 0 && alpha <= 1);
    std::vector<size_t > const & clusterSizes = {n, 2*n, 2*n, 3*n};

    double intraMean = -1 + alpha;
    double interMean = 1 - alpha;
    double sigma = sigma0 + alpha*(sigma1 - sigma0);

    std::normal_distribution<double> interDistribution(interMean, sigma);
    std::normal_distribution<double> intraDistribution(intraMean, sigma);

    Problem problem = Problem (
            clusterSizes,
            intraDistribution,
            interDistribution,
            seed
    );

    Scaler scaler(problem, edgeScale, tripletScale);

    Persistency persistency(scaler);

    // edge subgraph join - criterion 1
    std::filesystem::path findIndependentSubsetsPath = directory / std::filesystem::path(csvBaseFileName + "_findIndependentSubsets.csv");

    CSVWriter csvWriter = initializeCSVWriter(
            alpha,
            sigma0,
            sigma1,
            n,
            edgeScale,
            tripletScale,
            seed,
            findIndependentSubsetsPath
    );

    doInitialSnapshot(
            persistency,
            csvWriter,
            0
    );

    doFindIndependentSubsets(
            persistency,
            csvWriter,
            1
    );
    csvWriter.close();
}

// runs all criteria in succession with a join and a cut phase
inline
void
runPredefinedPartitionInstanceAllCriteria(
        double const & alpha,
        double const & sigma0,
        double const & sigma1,
        size_t const & n,
        double const & edgeScale = 1.0,
        double const & tripletScale = 1.0,
        size_t const & seed = 4242,
        std::string const & csvBaseFileName = "test",
        std::filesystem::path const & directory = "./test"
){
    if (!std::filesystem::exists(directory)) {
        std::filesystem::create_directories(directory);
    }

    assert(alpha >= 0 && alpha <= 1);
    std::vector<size_t > const & clusterSizes = {n, 2*n, 2*n, 3*n};

    double intraMean = -1 + alpha;
    double interMean = 1 - alpha;
    double sigma = sigma0 + alpha*(sigma1 - sigma0);

    std::normal_distribution<double> interDistribution(interMean, sigma);
    std::normal_distribution<double> intraDistribution(intraMean, sigma);

    Problem problem = Problem (
            clusterSizes,
            intraDistribution,
            interDistribution,
            seed
    );

    Scaler scaler(problem, edgeScale, tripletScale);

    Persistency persistency(scaler);


    std::filesystem::path fileName = directory / std::filesystem::path(csvBaseFileName + "_all.csv");

    CSVWriter csvWriter = initializeCSVWriter(
            alpha,
            sigma0,
            sigma1,
            n,
            edgeScale,
            tripletScale,
            seed,
            fileName
    );

    doInitialSnapshot(persistency, csvWriter, 0);

    // join phase
    bool didJoin = false;
    size_t oldNumberOfNodes = persistency.remainingNodes();
    int index = 1;

    do {
        didJoin = false;

        // cut phase
        doFindIndependentSubsets(persistency, csvWriter, index);
        index++;

        doSubsetJoinCriterion(persistency, csvWriter, index);
        // ensures, that the loop is only picked up again if a subsequent criterion changes the number of nodes
        oldNumberOfNodes = persistency.remainingNodes();
        index++;
        doEdgeJoinCriterion(persistency, csvWriter, index);
        index++;
        doTripletEdgeJoinCriterion(persistency, csvWriter, index);
        index++;
        doEdgeSubgraphCriterion(persistency, csvWriter, index);
        index++;
        doTripletSubgraphCriterion(persistency, csvWriter, index);
        index++;
        doTripletJoinCriterion(persistency, csvWriter, index);
        index++;

        if (persistency.remainingNodes() != oldNumberOfNodes){
            oldNumberOfNodes = persistency.remainingNodes();
            didJoin = true;
        }
    } while (didJoin);

    doEdgeCutCriterion(persistency, csvWriter, index);
    index++;
    doTripletCutCriterion(persistency, csvWriter, index);
    index++;
}



inline
void
runAllCriteriaExperiments(
        size_t const n
        ){
    double const & sigma0 = 0.1;
    double const & sigma1 = 0.4;

    std::vector<double> alphas = {
            0.0, 0.025, 0.05, 0.075,
            0.1, 0.125, 0.15, 0.175,
            0.2, 0.225, 0.25, 0.275,
            0.3, 0.325, 0.35, 0.375,
            0.4, 0.425, 0.45, 0.475,
            0.5, 0.525, 0.55, 0.575,
            0.6, 0.625, 0.65, 0.675,
            0.7, 0.725, 0.75, 0.775,
            0.8, 0.825, 0.85, 0.875, 0.9, 1.0
    };
    std::vector<double> betas = {0.0, 0.5, 1.0};

    size_t const & startSeed = 4242;
    size_t const & numSeeds = 30;
    std::vector<size_t > seeds(numSeeds);
    for (size_t seedIndex = 0; seedIndex < numSeeds; ++seedIndex){
        seeds[seedIndex] = startSeed + seedIndex;
    }

    auto totalStart = std::chrono::system_clock::now();

    std::string baseFileName = "experiment";

    std::string outDir = "partition_" + to_string_with_precision(8*n, 0) + "nodes_" + to_string_with_precision(numSeeds, 0) + "seeds";

    std::cout << "Running with n=" << n << std::endl;

    for (auto & alpha: alphas){
        for (auto & beta : betas){

            std::cout << "alpha: " << alpha << ", beta: " << beta << std::endl;
            auto innerStart = std::chrono::system_clock::now();

            for (auto & seed: seeds){
                std::string fileName = baseFileName;
                runPredefinedPartitionInstanceAllCriteria(alpha, sigma0, sigma1, n, 1 - beta, beta, seed, fileName, outDir);
        }

        auto innerEnd = std::chrono::system_clock::now();
        std::chrono::duration<double> innerDuration = innerEnd - innerStart;
        std::cout << "Time taken: " << innerDuration.count() << "s" << std::endl;
    }
    }

    auto totalEnd = std::chrono::system_clock::now();

    std::chrono::duration<double> totalDuration = totalEnd - totalStart;

    std::cout << "Total Runtime: " << totalDuration.count() / 60 << "min" << std::endl;
}


// runs only one criterion at a time
inline
void
runPredefinedPartitionInstanceIndividualCriteria(
        double const & alpha,
        double const & sigma0,
        double const & sigma1,
        size_t const & n,
        double const & edgeScale = 1.0,
        double const & tripletScale = 1.0,
        size_t const & seed = 4242,
        std::string const & csvBaseFileName = "test",
        std::filesystem::path const & directory = "./test"
){
    runPredefinedPartitionSubsetJoin(alpha, sigma0, sigma1, n, edgeScale, tripletScale, seed, csvBaseFileName, directory);
    runPredefinedPartitionEdgeSubgraph(alpha, sigma0, sigma1, n, edgeScale, tripletScale, seed, csvBaseFileName, directory);
    runPredefinedPartitionTripletEdgeJoin(alpha, sigma0, sigma1, n, edgeScale, tripletScale, seed, csvBaseFileName, directory);
    runPredefinedPartitionTripletSubgraph(alpha, sigma0, sigma1, n, edgeScale, tripletScale, seed, csvBaseFileName, directory);
    runPredefinedPartitionEdgeJoin(alpha, sigma0, sigma1, n, edgeScale, tripletScale, seed, csvBaseFileName, directory);
    runPredefinedPartitionTripletJoin(alpha, sigma0, sigma1, n, edgeScale, tripletScale, seed, csvBaseFileName, directory);
    runPredefinedPartitionEdgeCut(alpha, sigma0, sigma1, n, edgeScale, tripletScale, seed, csvBaseFileName, directory);
    runPredefinedPartitionTripletCut(alpha, sigma0, sigma1, n, edgeScale, tripletScale, seed, csvBaseFileName, directory);
    runPredefinedPartitionFindIndependentSubsets(alpha, sigma0, sigma1, n, edgeScale, tripletScale, seed, csvBaseFileName, directory);
}

inline
void
runIndividualExperiments(size_t const n){
    double const & sigma0 = 0.1;
    double const & sigma1 = 0.4;

    // alpha from 0 to 1
    std::vector<double> alphas = {
        0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
    };

    // triplet costs are scaled with beta
    // edge costs are scaled with 1 - beta
    std::vector<double> betas = {0.0, 0.5, 1.0};

    size_t const & startSeed = 4242;
    size_t const & numSeeds = 30;
    std::vector<size_t > seeds(numSeeds);
    for (size_t seedIndex = 0; seedIndex < numSeeds; ++seedIndex){
        seeds[seedIndex] = startSeed + seedIndex;
    }

    auto totalStart = std::chrono::system_clock::now();

    std::string baseFileName = "experiment_full_";

        std::string outDir = "partition_individual_" + to_string_with_precision(8*n, 0) + "nodes_" + to_string_with_precision(numSeeds, 0) + "seeds";

            for (auto & alpha: alphas){
                for (auto & beta : betas){

                    auto innerStart = std::chrono::system_clock::now();

                    std::cout << "alpha: " << alpha << ", beta: " << beta << std::endl;

                    for (auto & seed: seeds){


                        std::string fileName = baseFileName;
                        fileName += "n=" + std::to_string(n);

                        runPredefinedPartitionInstanceIndividualCriteria(alpha, sigma0, sigma1, n, 1 - beta, beta, seed, fileName, outDir);
                }

                    auto innerEnd = std::chrono::system_clock::now();
                    std::chrono::duration<double> innerDuration = innerEnd - innerStart;

                    std::cout << "Time taken: " << innerDuration.count() << "s" << std::endl;
                }
    }

    auto totalEnd = std::chrono::system_clock::now();

    std::chrono::duration<double> totalDuration = totalEnd - totalStart;

    std::cout << "Total Runtime: " << totalDuration.count() / 60 << "min" << std::endl;
}

int main(int argc, char** argv) {
    size_t defaultAllN = 6;
    size_t defaultIndividualN = 6;

    if (argc == 1){
        throw std::runtime_error("Use '-i' or '-a' to run the criteria individually or jointly. Supply a number n as first parameter. This will run with 8*n nodes for the partition instance. Default is n=6.");
    } else if (argc == 2){
        std::string flag = argv[1];
        if (flag == "-i" || flag == "--individual"){
            std::cout << "Running with n = " << defaultIndividualN << std::endl;
            std:: cout << "Running criteria separately. Use --individual or -i flag to run each criterion on its own." << std::endl;
            runIndividualExperiments(defaultIndividualN);
        } else if (flag == "-a" || flag == "--all") {
            std::cout << "Running with n = " << defaultAllN << std::endl;
            runAllCriteriaExperiments(defaultAllN);
        } else {
            throw std::runtime_error("Flag " + flag + " not allowed.");
        }
    } else if (argc == 3){
        std::string flag = argv[1];
        std::istringstream ss(argv[2]);
        size_t maxN;
        if (!(ss >> maxN)){
            throw std::runtime_error("Invalid number");
        } else if (!ss.eof()){
            throw std::runtime_error("Trailing characters after number.");
        } else {
            std::cout << "Running with n = " << maxN << std::endl;
        }

        if (flag == "-i" || flag == "--individual"){
            std:: cout << "Running criteria separately. Use --individual or -i flag to run each criterion on its own." << std::endl;
            runIndividualExperiments(maxN);
        } else if (flag == "-a" || flag == "--all") {
            runAllCriteriaExperiments(maxN);
        } else {
            throw std::runtime_error("Flag " + flag + " not allowed.");
        }
    }
    else {
        throw std::runtime_error("Cannot supply more than one argument.");
    }

    return 0;
}