// Applies all criteria jointly to partition instances of growing size.
// The code up to runPredefinedPartitionInstanceAllCriteria is copied from persistency.cxx.

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
runInstanceSizeExperiments(
        size_t const numSeeds,
        size_t const maxN
        ){
    double const & sigma0 = 0.1;
    double const & sigma1 = 0.4;

    std::vector<double> alphas = {0.0, 0.4, 0.5, 0.65};
    // triplet costs are scaled with beta
    // edge costs are scaled with 1 - beta
    std::vector<double> betas = {0.5};

    size_t const & startSeed = 4242;

    std::vector<size_t > seeds(numSeeds);
    for (size_t seedIndex = 0; seedIndex < numSeeds; ++seedIndex){
        seeds[seedIndex] = startSeed + seedIndex;
    }

    auto totalStart = std::chrono::system_clock::now();

    std::string outDir = "partition_instance_sizes_" + to_string_with_precision(numSeeds, 0) + "seeds";

    for (size_t n = 1; n <= maxN; ++n){
        std::cout << "Running with numNodes: " << 8*n << std::endl;

        std::string fileName = to_string_with_precision(8*n, 0) + "nodes";

        for (auto & alpha: alphas){
            for (auto & beta : betas){
                std::cout << "alpha: " << alpha << ", beta: " << beta << std::endl;
                auto innerStart = std::chrono::system_clock::now();

                for (auto & seed: seeds){
                    runPredefinedPartitionInstanceAllCriteria(alpha, sigma0, sigma1, n, 1 - beta, beta, seed, fileName, outDir);
                }

                auto innerEnd = std::chrono::system_clock::now();
                std::chrono::duration<double> innerDuration = innerEnd - innerStart;
                std::cout << "Time taken: " << innerDuration.count() << "s" << std::endl;
            }
        }
    }

    auto totalEnd = std::chrono::system_clock::now();
    std::chrono::duration<double> totalDuration = totalEnd - totalStart;
    std::cout << "Total Runtime: " << totalDuration.count() / 60 << "min" << std::endl;
}

inline
size_t
parseNumber(std::string const & text){
    std::istringstream ss(text);
    size_t number;
    if (!(ss >> number)){
        throw std::runtime_error("Invalid number");
    } else if (!ss.eof()){
        throw std::runtime_error("Trailing characters after number.");
    }
    return number;
}


int main(int argc, char** argv) {
    size_t numSeeds = 30;
    size_t maxN = 13;

    for (int i = 1; i < argc; i += 2){
        std::string flag = argv[i];
        if (i + 1 >= argc){
            throw std::runtime_error("Flag " + flag + " needs a number.");
        }
        if (flag == "-s" || flag == "--seeds"){
            numSeeds = parseNumber(argv[i + 1]);
        } else if (flag == "-n" || flag == "--max-n"){
            maxN = parseNumber(argv[i + 1]);
        } else {
            throw std::runtime_error("Flag " + flag + " not allowed. Use --seeds <number> and --max-n <n>.");
        }
    }

    if (numSeeds == 0 || maxN < 1){
        throw std::runtime_error("Use at least one seed and n >= 1, i.e. instances with 8*n nodes.");
    }

    std::cout << "Running with " << numSeeds << " seeds and n = 1, ..., " << maxN
              << ", i.e. up to " << 8*maxN << " nodes." << std::endl;
    runInstanceSizeExperiments(numSeeds, maxN);

    return 0;
}
