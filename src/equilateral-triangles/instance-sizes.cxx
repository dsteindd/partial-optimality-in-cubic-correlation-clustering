// Applies all criteria jointly to geometric instances of growing size.
// The code up to runEquilateralTrianglesInstanceAllCriteria is copied from persistency.cxx.

#include <andres/graph/multicut-cubic/persistency.hxx>
#include <random>
#include <iostream>
#include "dssc/equilateral-triangles/problem.hxx"
#include "dssc/equilateral-triangles/problem-adaptor.hxx"
#include "dssc/csv-writer.hxx"
#include "chrono"
#include "filesystem"

typedef dssc::ProblemEquilateralTriangles<double> Problem;
typedef dssc::AngleCostMode AngleCostMode;
template<AngleCostMode AngleCost>
using Adaptor = dssc::ProblemEquilateralTrianglesAdaptor<
        double,
        AngleCost
>;


typedef andres::graph::multicut_cubic::Persistency<double> Persistency;
typedef dssc::CSVWriter CSVWriter;
typedef Problem::Point2 Point2;


template<typename T>
std::string to_string_with_precision(const T a_value, const int n = 3) {
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}


inline
CSVWriter
initializeCSVWriter(
        double const &sigma,
        size_t const &numberOfNodes,
        size_t const &numberOfExpectedRemainingNodes,
        double const &costAtDistanceZero,
        double const &distanceAtCostZero,
        double const &edgeMaxLengthThreshold,
        double const &tripletMaxEdgeLengthThreshold,
        size_t const &seed = 4242,
        std::string const &fileName = "test.csv"
) {
    CSVWriter csvWriter(fileName, {
            "seed",
            "sigma",
            "numberOfNodes",
            "numberOfExpectedRemainingNodes",
            "costAtDistanceZero",
            "distanceAtCostZero",
            "edgeMaxLengthThreshold",
            "tripletMaxEdgeLengthThreshold",
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
                                to_string_with_precision(sigma, 5),
                                std::to_string(numberOfNodes),
                                std::to_string(numberOfExpectedRemainingNodes),
                                to_string_with_precision(costAtDistanceZero, 5),
                                to_string_with_precision(distanceAtCostZero, 5),
                                to_string_with_precision(edgeMaxLengthThreshold, 5),
                                to_string_with_precision(tripletMaxEdgeLengthThreshold, 5)
                        });

    return csvWriter;
}


inline
void
doInitialSnapshot(
        Persistency &persistency,
        CSVWriter &csvWriter,
        int const &index = 0
) {
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
        Persistency &persistency,
        CSVWriter &csvWriter,
        int const &index = 1
) {

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
        Persistency &persistency,
        CSVWriter &csvWriter,
        int const &index = 2
) {
    // join criteria
//    std::cout << "Run: edgeSubgraphCriterion" << std::endl;
    persistency.edgeSubgraphCriterion();
//    std::cout << "Done: edgeSubgraphCriterion" << std::endl;

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
        Persistency &persistency,
        CSVWriter &csvWriter,
        int const &index = 3
) {
//    std::cout << "Run: tripletSubgraphCriterion" << std::endl;
    persistency.tripletSubgraphCriterion();
//    std::cout << "Done: tripletSubgraphCriterion" << std::endl;

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
        Persistency &persistency,
        CSVWriter &csvWriter,
        int const &index = 4
) {
//    std::cout << "Run: tripletEdgeJoinCriterion" << std::endl;
    persistency.tripletEdgeJoinCriterion();
//    std::cout << "Done: tripletEdgeJoinCriterion" << std::endl;

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
        Persistency &persistency,
        CSVWriter &csvWriter,
        int const &index = 5
) {
//    std::cout << "Run: edgeJoinCriterion" << std::endl;
    persistency.edgeJoinCriterion();
//    std::cout << "Done: edgeJoinCriterion" << std::endl;

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
        Persistency &persistency,
        CSVWriter &csvWriter,
        int const &index = 6
) {
//    std::cout << "Run: tripletJoinCriterion" << std::endl;
    persistency.tripletJoinCriterion();
//    std::cout << "Done: tripletJoinCriterion" << std::endl;

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
        Persistency &persistency,
        CSVWriter &csvWriter,
        int const &index = 7
) {
    // cut criteria
//    std::cout << "Run: edgeCutCriterion" << std::endl;
    persistency.edgeCutCriterion();
//    std::cout << "Done: edgeCutCriterion" << std::endl;

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
        Persistency &persistency,
        CSVWriter &csvWriter,
        int const &index = 8
) {
//    std::cout << "Run: tripletCutCriterion" << std::endl;
    persistency.tripletCutCriterion();
//    std::cout << "Done: tripletCutCriterion" << std::endl;

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
bool
doFindIndependentSubsets(
        Persistency &persistency,
        CSVWriter &csvWriter,
        int const &index = 8
) {
//    std::cout << "Run: tripletCutCriterion" << std::endl;
//    std::cout << "Done: tripletCutCriterion" << std::endl;
    bool foundIndependentSubsets = persistency.findIndependentSubProblems();

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

    return foundIndependentSubsets;
}

template<class PROBLEM_VIEW>
inline
void
runEquilateralTrianglesInstanceAllCriteria(
        PROBLEM_VIEW &adaptor,
        CSVWriter &csvWriter,
        std::filesystem::path const &directory = "./test"
) {
    if (!std::filesystem::exists(directory)) {
        std::filesystem::create_directories(directory);
    }

    Persistency persistency(adaptor);

    doInitialSnapshot(persistency, csvWriter, 0);

    // join phase
    bool didJoin = false;
    size_t oldNumberOfNodes = persistency.remainingNodes();
    int index = 1;

    do {
        didJoin = false;

        // all joins applied
        // try to find independent subsets
        doFindIndependentSubsets(persistency, csvWriter, index);
        index++;

        doSubsetJoinCriterion(persistency, csvWriter, index);

        // this ensures restarts of join phase only if a subsequent criterion leads to another join
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


        if (persistency.remainingNodes() != oldNumberOfNodes) {
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
runInstanceSizeExperiments(size_t const numSeeds, size_t const maxN) {

    double costAtDistanceZero = -3;
//     30° dev -> fixed
    double const &distanceAtCostZero = M_PI / 6;

    double const thresholdFactor = 4;

    size_t const &startSeed = 4242;

    std::vector<double> sigmas = {0.01, 0.06, 0.1};

    std::vector<Point2> midPoints = {
            Point2({0, 1}),
            Point2({0.5, 0.55}),
            Point2({1.5, -0.7})
    };
    std::vector<double> radii = {
            1.1 / sqrt(3),
            0.9 / sqrt(3),
            1.0 / sqrt(3)
    };

    std::vector<double> angleOffsets = {
            80 * M_PI / 180,
            -110 * M_PI / 180,
            4 * M_PI_2 / 6
    };

    size_t const numberOfExpectedRemainingNodes = midPoints.size();

    std::vector<size_t> seeds(numSeeds);
    for (size_t seedIndex = 0; seedIndex < numSeeds; ++seedIndex) {
        seeds[seedIndex] = startSeed + seedIndex;
    }

    std::filesystem::path outDir = "equilateral_instance_sizes_"
                                   + std::to_string(numSeeds) + "seeds_"
                                   + to_string_with_precision(distanceAtCostZero * 180 / M_PI, 0) +
                                   "distanceOfCostZero_"
                                   + to_string_with_precision(thresholdFactor, 0) + "sigma";

    // create the directory before the first CSV writer opens a file in it, otherwise the first seed is not written
    std::filesystem::create_directories(outDir);

    auto totalStart = std::chrono::system_clock::now();

    for (size_t n = 2; n <= maxN; ++n) {

        std::vector<std::vector<size_t>> pointsAtTriangleNodes = {
                {n - 2, n - 1, n},
                {n - 1, n,     n - 1},
                {n - 2, n - 2, n}
        };

        size_t numberOfNodes = 9 * n - 9;

        std::cout << "Running with numNodes: " << numberOfNodes << std::endl;

        std::filesystem::path csvFilePath = outDir / (to_string_with_precision(numberOfNodes, 0) + "nodes.csv");

        for (auto &sigma: sigmas) {

            std::cout << "Running with sigma = " << sigma << "... " << std::endl;

            double const &tripletMaxEdgeLengthThreshold = thresholdFactor * sigma;
            double const &edgeLengthThreshold = thresholdFactor * sigma;

            auto innerStart = std::chrono::system_clock::now();
            for (auto &seed: seeds) {

                CSVWriter csvWriter = initializeCSVWriter(
                        sigma,
                        numberOfNodes,
                        numberOfExpectedRemainingNodes,
                        costAtDistanceZero,
                        distanceAtCostZero,
                        edgeLengthThreshold,
                        tripletMaxEdgeLengthThreshold,
                        seed,
                        csvFilePath
                );

                Problem problem(
                        midPoints,
                        radii,
                        angleOffsets,
                        pointsAtTriangleNodes,
                        sigma,
                        seed
                );

                Adaptor<AngleCostMode::SUM_ANGLE_DIFFERENCES_ABS> adaptor(
                        problem,
                        distanceAtCostZero,
                        costAtDistanceZero,
                        costAtDistanceZero,
                        tripletMaxEdgeLengthThreshold,
                        edgeLengthThreshold
                );

                runEquilateralTrianglesInstanceAllCriteria(
                        adaptor,
                        csvWriter,
                        outDir
                );
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

inline
size_t
parseNumber(std::string const &text) {
    std::istringstream ss(text);
    size_t number;
    if (!(ss >> number)) {
        throw std::runtime_error("Invalid number");
    } else if (!ss.eof()) {
        throw std::runtime_error("Trailing characters after number.");
    }
    return number;
}


int main(int argc, char** argv) {
    size_t numSeeds = 30;
    size_t maxN = 13;

    for (int i = 1; i < argc; i += 2) {
        std::string flag = argv[i];
        if (i + 1 >= argc) {
            throw std::runtime_error("Flag " + flag + " needs a number.");
        }
        if (flag == "-s" || flag == "--seeds") {
            numSeeds = parseNumber(argv[i + 1]);
        } else if (flag == "-n" || flag == "--max-n") {
            maxN = parseNumber(argv[i + 1]);
        } else {
            throw std::runtime_error("Flag " + flag + " not allowed. Use --seeds <number> and --max-n <n>.");
        }
    }

    if (numSeeds == 0 || maxN < 2) {
        throw std::runtime_error("Use at least one seed and n >= 2, i.e. instances with 9*n - 9 nodes.");
    }

    std::cout << "Running with " << numSeeds << " seeds and n = 2, ..., " << maxN
              << ", i.e. up to " << 9 * maxN - 9 << " nodes." << std::endl;
    runInstanceSizeExperiments(numSeeds, maxN);

    return 0;
}
