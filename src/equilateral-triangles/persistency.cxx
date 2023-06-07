#include <andres/graph/multicut-cubic/persistency.hxx>
#include <random>
#include <iostream>
#include "dssc/equilateral-triangles/problem.hxx"
#include "dssc/equilateral-triangles/problem-adaptor.hxx"
#include "dssc/csv-writer.hxx"
#include "chrono"
#include "filesystem"
#include "dssc/equilateral-triangles/latex.hxx"
#include "dssc/cost-adjustment-adaptors/problem-rescaling-adaptor.hxx"

typedef dssc::ProblemEquilateralTriangles<double> Problem;
typedef dssc::AngleCostMode AngleCostMode;
template<AngleCostMode AngleCost>
using Adaptor = dssc::ProblemEquilateralTrianglesAdaptor<
        double,
        AngleCost
>;
typedef dssc::CostRescalingAdaptor<double, dssc::RescalingCostType::TRIPLES> Rescaler;


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
std::string
currentDateTime() {
    std::time_t t = std::time(nullptr);
    std::tm *now = std::localtime(&t);

    char buffer[128];
    strftime(buffer, sizeof(buffer), "%m-%d-%Y-%H-%M-%S", now);
    return buffer;
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
runAdjustedEquilateralTrianglesSubsetJoin(
        PROBLEM_VIEW &adaptor,
        CSVWriter &csvWriter,
        std::filesystem::path const &directory = "./test"
) {
    if (!std::filesystem::exists(directory)) {
        std::filesystem::create_directories(directory);
    }

    // subset join - criterion 1
    Persistency persistency(adaptor);

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

template<class PROBLEM_VIEW>
inline
void
runAdjustedEquilateralTrianglesEdgeSubgraph(
        PROBLEM_VIEW &adaptor,
        CSVWriter &csvWriter,
        std::filesystem::path const &directory = "./test"
) {
    if (!std::filesystem::exists(directory)) {
        std::filesystem::create_directories(directory);
    }

    // subset join - criterion 1
    Persistency persistency(adaptor);

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

template<class PROBLEM_VIEW>
inline
void
runAdjustedEquilateralTrianglesTripletSubgraph(
        PROBLEM_VIEW &adaptor,
        CSVWriter &csvWriter,
        std::filesystem::path const &directory = "./test"
) {
    if (!std::filesystem::exists(directory)) {
        std::filesystem::create_directories(directory);
    }

    // subset join - criterion 1
    Persistency persistency(adaptor);

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


template<class PROBLEM_VIEW>
inline
void
runAdjustedEquilateralTrianglesTripletEdgeJoin(
        PROBLEM_VIEW &adaptor,
        CSVWriter &csvWriter,
        std::filesystem::path const &directory = "./test"
) {
    if (!std::filesystem::exists(directory)) {
        std::filesystem::create_directories(directory);
    }

    // subset join - criterion 1
    Persistency persistency(adaptor);

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

template<class PROBLEM_VIEW>
inline
void
runAdjustedEquilateralTrianglesEdgeJoin(
        PROBLEM_VIEW &adaptor,
        CSVWriter &csvWriter,
        std::filesystem::path const &directory = "./test"
) {
    if (!std::filesystem::exists(directory)) {
        std::filesystem::create_directories(directory);
    }

    // subset join - criterion 1
    Persistency persistency(adaptor);


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

template<class PROBLEM_VIEW>
inline
void
runAdjustedEquilateralTrianglesTripletJoin(
        PROBLEM_VIEW &adaptor,
        CSVWriter &csvWriter,
        std::filesystem::path const &directory = "./test"
) {
    if (!std::filesystem::exists(directory)) {
        std::filesystem::create_directories(directory);
    }

    // subset join - criterion 1
    Persistency persistency(adaptor);

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


template<class PROBLEM_VIEW>
inline
void
runAdjustedEquilateralTrianglesEdgeCut(
        PROBLEM_VIEW &adaptor,
        CSVWriter &csvWriter,
        std::filesystem::path const &directory = "./test"
) {
    if (!std::filesystem::exists(directory)) {
        std::filesystem::create_directories(directory);
    }

    // subset join - criterion 1
    Persistency persistency(adaptor);

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

template<class PROBLEM_VIEW>
inline
void
runAdjustedEquilateralTrianglesTripletCut(
        PROBLEM_VIEW &adaptor,
        CSVWriter &csvWriter,
        std::filesystem::path const &directory = "./test"
) {
    if (!std::filesystem::exists(directory)) {
        std::filesystem::create_directories(directory);
    }

    // subset join - criterion 1
    Persistency persistency(adaptor);

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

template<class PROBLEM_VIEW>
inline
void
runAdjustedEquilateralTrianglesFindIndependentSubProblems(
        PROBLEM_VIEW &adaptor,
        CSVWriter &csvWriter,
        std::filesystem::path const &directory = "./test"
) {
    if (!std::filesystem::exists(directory)) {
        std::filesystem::create_directories(directory);
    }

    // subset join - criterion 1
    Persistency persistency(adaptor);

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

template<
        AngleCostMode AngleCost = AngleCostMode::SUM_ANGLE_DIFFERENCES_ABS
>
inline
void
runAllCriteriaExperiments(size_t const n) {

    double costAtDistanceZero = -3;
//     30° dev -> fixed
    double const &distanceAtCostZero = M_PI / 6;

    double maxPossibleCost = 0;
    if (AngleCost == AngleCostMode::SUM_ANGLE_DIFFERENCES_ABS) {
        maxPossibleCost = -costAtDistanceZero / distanceAtCostZero * (4 * M_PI / 3 - distanceAtCostZero);
    } else {
        throw std::runtime_error("unknown angle cost mode.");
    }

    double const thresholdFactor = 4;

    size_t const &startSeed = 4242;
    size_t const &numSeeds = 30;

    std::vector<double> sigmas;

    double const & sigmaMin = 0.00625;
    double const & sigmaStep = 0.00625;
    double const & sigmaMax = 0.1;

    double sigma = sigmaMin;
    while (sigma <= sigmaMax){
        sigmas.push_back(sigma);
        sigma += sigmaStep;
    }

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

    auto totalStart = std::chrono::system_clock::now();

    std::vector<std::vector<size_t>> pointsAtTriangleNodes = {
            {n - 2, n - 1, n},
            {n - 1, n,     n - 1},
            {n - 2, n - 2, n}
    };

    size_t numberOfNodes = 9 * n - 9;

    std::cout << "Running with numNodes: " << numberOfNodes << std::endl;

    std::filesystem::path outDir = "equilateral_"
                                   + to_string_with_precision(numberOfNodes, 0) + "nodes_"
                                   + std::to_string(numSeeds) + "seeds_"
                                   + to_string_with_precision(distanceAtCostZero * 180 / M_PI, 0) +
                                   "distanceOfCostZero_"
                                   + to_string_with_precision(thresholdFactor, 0) + "sigma";

    std::filesystem::path fileName = "experiment.csv";

    std::filesystem::path csvFilePath = outDir / fileName;

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

            Adaptor<AngleCost> baseAdaptor(
                    problem,
                    distanceAtCostZero,
                    costAtDistanceZero,
                    costAtDistanceZero,
                    tripletMaxEdgeLengthThreshold,
                    edgeLengthThreshold
            );

            Rescaler adaptor(
                    baseAdaptor,
                    -costAtDistanceZero,
                    maxPossibleCost,
                    costAtDistanceZero,
                    costAtDistanceZero
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

    auto totalEnd = std::chrono::system_clock::now();

    std::chrono::duration<double> totalDuration = totalEnd - totalStart;

    std::cout << "Total Runtime: " << totalDuration.count() / 60 << "min" << std::endl;
}

inline
size_t
calculateSum(std::vector<std::vector<size_t>> const &array) {
    size_t sum = 0;
    for (auto const &arr: array) {
        for (auto const &el: arr) {
            sum += el;
        }
    }
    return sum;
}


//// runs only one criterion at a time
template<
        AngleCostMode AngleCost = AngleCostMode::SUM_ANGLE_DIFFERENCES_ABS
>
inline
void
runEquilateralTrianglesInstanceIndividualCriteria(
        std::vector<Point2> const &midPoints,
        std::vector<double> const &radii,
        std::vector<double> const &angleOffsets,
        std::vector<std::vector<size_t >> const &pointsAtTriangleNodes,
        double const &sigma,
        double const &distanceAtCostZero,
        double const &costAtDistanceZero,
        double const &tripletMaxEdgeLengthThreshold,
        double const &edgeLengthThreshold,
        size_t const &seed = 4242,
        std::string const &csvBaseFileName = "test",
        std::filesystem::path const &directory = "./test"
) {
    Problem problem(
            midPoints,
            radii,
            angleOffsets,
            pointsAtTriangleNodes,
            sigma,
            seed
    );

    Adaptor<AngleCost> baseAdaptor(
            problem,
            distanceAtCostZero,
            costAtDistanceZero,
            costAtDistanceZero,
            tripletMaxEdgeLengthThreshold,
            edgeLengthThreshold
    );

    double maxPossibleCost = costAtDistanceZero / distanceAtCostZero * (4 * M_PI / 3 - distanceAtCostZero);

    Rescaler adaptor(
            baseAdaptor,
            -costAtDistanceZero,
            maxPossibleCost,
            costAtDistanceZero,
            costAtDistanceZero
    );

    std::string subsetJoinFileName = csvBaseFileName + "_subsetJoin.csv";
    std::filesystem::path subsetJoinPath = directory / subsetJoinFileName;
    CSVWriter csvWriter = initializeCSVWriter(
            sigma,
            calculateSum(pointsAtTriangleNodes),
            midPoints.size(),
            costAtDistanceZero,
            distanceAtCostZero,
            edgeLengthThreshold,
            tripletMaxEdgeLengthThreshold,
            seed,
            subsetJoinPath
    );
    runAdjustedEquilateralTrianglesSubsetJoin(adaptor, csvWriter, directory);

    std::string edgeSubgraphFilename = csvBaseFileName + "_edgeSubgraph.csv";
    std::filesystem::path edgeSubgraphPath = directory / edgeSubgraphFilename;
    csvWriter = initializeCSVWriter(
            sigma,
            calculateSum(pointsAtTriangleNodes),
            midPoints.size(),
            costAtDistanceZero,
            distanceAtCostZero,
            edgeLengthThreshold,
            tripletMaxEdgeLengthThreshold,
            seed,
            edgeSubgraphPath
    );
    runAdjustedEquilateralTrianglesEdgeSubgraph(adaptor, csvWriter, directory);

    std::string tripletEdgeJoinFileName = csvBaseFileName + "_tripletEdgeJoin.csv";
    std::filesystem::path tripletEdgeJoinPath = directory / tripletEdgeJoinFileName;
    csvWriter = initializeCSVWriter(
            sigma,
            calculateSum(pointsAtTriangleNodes),
            midPoints.size(),
            costAtDistanceZero,
            distanceAtCostZero,
            edgeLengthThreshold,
            tripletMaxEdgeLengthThreshold,
            seed,
            tripletEdgeJoinPath
    );
    runAdjustedEquilateralTrianglesTripletEdgeJoin(adaptor, csvWriter, directory);

    std::string tripletSubgraphFilename = csvBaseFileName + "_tripletSubgraph.csv";
    std::filesystem::path tripletSubgraphPath = directory / tripletSubgraphFilename;
    csvWriter = initializeCSVWriter(
            sigma,
            calculateSum(pointsAtTriangleNodes),
            midPoints.size(),
            costAtDistanceZero,
            distanceAtCostZero,
            edgeLengthThreshold,
            tripletMaxEdgeLengthThreshold,
            seed,
            tripletSubgraphPath
    );
    runAdjustedEquilateralTrianglesTripletSubgraph(adaptor, csvWriter, directory);


    std::string edgeJoinFileName = csvBaseFileName + "_edgeJoin.csv";
    std::filesystem::path edgeJoinPath = directory / edgeJoinFileName;
    csvWriter = initializeCSVWriter(
            sigma,
            calculateSum(pointsAtTriangleNodes),
            midPoints.size(),
            costAtDistanceZero,
            distanceAtCostZero,
            edgeLengthThreshold,
            tripletMaxEdgeLengthThreshold,
            seed,
            edgeJoinPath
    );
    runAdjustedEquilateralTrianglesEdgeJoin(adaptor, csvWriter, directory);

    std::string tripletJoinFileName = csvBaseFileName + "_tripletJoin.csv";
    std::filesystem::path tripletJoinPath = directory / tripletJoinFileName;
    csvWriter = initializeCSVWriter(
            sigma,
            calculateSum(pointsAtTriangleNodes),
            midPoints.size(),
            costAtDistanceZero,
            distanceAtCostZero,
            edgeLengthThreshold,
            tripletMaxEdgeLengthThreshold,
            seed,
            tripletJoinPath
    );
    runAdjustedEquilateralTrianglesTripletJoin(adaptor, csvWriter, directory);

    std::string edgeCutFileName = csvBaseFileName + "_edgeCut.csv";
    std::filesystem::path edgeCutPath = directory / edgeCutFileName;
    csvWriter = initializeCSVWriter(
            sigma,
            calculateSum(pointsAtTriangleNodes),
            midPoints.size(),
            costAtDistanceZero,
            distanceAtCostZero,
            edgeLengthThreshold,
            tripletMaxEdgeLengthThreshold,
            seed,
            edgeCutPath
    );
    runAdjustedEquilateralTrianglesEdgeCut(adaptor, csvWriter, directory);

    std::string tripletCutFileName = csvBaseFileName + "_tripletCut.csv";
    std::filesystem::path tripletCutPath = directory / tripletCutFileName;
    csvWriter = initializeCSVWriter(
            sigma,
            calculateSum(pointsAtTriangleNodes),
            midPoints.size(),
            costAtDistanceZero,
            distanceAtCostZero,
            edgeLengthThreshold,
            tripletMaxEdgeLengthThreshold,
            seed,
            tripletCutPath
    );
    runAdjustedEquilateralTrianglesTripletCut(adaptor, csvWriter, directory);

    std::string independentSubsetsFileName = csvBaseFileName + "_findIndependentSubsets.csv";
    std::filesystem::path independentSubsetsPath = directory / independentSubsetsFileName;
    csvWriter = initializeCSVWriter(
            sigma,
            calculateSum(pointsAtTriangleNodes),
            midPoints.size(),
            costAtDistanceZero,
            distanceAtCostZero,
            edgeLengthThreshold,
            tripletMaxEdgeLengthThreshold,
            seed,
            independentSubsetsPath
    );
    runAdjustedEquilateralTrianglesFindIndependentSubProblems(adaptor, csvWriter, directory);
}

inline
void
runIndividualExperiments(size_t const &n) {
    double costAtDistanceZero = -3;
    // 30° dev -> fixed
    double const &distanceAtCostZero = M_PI / 6;

    double const thresholdFactor = 4;

    std::cout << "Running with " << std::to_string(9 * n - 9) << " nodes..." << std::endl;

    size_t const &startSeed = 4242;
    size_t const &numSeeds = 30;

    double const &sigmaMin = 0.00625;
    double const &sigmaStep = 0.00625;
    double const &sigmaMax = 0.15;

    int const &sigmaCount = (sigmaMax - sigmaMin) / sigmaStep + 1;

    std::vector<double> sigmas(sigmaCount);
    for (int i = 0; i < sigmaCount; ++i) {
        sigmas[i] = sigmaMin + i * sigmaStep;
    }
    std::vector<Point2> midPoints = {
            Point2({0, 1}),
            Point2({0.5, 0.55}),
            Point2({1.5, -0.7})
    };
    std::vector<double> radii = {
            1.1 / sqrt(3),
            1.0 / sqrt(3),
            0.9 / sqrt(3)
    };

    std::vector<double> angleOffsets = {
            80 * M_PI / 180,
            -110 * M_PI / 180,
            4 * M_PI_2 / 6
    };

    std::vector<std::vector<size_t>> pointsAtTriangleNodes = {
            {n - 2, n - 1, n},
            {n - 1, n,     n - 1},
            {n - 2, n - 2, n}
    };

    std::vector<size_t> seeds(numSeeds);
    for (size_t seedIndex = 0; seedIndex < numSeeds; ++seedIndex) {
        seeds[seedIndex] = startSeed + seedIndex;
    }

    auto totalStart = std::chrono::system_clock::now();

    std::filesystem::path outDir = "equilateral_individual_"
                                   + to_string_with_precision(9 * n - 9, 0) + "nodes_"
                                   + std::to_string(numSeeds) + "seeds_"
                                   + to_string_with_precision(distanceAtCostZero * 180 / M_PI, 0) +
                                   "distanceOfCostZero_"
                                   + to_string_with_precision(thresholdFactor, 0) + "sigma";
    std::string baseFileName = "experiment";

    for (auto &sigma: sigmas) {
        std::cout << "Running with sigma = " << sigma << "... " << std::endl;

        double const &tripletMaxEdgeLengthThreshold = thresholdFactor * sigma;
        double const &edgeLengthThreshold = thresholdFactor * sigma;

        auto innerStart = std::chrono::system_clock::now();

        for (auto &seed: seeds) {


//                    fileName += ".csv";
            runEquilateralTrianglesInstanceIndividualCriteria(
                    midPoints,
                    radii,
                    angleOffsets,
                    pointsAtTriangleNodes,
                    sigma,
                    distanceAtCostZero,
                    costAtDistanceZero,
                    tripletMaxEdgeLengthThreshold,
                    edgeLengthThreshold,
                    seed,
                    baseFileName,
                    outDir);
        }

        auto innerEnd = std::chrono::system_clock::now();
        std::chrono::duration<double> innerDuration = innerEnd - innerStart;

        std::cout << "Time taken: " << innerDuration.count() << "s" << std::endl;
    }

    auto totalEnd = std::chrono::system_clock::now();

    std::chrono::duration<double> totalDuration = totalEnd - totalStart;

    std::cout << "Total Runtime: " << totalDuration.count() / 60 << "min" << std::endl;
}


int main(int argc, char** argv) {
    size_t const defaultAllN = 6;
    size_t const defaultIndividualN = 6;

    if (argc == 1){
        throw std::runtime_error("Use '-i' or '-a' to run the criteria individually or jointly. "
                                 "Supply a number n as second parameter. "
                                 "This will run the instance with 9*n - 9 nodes. Default is n=6.");
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