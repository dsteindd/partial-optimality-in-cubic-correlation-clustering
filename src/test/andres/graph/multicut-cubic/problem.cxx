#include <stdexcept>

#include "andres/graph/multicut-cubic/problem.hxx"

inline void test(bool const condition) {
    if(!condition) throw std::logic_error("test failed.");
}

template<bool ASSUME_STRICTLY_DECREASING_TRIPLES>
void
testConstructionAndAccess() {
    typedef andres::graph::multicut_cubic::Problem<size_t, ASSUME_STRICTLY_DECREASING_TRIPLES> Problem;

    // define problem
    size_t const numberOfElements = 10;
    Problem problem(numberOfElements);
    test(problem.numberOfElements() == numberOfElements);
    {
        size_t cost = 1;
        for(size_t j = 0; j < numberOfElements; ++j)
        for(size_t k = 0; k < j; ++k)
        for(size_t l = 0; l < k; ++l) {
            problem.costOfTriple({j, k, l}) = cost;
            ++cost;
        }
        for (size_t j = 0; j < numberOfElements; ++j)
        for (size_t k = 0; k < j; ++k) {
            problem.costOfEdge({j, k}) = cost;
            ++cost;
        }
    }

    // test costs
    size_t cost = 1;
    for(size_t j = 0; j < numberOfElements; ++j)
    for(size_t k = 0; k < j; ++k)
    for(size_t l = 0; l < k; ++l) {
        test(problem.costOfTriple({j, k, l}) == cost);
        ++cost;
    }
    for (size_t j = 0; j < numberOfElements; ++j)
    for (size_t k = 0; k < j; ++k) {
        test(problem.costOfEdge({j, k}) == cost);
        ++cost;
    }
}

template<bool ASSUME_STRICTLY_DECREASING_TRIPLES>
void
testContraction() {
    typedef andres::graph::multicut_cubic::Problem<size_t, ASSUME_STRICTLY_DECREASING_TRIPLES> Problem;

    Problem problem(4);
    problem.costOfTriple({2, 1, 0}) = 2;
    problem.costOfTriple({3, 1, 0}) = 3;
    problem.costOfTriple({3, 2, 0}) = 5;
    problem.costOfTriple({3, 2, 1}) = 7;
    problem.costOfEdge({1, 0}) = 9;
    problem.costOfEdge({2, 0}) = 11;
    problem.costOfEdge({2, 1}) = 13;
    problem.costOfEdge({3, 0}) = 15;
    problem.costOfEdge({3, 1}) = 17;
    problem.costOfEdge({3, 2}) = 19;
    problem.contract(2, 1);

    // 0/0 1/{1,2} 2/3

    // 0/0 1/3 2/{1, 2}
    // c_{12}30 = c_{130} + c_{230} = 3 + 5
    // c_{12}3 = c_13 + c_23 + c_123 = 17 + 7 + 19 = 43
    // c_{12}0 = c_{120} + c_{10} + c_{20} = 2 + 9 + 11 = 22


    test(problem.numberOfElements() == 3);

    test(problem.costOfTriple({2, 1, 0}) == 8);
    test(problem.costOfEdge({2, 1}) == 43);
    test(problem.costOfEdge({1, 0}) == 22);
    test(problem.costOfEdge({2, 0}) == 15);
    test(problem.costConstant() == 13);


//    test(problem.costOfTriple({2, 1, 0}) == 17);

    const andres::Partition<size_t> partition = problem.partition();
    test(partition.numberOfSets() == 3);
    test(partition.numberOfElements() == 4);

    std::vector<int> r(static_cast<std::size_t>(partition.numberOfElements()));

    partition.elementLabeling(r.begin());

    const std::array<size_t, 4> expectedLabels = {0, 1, 1, 2};

    for (size_t i = 0; i < expectedLabels.size(); i++){
        test(r[i] == expectedLabels[i]);
    }
}

template<bool ASSUME_STRICTLY_DECREASING_TRIPLES>
void
testContractionLastNode() {
    typedef andres::graph::multicut_cubic::Problem<size_t, ASSUME_STRICTLY_DECREASING_TRIPLES> Problem;

    Problem problem(4);
    problem.costOfTriple({2, 1, 0}) = 2;
    problem.costOfTriple({3, 1, 0}) = 3;
    problem.costOfTriple({3, 2, 0}) = 5;
    problem.costOfTriple({3, 2, 1}) = 7;
    problem.costOfEdge({1, 0}) = 9;
    problem.costOfEdge({2, 0}) = 11;
    problem.costOfEdge({2, 1}) = 13;
    problem.costOfEdge({3, 0}) = 15;
    problem.costOfEdge({3, 1}) = 17;
    problem.costOfEdge({3, 2}) = 19;
    problem.contract(3, 1);

    test(problem.numberOfElements() == 3);

    // 0/0 1/{1, 3} 2/2
    // c_{1,3}20 = c_210 + c_320 = 2 + 5
    // c_{1, 3}2 = c_321 + c_21 + c_32 = 7 + 13 + 19 = 39
    // c_{1, 3}0 = c_310 + c_30 + c_10 = 3 + 15 + 9 = 27;

    test(problem.costOfTriple({2, 1, 0}) == 7);
    test(problem.costOfEdge({1, 0}) == 27);
    test(problem.costOfEdge({2, 1}) == 39);
    test(problem.costOfEdge({2, 0}) == 11);
    test(problem.costConstant() == 17);
    // sum all = 7 + 27 + 39 + 17 = 101

//    test(problem.costOfTriple({2, 1, 0}) == 17);

    const andres::Partition<size_t> partition = problem.partition();
    test(partition.numberOfSets() == 3);
    test(partition.numberOfElements() == 4);

    std::vector<int> r(static_cast<std::size_t>(partition.numberOfElements()));

    partition.elementLabeling(r.begin());

    const std::array<size_t, 4> expectedLabels = {0, 1, 2, 1};

    test(r[0] != r[1]);
    test(r[0] != r[2]);
    test(r[0] != r[3]);
    test(r[1] != r[2]);
    test(r[1] == r[3]);
    test(r[2] != r[3]);
}

template<bool ASSUME_STRICTLY_DECREASING_TRIPLES>
void testContractionTwoByTwoNodes(){
    typedef andres::graph::multicut_cubic::Problem<size_t, ASSUME_STRICTLY_DECREASING_TRIPLES> Problem;

    Problem problem(4);
    problem.costOfTriple({2, 1, 0}) = 2;
    problem.costOfTriple({3, 1, 0}) = 3;
    problem.costOfTriple({3, 2, 0}) = 5;
    problem.costOfTriple({3, 2, 1}) = 7;
    problem.costOfEdge({1, 0}) = 9;
    problem.costOfEdge({2, 0}) = 11;
    problem.costOfEdge({2, 1}) = 13;
    problem.costOfEdge({3, 0}) = 15;
    problem.costOfEdge({3, 1}) = 17;
    problem.costOfEdge({3, 2}) = 19;
    problem.contract(1, 0);
    // contraction makes node s/0 -> last, and node t/1 to the merge node
    // {0, 1}/0 3/1 2/2
    // c_32{01} = c_320 + c_321 = 5 + 7 = 12 != c'_210
    // c_2{0,1} = c_20 + c_21 + c_210 = 11 + 13 + 2 = 26 = c'_21
    // c_3{0,1} = c_310 + c_30 + c_31 = 3 + 15 + 17 = 35 = c'_10

    test(problem.numberOfElements() == 3);
    test(problem.costOfTriple({2, 1, 0}) == 12);
    test(problem.costOfEdge({2, 0}) == 26);
    test(problem.costOfEdge({1, 0}) == 35);
    test(problem.costOfEdge({2, 1}) == 19);
    test(problem.costConstant() == 9);

    // check costs...


    problem.contract(2, 1);

    // 0/0 {2, 1}/1
    // c_{2, 1}0 = c_20 + c_10 + c_210 = 26 + 35 + 12 = 73
    test(problem.numberOfElements() == 2);
    test(problem.costOfEdge({1, 0}) == 73);
    test(problem.costConstant() == 28);

    const andres::Partition<size_t> partition = problem.partition();
    test(partition.numberOfSets() == 2);
    test(partition.numberOfElements() == 4);

    std::vector<int> r(static_cast<std::size_t>(partition.numberOfElements()));

    partition.elementLabeling(r.begin());

    const std::array<size_t, 4> expectedLabels = {0, 0, 1, 1};

    test(r[0] == r[1]);
    test(r[2] == r[3]);
    test(r[0] != r[2]);

}

template<bool ASSUME_STRICTLY_DECREASING_TRIPLES>
void testContractionThreeNodes(){
    typedef andres::graph::multicut_cubic::Problem<size_t, ASSUME_STRICTLY_DECREASING_TRIPLES> Problem;

    Problem problem(4);
    problem.costOfTriple({2, 1, 0}) = 2;
    problem.costOfTriple({3, 1, 0}) = 3;
    problem.costOfTriple({3, 2, 0}) = 5;
    problem.costOfTriple({3, 2, 1}) = 7;
    problem.costOfEdge({1, 0}) = 9;
    problem.costOfEdge({2, 0}) = 11;
    problem.costOfEdge({2, 1}) = 13;
    problem.costOfEdge({3, 0}) = 15;
    problem.costOfEdge({3, 1}) = 17;
    problem.costOfEdge({3, 2}) = 19;
    problem.contract(1, 0);
    // contraction makes node s/0 -> last, and node t/1 to the merge node
    // {0, 1}/0 3/1 2/2
    // c_32{01} = c_320 + c_321 = 5 + 7 = 12 != c'_210
    // c_2{0,1} = c_20 + c_21 + c_210 = 11 + 13 + 2 = 26 = c'_21
    // c_3{0,1} = c_310 + c_30 + c_31 = 3 + 15 + 17 = 35 = c'_10

    test(problem.numberOfElements() == 3);
    test(problem.costOfTriple({2, 1, 0}) == 12);
    test(problem.costOfEdge({2, 0}) == 26);
    test(problem.costOfEdge({1, 0}) == 35);
    test(problem.costOfEdge({2, 1}) == 19);
    test(problem.costConstant() == 9);
    // contraction makes node s/0 -> 3, and node t/1 to the merge node
    problem.contract(2, 0);
    // 0/{0, 2} 1/3
    // c_{0, 2}1 = c_21 + c_10 + c_210 = 12 + 26 + 35 = 73
    test(problem.costOfEdge({1, 0}) == 66);
    test(problem.costConstant() == 35);


    test(problem.numberOfElements() == 2);

    const andres::Partition<size_t> partition = problem.partition();
    test(partition.numberOfSets() == 2);
    test(partition.numberOfElements() == 4);

    std::vector<int> r(static_cast<std::size_t>(partition.numberOfElements()));

    partition.elementLabeling(r.begin());

    const std::array<size_t, 4> expectedLabels = {0, 0, 0, 1};

    test(r[0] == r[1] && r[1] == r[2]);
    test(r[2] != r[3]);
}

template<bool ASSUME_STRICTLY_DECREASING_TRIPLES>
void testContractionThreeNodesWithLastNodeFirst(){
    typedef andres::graph::multicut_cubic::Problem<size_t, ASSUME_STRICTLY_DECREASING_TRIPLES> Problem;

    Problem problem(5);
    // 5*4*3 / 6 = 10 costs
    problem.costOfTriple({2, 1, 0}) = 2;
    problem.costOfTriple({3, 1, 0}) = 3;
    problem.costOfTriple({3, 2, 0}) = 5;
    problem.costOfTriple({3, 2, 1}) = 7;
    problem.costOfTriple({4, 1, 0}) = 9;
    problem.costOfTriple({4, 2, 0}) = 11;
    problem.costOfTriple({4, 2, 1}) = 13;
    problem.costOfTriple({4, 3, 0}) = 15;
    problem.costOfTriple({4, 3, 1}) = 17;
    problem.costOfTriple({4, 3, 2}) = 19;
    problem.costOfEdge({1, 0}) = 21;
    problem.costOfEdge({2, 0}) = 23;
    problem.costOfEdge({2, 1}) = 25;
    problem.costOfEdge({3, 0}) = 27;
    problem.costOfEdge({3, 1}) = 29;
    problem.costOfEdge({3, 2}) = 31;
    problem.costOfEdge({4, 0}) = 33;
    problem.costOfEdge({4, 1}) = 35;
    problem.costOfEdge({4, 2}) = 37;
    problem.costOfEdge({4, 3}) = 39;
    // 0/0, 1/1, 2/2, 3/3 4/4
    problem.contract(4, 0);
    // 0/{0, 4} 1/1 2/2 3/3
    // c_{0,4}32 = c_432 + c_320 = 19 + 5 = 24
    // c_{0,4}31 = c_431 + c_310 = 17 + 3 = 20
    // c_{0,4}21 = c_{421} + c_210 = 13 + 2 = 15
    // c_{0,4}3 = c_43 + c_30 + c_430 = 39 + 27 + 15 = 81
    // c_{0,4}2 = c_42 + c_20 + c_420 = 37 + 23 + 11 = 71
    // c_{0, 4}1 = c_41 + c_10 + c_410 = 35 + 21 + 9 = 65

    test(problem.costOfTriple({3, 2, 1}) == 7);
    test(problem.costOfTriple({3, 2, 0}) == 24);
    test(problem.costOfTriple({3, 1, 0}) == 20);
    test(problem.costOfTriple({2, 1, 0}) == 15);
    test(problem.costOfEdge({3, 0}) == 81);
    test(problem.costOfEdge({2, 0}) == 71);
    test(problem.costOfEdge({1, 0}) == 65);
    test(problem.costOfEdge({3, 2}) == 31);
    test(problem.costOfEdge({3, 1}) == 29);
    test(problem.costOfEdge({2, 1}) == 25);
    test(problem.costConstant() == 33);


    problem.contract(3, 0);
//    // 0/{0, 3} 1/1 2/2
    // c_21{0, 3} = c_210 + c_321 = 15 + 7 = 22
    // c_2{0, 3} = c_20 + c_32 + c_320 = 71 + 31 + 24 = 126
    // c_1{0, 3} = c_31 + c_10 + c_310 = 29 + 65 + 20 = 114
    test(problem.costOfTriple({2, 1, 0}) == 22);
    test(problem.costOfEdge({2, 0}) == 126);
    test(problem.costOfEdge({1, 0}) == 114);
    test(problem.costOfEdge({2, 1}) == 25);
    test(problem.costConstant() == 114);


    test(problem.numberOfElements() == 3);

    const andres::Partition<size_t> partition = problem.partition();
    test(partition.numberOfSets() == 3);
    test(partition.numberOfElements() == 5);

    std::vector<int> r(static_cast<std::size_t>(partition.numberOfElements()));

    partition.elementLabeling(r.begin());

    const std::array<size_t, 5> expectedLabels = {0, 1, 2, 0, 0};

    test((r[0] == r[3]) && (r[3] == r[4]));
    test(r[0] != r[1]);
    test(r[0] != r[2]);
    test(r[1] != r[2]);
}


template<bool ASSUME_STRICTLY_DECREASING_TRIPLES>
void testTripleContraction(){
    typedef andres::graph::multicut_cubic::Problem<size_t, ASSUME_STRICTLY_DECREASING_TRIPLES> Problem;

    Problem problem(5);
    // 5*4*3 / 6 = 10 costs
    problem.costOfTriple({2, 1, 0}) = 2;
    problem.costOfTriple({3, 1, 0}) = 3;
    problem.costOfTriple({3, 2, 0}) = 5;
    problem.costOfTriple({3, 2, 1}) = 7;
    problem.costOfTriple({4, 1, 0}) = 9;
    problem.costOfTriple({4, 2, 0}) = 11;
    problem.costOfTriple({4, 2, 1}) = 13;
    problem.costOfTriple({4, 3, 0}) = 15;
    problem.costOfTriple({4, 3, 1}) = 17;
    problem.costOfTriple({4, 3, 2}) = 19;
    problem.costOfEdge({1, 0}) = 21;
    problem.costOfEdge({2, 0}) = 23;
    problem.costOfEdge({2, 1}) = 25;
    problem.costOfEdge({3, 0}) = 27;
    problem.costOfEdge({3, 1}) = 29;
    problem.costOfEdge({3, 2}) = 31;
    problem.costOfEdge({4, 0}) = 33;
    problem.costOfEdge({4, 1}) = 35;
    problem.costOfEdge({4, 2}) = 37;
    problem.costOfEdge({4, 3}) = 39;


    problem.contract(4, 3, 0);
    // 0/{0, 3, 4} 1/1 2/2
    // c_210 = c_421 + c_321 + c_210 = 13 + 7 + 2 = 22
    test(problem.costOfTriple({2, 1, 0}) == 22);
    // c_20 = c_42 + c_32 + c_20 + c_432 + c_402 + c_302 = 37 + 31 + 23 + 19 + 11 + 5 = 126
    test(problem.costOfEdge({2, 0}) == 126);
    test(problem.costOfEdge({1, 0}) == 114);
    test(problem.costOfEdge({2, 1}) == 25);

    test(problem.costConstant() == 114);

    test(problem.numberOfElements() == 3);

    const andres::Partition<size_t> partition = problem.partition();
    test(partition.numberOfSets() == 3);
    test(partition.numberOfElements() == 5);

    std::vector<int> r(static_cast<std::size_t>(partition.numberOfElements()));

    partition.elementLabeling(r.begin());

    const std::array<size_t, 5> expectedLabels = {0, 1, 2, 0, 0};

    test((r[0] == r[3]) && (r[3] == r[4]));
    test(r[0] != r[1]);
    test(r[1] != r[2]);
    test(r[0] != r[2]);
}

template<bool ASSUME_STRICTLY_DECREASING_TRIPLES>
void testNegativePartsProblem(){
    typedef andres::graph::multicut_cubic::Problem<double, ASSUME_STRICTLY_DECREASING_TRIPLES> Problem;

    Problem problem(4);
    problem.costOfTriple({2, 1, 0}) = -2;
    problem.costOfTriple({3, 1, 0}) = 3;
    problem.costOfTriple({3, 2, 0}) = -5;
    problem.costOfTriple({3, 2, 1}) = 7;
    problem.costOfEdge({1, 0}) = 9;
    problem.costOfEdge({2, 0}) = -11;
    problem.costOfEdge({2, 1}) = 13;
    problem.costOfEdge({3, 0}) = -15;
    problem.costOfEdge({3, 1}) = 17;
    problem.costOfEdge({3, 2}) = -19;


    Problem onlyNegatives = problem.onlyNegativePartCosts();

    test(onlyNegatives.costOfTriple({2, 1, 0}) == 2);
    test(onlyNegatives.costOfTriple({3, 1, 0}) == 0);
    test(onlyNegatives.costOfTriple({3, 2, 0}) == 5);
    test(onlyNegatives.costOfTriple({3, 2, 1}) == 0);
    test(onlyNegatives.costOfEdge({1, 0}) == 0);
    test(onlyNegatives.costOfEdge({2, 0}) == 11);
    test(onlyNegatives.costOfEdge({2, 1}) == 0);
    test(onlyNegatives.costOfEdge({3, 0}) == 15);
    test(onlyNegatives.costOfEdge({3, 1}) == 0);
    test(onlyNegatives.costOfEdge({3, 2}) == 19);

}

template<bool ASSUME_STRICTLY_DECREASING_TRIPLES>
void testAbsoluteCostsProblem(){
    typedef andres::graph::multicut_cubic::Problem<double, ASSUME_STRICTLY_DECREASING_TRIPLES> Problem;

    Problem problem(4);
    problem.costOfTriple({2, 1, 0}) = -2;
    problem.costOfTriple({3, 1, 0}) = 3;
    problem.costOfTriple({3, 2, 0}) = -5;
    problem.costOfTriple({3, 2, 1}) = 7;
    problem.costOfEdge({1, 0}) = 9;
    problem.costOfEdge({2, 0}) = -11;
    problem.costOfEdge({2, 1}) = 13;
    problem.costOfEdge({3, 0}) = -15;
    problem.costOfEdge({3, 1}) = 17;
    problem.costOfEdge({3, 2}) = -19;


    Problem onlyNegatives = problem.absoluteCosts();

    test(onlyNegatives.costOfTriple({2, 1, 0}) == 2);
    test(onlyNegatives.costOfTriple({3, 1, 0}) == 3);
    test(onlyNegatives.costOfTriple({3, 2, 0}) == 5);
    test(onlyNegatives.costOfTriple({3, 2, 1}) == 7);
    test(onlyNegatives.costOfEdge({1, 0}) == 9);
    test(onlyNegatives.costOfEdge({2, 0}) == 11);
    test(onlyNegatives.costOfEdge({2, 1}) == 13);
    test(onlyNegatives.costOfEdge({3, 0}) == 15);
    test(onlyNegatives.costOfEdge({3, 1}) == 17);
    test(onlyNegatives.costOfEdge({3, 2}) == 19);

}

template<bool ASSUME_STRICTLY_DECREASING_TRIPLES>
void testProblemProjection(){
    typedef andres::graph::multicut_cubic::Problem<double, ASSUME_STRICTLY_DECREASING_TRIPLES> Problem;

    {
        Problem problem(4);
        problem.costOfTriple({2, 1, 0}) = -2;
        problem.costOfTriple({3, 1, 0}) = 3;
        problem.costOfTriple({3, 2, 0}) = -5;
        problem.costOfTriple({3, 2, 1}) = 7;
        problem.costOfEdge({1, 0}) = 9;
        problem.costOfEdge({2, 0}) = -11;
        problem.costOfEdge({2, 1}) = 13;
        problem.costOfEdge({3, 0}) = -15;
        problem.costOfEdge({3, 1}) = 17;
        problem.costOfEdge({3, 2}) = -19;


        Problem projectedProblem = problem.orderedProjectTo({1, 0});

        test(projectedProblem.numberOfElements() == 2);
        test(projectedProblem.costOfEdge({1, 0}) == problem.costOfEdge({1, 0}));

        Problem unorderedProjectionProblem = problem.unorderedProjectTo({1, 0});

        test(unorderedProjectionProblem.numberOfElements() == 2);
        test(unorderedProjectionProblem.costOfEdge({1, 0}) == 9);

    }

    {
        Problem problem(4);
        problem.costOfTriple({2, 1, 0}) = -2;
        problem.costOfTriple({3, 1, 0}) = 3;
        problem.costOfTriple({3, 2, 0}) = -5;
        problem.costOfTriple({3, 2, 1}) = 7;
        problem.costOfEdge({1, 0}) = 9;
        problem.costOfEdge({2, 0}) = -11;
        problem.costOfEdge({2, 1}) = 13;
        problem.costOfEdge({3, 0}) = -15;
        problem.costOfEdge({3, 1}) = 17;
        problem.costOfEdge({3, 2}) = -19;


        Problem projectedProblem = problem.orderedProjectTo({3, 1, 0});

        test(projectedProblem.numberOfElements() == 3);
        test(projectedProblem.costOfTriple({2, 1, 0}) == 3);
        test(projectedProblem.costOfEdge({2, 1}) == 9);
        test(projectedProblem.costOfEdge({2, 0}) == -15);
        test(projectedProblem.costOfEdge({1, 0}) == 17);

        Problem unorderedProjectionProblem = problem.unorderedProjectTo({3, 1, 0});
        test(unorderedProjectionProblem.numberOfElements() == 3);
        test(unorderedProjectionProblem.costOfTriple({2, 1, 0}) == 3);

        double edge10 = unorderedProjectionProblem.costOfEdge({1, 0});
        double edge20 = unorderedProjectionProblem.costOfEdge({2, 0});
        double edge21 = unorderedProjectionProblem.costOfEdge({2, 1});

        test(edge10 != edge20 && edge20 != edge21 && edge10 != edge21);

        test(edge21 == 9 || edge21 == -15 || edge21 == 17);
        test(edge20 == 9 || edge20 == -15 || edge20 == 17);
        test(edge10 == 9 || edge10 == -15 || edge10 == 17);
    }
}


template<bool ASSUME_STRICTLY_DECREASING>
inline
void testSubsetContraction(){
    typedef andres::graph::multicut_cubic::Problem<double, ASSUME_STRICTLY_DECREASING> Problem;

    {
        Problem problem(4);
        problem.costOfTriple({2, 1, 0}) = 2;
        problem.costOfTriple({3, 1, 0}) = 3;
        problem.costOfTriple({3, 2, 0}) = 5;
        problem.costOfTriple({3, 2, 1}) = 7;
        problem.costOfEdge({1, 0}) = 9;
        problem.costOfEdge({2, 0}) = 11;
        problem.costOfEdge({2, 1}) = 13;
        problem.costOfEdge({3, 0}) = 15;
        problem.costOfEdge({3, 1}) = 17;
        problem.costOfEdge({3, 2}) = 19;
        problem.contract({1, 0});

        test(problem.numberOfElements() == 3);

        // before swap: 0 1 2 3
        // after swap 1 <-> last: 0 3 2 1
        // after swap 0 <-> last - 1: 2 3 0 1
        // after merge: 2 3 {0, 1}

        // expected: c(1, 0)
        test(problem.costConstant() == 9);
        // expected: c(1, 0) = c(3, 2)
        test(problem.costOfEdge({1, 0}) == 19);
        // expected: (U, 2) = c210 + c20 + c21 = 2 + 11 + 13 = 26
        test(problem.costOfEdge({2, 0}) == 26);
        // expeected: (U, 3) = c310 + c31 + c30 = 3 + 17 + 15 = 35
        test(problem.costOfEdge({2, 1}) == 35);
        // expected: (U, 2, 3) = c320 + c321 = 5 + 7
        test(problem.costOfTriple({2, 1, 0}) == 12);

        // contract again
        problem.contract({1, 0});

        test(problem.numberOfElements() == 2);

        // before swap: 0 1 2
        // after swap 1 <-> last: 2 1 0
        // after swap 1 <-> last - 1: 2 1 0
        // after merge: 2 {1, 0}

        // expected additive constant c(1, 0)
        test(problem.costConstant() == 28);
        // expected c(1, 0) = c210 + c21 + c20 = 12 + 35 + 26 = 73
        test(problem.costOfEdge({1, 0}) == 73);

        // contract again
        problem.contract({1, 0});
        test(problem.numberOfElements() == 1);
        test(problem.costConstant() == 101);
    }

    {
        Problem problem(4);
        problem.costOfTriple({2, 1, 0}) = 2;
        problem.costOfTriple({3, 1, 0}) = 3;
        problem.costOfTriple({3, 2, 0}) = 5;
        problem.costOfTriple({3, 2, 1}) = 7;
        problem.costOfEdge({1, 0}) = 9;
        problem.costOfEdge({2, 0}) = 11;
        problem.costOfEdge({2, 1}) = 13;
        problem.costOfEdge({3, 0}) = 15;
        problem.costOfEdge({3, 1}) = 17;
        problem.costOfEdge({3, 2}) = 19;
        problem.contract({3, 2, 1, 0});

        test(problem.numberOfElements() == 1);
        test(problem.costConstant() == 101);
    }

    {
        Problem problem(4);
        problem.costOfTriple({2, 1, 0}) = 2;
        problem.costOfTriple({3, 1, 0}) = 3;
        problem.costOfTriple({3, 2, 0}) = 5;
        problem.costOfTriple({3, 2, 1}) = 7;
        problem.costOfEdge({1, 0}) = 9;
        problem.costOfEdge({2, 0}) = 11;
        problem.costOfEdge({2, 1}) = 13;
        problem.costOfEdge({3, 0}) = 15;
        problem.costOfEdge({3, 1}) = 17;
        problem.costOfEdge({3, 2}) = 19;
        problem.contract({1, 0});

        test(problem.numberOfElements() == 3);

        // before swap: 0 1 2 3
        // after swap 1 <-> last: 0 3 2 1
        // after swap 0 <-> last - 1: 2 3 0 1
        // after merge: 2 3 {0, 1}

        // expected: c(1, 0)
        test(problem.costConstant() == 9);
        // expected: c(1, 0) = c(3, 2)
        test(problem.costOfEdge({1, 0}) == 19);
        // expected: (U, 2) = c210 + c20 + c21 = 2 + 11 + 13 = 26
        test(problem.costOfEdge({2, 0}) == 26);
        // expeected: (U, 3) = c310 + c31 + c30 = 3 + 17 + 15 = 35
        test(problem.costOfEdge({2, 1}) == 35);
        // expected: (U, 2, 3) = c320 + c321 = 5 + 7
        test(problem.costOfTriple({2, 1, 0}) == 12);

        // contract again
        problem.contract({2, 1});

        test(problem.numberOfElements() == 2);

        // before swap: 0 1 2
        // after swap 2 <-> last: 0 1 2
        // after swap 1 <-> last - 1: 0 1 2
        // after merge: 0 {1, 2}
        // back track: 2 {3, 0, 1}

        // expected additive constant c(1, 2)
        test(problem.costConstant() == 9 + 35);
        // expected c(1, 0) = c210 + c10 + c20 = 12 + 19 + 26 = 57
        test(problem.costOfEdge({1, 0}) == 57);

        // contract again
        problem.contract({1, 0});
        test(problem.numberOfElements() == 1);
        test(problem.costConstant() == 101);
    }

    {
        Problem problem(4);
        problem.costOfTriple({2, 1, 0}) = 2;
        problem.costOfTriple({3, 1, 0}) = 3;
        problem.costOfTriple({3, 2, 0}) = 5;
        problem.costOfTriple({3, 2, 1}) = 7;
        problem.costOfEdge({1, 0}) = 9;
        problem.costOfEdge({2, 0}) = 11;
        problem.costOfEdge({2, 1}) = 13;
        problem.costOfEdge({3, 0}) = 15;
        problem.costOfEdge({3, 1}) = 17;
        problem.costOfEdge({3, 2}) = 19;
        problem.contract({3, 0, 1});

        test(problem.numberOfElements() == 2);

        // before swap: 0 1 2 3
        // after swap 3 <-> last: 0 1 2 3
        // after swap 0 <-> last - 1: 2 3 0 1
        // after swap 1 <-> last - 2: 2 3 0 1
        // after merge: 2 {3, 0, 1}

        // new ordered --> same costs
        // before swap: 0 1 2 3
        // after swap 3 <-> last: 0 1 2 3
        // after swap 1 <-> last - 1: 0 2 1 3
        // after swap 0 <-> last - 2: 2 0 1 3
        // after merge: 2 {0, 1, 3}

        // expected additive constant c(1, 2)
        test(problem.costConstant() == 9 + 35);
        // expected c(1, 0) = c210 + c10 + c20 = 12 + 19 + 26 = 57
        test(problem.costOfEdge({1, 0}) == 57);

        // contract again
        problem.contract({1, 0});
        test(problem.numberOfElements() == 1);
        test(problem.costConstant() == 101);
    }

    {
        Problem problem(4);
        problem.costOfTriple({2, 1, 0}) = 2;
        problem.costOfTriple({3, 1, 0}) = 3;
        problem.costOfTriple({3, 2, 0}) = 5;
        problem.costOfTriple({3, 2, 1}) = 7;
        problem.costOfEdge({1, 0}) = 9;
        problem.costOfEdge({2, 0}) = 11;
        problem.costOfEdge({2, 1}) = 13;
        problem.costOfEdge({3, 0}) = 15;
        problem.costOfEdge({3, 1}) = 17;
        problem.costOfEdge({3, 2}) = 19;

        problem.contract({1, 3});

        test(problem.numberOfElements() == 3);

        // before swap: 0 1 2 3 (NEW - correct with ordering in swapping in descending order)
        // after swap 3 <-> last: 0 1 2 3
        // after swap 1 <-> last - 1: 0 2 1 3
        // after merge: 0 2 {1, 3}

        // expected additive constant c(1, 3)
        test(problem.costConstant() == 17);
        // expected c(1, 0) = c02 = 11
        test(problem.costOfEdge({1, 0}) == 11);
        // expected c(2, 0) = c310 + c30 + c10 = 3 + 15 + 9 = 27
        test(problem.costOfEdge({2, 0}) == 27);
        // expected c(2, 1) = c321 + c32 + c21 = 7 + 19 + 13 = 39
        test(problem.costOfEdge({2, 1}) == 39);
        // expected c(2, 1, 0) = c320 + c210 = 5 + 2 = 7
        test(problem.costOfTriple({2, 1, 0}) == 7);
    }
}


int main() {
    testConstructionAndAccess<false>();
    testConstructionAndAccess<true>();
    testContraction<false>();
    testContraction<true>();
    testContractionLastNode<false>();
    testContractionLastNode<true>();
    testContractionTwoByTwoNodes<true>();
    testContractionTwoByTwoNodes<false>();
    testContractionThreeNodes<true>();
    testContractionThreeNodes<false>();
    testContractionThreeNodesWithLastNodeFirst<true>();
    testContractionThreeNodesWithLastNodeFirst<false>();
    testTripleContraction<false>();
    testTripleContraction<true>();
    testNegativePartsProblem<false>();
    testNegativePartsProblem<true>();
    testAbsoluteCostsProblem<true>();
    testAbsoluteCostsProblem<false>();
    testSubsetContraction<true>();
    testSubsetContraction<false>();
    testProblemProjection<true>();
    testProblemProjection<false>();

    return 0;
}

