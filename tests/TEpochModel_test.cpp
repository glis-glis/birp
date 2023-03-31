//
// Created by liam on 12/12/22.
//

#include "../stattools/tests/commonTestingFunctions.h"
#include "TBirpTests.h"
#include "TEpochModel.h"
#include "gtest/gtest.h"

//----------------------------------
// TPoissonEpochModel
//----------------------------------

class TPoissonEpochModelBridge : public TPoissonEpochModel {
    // Dummy class inheriting from TEpochModel: allows for testing protected functions
public:
    TPoissonEpochModelBridge(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                             stattools::TParameterTyped<TypeTime, 1> *TimesOfChange)
            : TPoissonEpochModel(Gamma, TimesOfChange) {}

    std::vector<size_t> _getEpochStarts(const TUniqueContainer<size_t> &Timepoints) {
        return TPoissonEpochModel::_getEpochStarts(Timepoints);
    }

    std::pair<std::vector<double>, double> _calculateLLandJP(const std::vector<TMethods> &Methods) const {
        return TPoissonEpochModel::_calculateLLandJP(Methods);
    }
};

class TEpochModelPoissonTest : public testing::Test {
    // class for google test suite: define variables needed for testing TBirpPrior
public:
    // define parameters for test
    std::unique_ptr<stattools::TParameterTyped<TypeGamma, 1>> gamma;
    std::shared_ptr<TPoissonEpochModelBridge> epochModel;
    std::unique_ptr<stattools::TParameterTyped<TypeTime, 1>> timesOfChange;

    void initialize(const std::vector<size_t> &TimesOfChange, const TUniqueContainer<size_t> &Timepoints,
                    const std::vector<TMethods> &Methods,
                    std::vector<double> NewGamma = {}, std::vector<double> OldGamma = {}) {
        // this function is always called whenever a test is called TEST_F(TBirpPriorTest, ...)
        using namespace stattools::prior;

        const size_t numEpochs = TimesOfChange.size() + 1;

        // create a gamma
        coretools::instances::parameters().addParameter("timesOfChange",
                                                        coretools::str::concatenateString(TimesOfChange, ","));
        if (OldGamma.empty()) { OldGamma.resize(numEpochs, 0.0); }
        if (NewGamma.empty()) { NewGamma.resize(numEpochs, 0.0); }
        testing::initializeStorageParameter(gamma, "gamma", {numEpochs}, OldGamma, NewGamma);
        testing::initializeStorageParameter(timesOfChange, "timesOfChange", {numEpochs - 1}, TimesOfChange,
                                            TimesOfChange);

        epochModel =
                std::make_shared<TPoissonEpochModelBridge>(gamma.get(),
                                                           timesOfChange.get()); // get(): unique pointer gets raw pointer
        epochModel->initializeInferred(Timepoints, {}, Methods);
    }
};

TEST_F(TEpochModelPoissonTest, _fillRho_easy1) {
    TBirpCoreBridge test;
    test._readFile("../tests/test_files/method_lookingForWolves.csv");
    test._sortTimepoints();
    auto times = test.timepoints();
    std::vector<size_t> timesOfChange = {2018};
    initialize(timesOfChange, times, test.methods());

    // define expected 1D output vector (one row is one epoch)
    std::vector<std::vector<size_t>> outExp = {
            {2015, 2016, 2017, 2018, 2018},
            {0,    0,    0,    1,    2}};
    // test if result is as expected
    for (size_t i = 0; i < times.size(); i++) {
        for (size_t j = 0; j < timesOfChange.size() + 1; j++) { EXPECT_EQ(epochModel->rho(i, j), outExp.at(j).at(i)); }
    }
}


TEST_F(TEpochModelPoissonTest, _fillRho_easy2) {
    TBirpCoreBridge test;
    test._readFile("../tests/test_files/method_lookingForFloesselhecht.csv");
    test._sortTimepoints();
    auto times = test.timepoints();
    std::vector<size_t> timesOfChange = {2002};
    initialize(timesOfChange, times, test.methods());

    // define expected 1D output vector (one row is one epoch)
    std::vector<std::vector<size_t>> outExp = {{2000, 2001, 2002, 2002, 2002},
                                               {0,    0,    0,    2,    3}};
    // test if result is as expected
    for (size_t i = 0; i < times.size(); i++) {
        for (size_t j = 0; j < timesOfChange.size() + 1; j++) { EXPECT_EQ(epochModel->rho(i, j), outExp.at(j).at(i)); }
    }
}

TEST_F(TEpochModelPoissonTest, _fillRho_TimesOfChange_same_as_times) {
    TBirpCoreBridge test;
    test._readFile("../tests/test_files/method_fossilExtraction_TOC_same_as_times.csv");
    auto times = test.timepoints();
    test._sortTimepoints();
    std::vector<size_t> timesOfChange = {3, 7};
    initialize(timesOfChange, times, test.methods());
    // define expected 1D output vector (one row is one epoch)
    std::vector<std::vector<size_t>> outExp = {{1, 2, 3, 3, 3, 3, 3, 3, 3},
                                               {0, 0, 0, 1, 2, 3, 4, 4, 4},
                                               {0, 0, 0, 0, 0, 0, 0, 1, 2}};
    // test if result is as expected
    for (size_t i = 0; i < times.size(); i++) {
        for (size_t j = 0; j < timesOfChange.size(); j++) { EXPECT_EQ(epochModel->rho(i, j), outExp.at(j).at(i)); }
    }
}


TEST_F(TEpochModelPoissonTest, _fillRho_TimesOfChange_empty) {
    TBirpCoreBridge test;
    test._readFile("../tests/test_files/method_lureWithCake_TOC_empty.csv");
    auto times = test.timepoints();
    std::vector<size_t> timesOfChange = {};
    initialize(timesOfChange, times, test.methods());

    EXPECT_EQ(timesOfChange.size(), 0);
}


// testing _getEpochStarts
TEST_F(TEpochModelPoissonTest, _getEpochStarts_easy1) {
    TBirpCoreBridge test;
    test._readFile("../tests/test_files/method_loudYelling2.csv");
    auto times = test.timepoints();
    std::vector<size_t> timesOfChange = {3, 8};
    initialize(timesOfChange, times, test.methods());

    // call function
    auto epochStarts = epochModel->_getEpochStarts(times);

    // define expected 1D output vector (rows are columns)
    std::vector<size_t> outExp = {0, 3, 8, 10};

    // test if result is as expected
    EXPECT_EQ(epochStarts.size(), outExp.size());
    for (size_t i = 0; i < timesOfChange.size() + 2; i++) { EXPECT_EQ(epochStarts[i], outExp[i]); }
}

// testing _calculateLLandJP
TEST_F(TEpochModelPoissonTest, _calculateLLandJP_1) {
    TBirpCoreBridge test;
    test._readFile("../tests/test_files/method_diggingInDirt_LLandJP1.csv");
    auto methods = test.methods();
    auto times = test.timepoints();
    std::vector<size_t> timesOfChange = {2};
    std::vector<double> NewGamma = {-2, 1.5};

    initialize(timesOfChange, times, methods, NewGamma);

    // define expected result
    double expectedLL = -52.768731;
    double expectedJP = 4.420733;
    auto [LL, JP] = epochModel->_calculateLLandJP(methods);
    // call function
    EXPECT_FLOAT_EQ(expectedLL, coretools::containerSum(LL));
    EXPECT_FLOAT_EQ(expectedJP, JP);
}


// testing _calculateLLandJP
TEST_F(TEpochModelPoissonTest, _calculateLLandJP_2) {
    TBirpCoreBridge test;
    test._readFile("../tests/test_files/method_diggingInDirt_LLandJP2.csv");
    auto methods = test.methods();
    auto times = test.timepoints();
    std::vector<size_t> timesOfChange = {4};
    std::vector<double> NewGamma = {-1, 1};

    initialize(timesOfChange, times, methods, NewGamma);

    // define expected result
    double expectedLL = -1422.29668;
    double expectedJP = 12.89151;
    auto [LL, JP] = epochModel->_calculateLLandJP(methods);
    // call function
    EXPECT_FLOAT_EQ(expectedLL, coretools::containerSum(LL));
    EXPECT_NEAR(expectedJP, JP, 0.0001); // determinant is not exact in armadillo library
}

//----------------------------------
// TNegBinEpochModel
//----------------------------------

class TNegBinEpochModelBridge : public TNegBinEpochModel {
    // Dummy class inheriting from TEpochModel: allows for testing protected functions
public:
    TNegBinEpochModelBridge(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
                            stattools::TParameterTyped<TypeTime, 1> *TimesOfChange,
                            stattools::TParameterTyped<TypeMu, 2> *Mu, stattools::TParameterTyped<TypeLogA, 1> *LogA)
            : TNegBinEpochModel(Gamma, TimesOfChange, Mu, LogA) {}

    std::pair<std::vector<double>, double> _calculateLLandJP(const std::vector<TMethods> &Methods) const {
        return TNegBinEpochModel::_calculateLLandJP(Methods);
    }

    std::pair<size_t, size_t> _getMethAndLocIndex(size_t mu_index, const std::vector<TMethods> &Methods) {
        return TNegBinEpochModel::_getMethAndLocIndex(mu_index, Methods);
    };


};

class TEpochModelNegBinTest : public testing::Test {
    // class for google test suite: define variables needed for testing TBirpPrior
public:
    // define parameters for test
    std::unique_ptr<stattools::TParameterTyped<TypeGamma, 1>> gamma;
    std::unique_ptr<stattools::TParameterTyped<TypeTime, 1>> timesOfChange;
    std::unique_ptr<stattools::TParameterTyped<TypeMu, 2>> mu;
    std::unique_ptr<stattools::TParameterTyped<TypeLogA, 1>> logA;
    std::shared_ptr<TNegBinEpochModelBridge> epochModel;


    void initialize(const std::vector<size_t> &TimesOfChange, const TUniqueContainer<size_t> &Timepoints,
                    const std::vector<TMethods> &Methods, const TUniqueContainer<std::string> &Locations,
                    std::vector<double> NewGamma = {}, std::vector<double> NewMu = {}, std::vector<double> NewLogA = {},
                    std::vector<double> OldGamma = {}, std::vector<double> OldMu = {},
                    std::vector<double> OldLogA = {}) {
        using namespace stattools::prior;

        const size_t numEpochs = TimesOfChange.size() + 1;
        coretools::instances::parameters().addParameter("timesOfChange",
                                                        coretools::str::concatenateString(TimesOfChange, ","));

        // create a gamma
        if (OldGamma.empty()) { OldGamma.resize(numEpochs, 0.0); }
        if (NewGamma.empty()) { NewGamma.resize(numEpochs, 0.0); }
        testing::initializeStorageParameter(gamma, "gamma", {numEpochs}, OldGamma, NewGamma);

        // create a A
        if (OldLogA.empty()) { OldLogA.resize(1, 0.0); }
        if (NewLogA.empty()) { NewLogA.resize(1, 0.0); }
        testing::initializeStorageParameter(logA, "logA", {1}, OldLogA, NewLogA);

        // create times of change
        testing::initializeStorageParameter(timesOfChange, "timesOfChange", {numEpochs - 1}, TimesOfChange,
                                            TimesOfChange);

        // create mu
        if (OldMu.empty()) { OldMu.resize(Methods.size() * Locations.size(), 1.0); }
        if (NewMu.empty()) { NewMu.resize(Methods.size() * Locations.size(), 1.0); }
        std::array<size_t, 2> dimensionsMu = {Methods.size(), Locations.size()};
        testing::initializeStorageParameter(mu, "mu", dimensionsMu, OldMu, NewMu);

        epochModel = std::make_shared<TNegBinEpochModelBridge>(gamma.get(), timesOfChange.get(),
                                                               mu.get(),
                                                               logA.get()); // get(): unique pointer gets raw pointer // TODO
        epochModel->initializeInferred(Timepoints, Locations, Methods);
    }
};

// testing _calculateLLandJP
TEST_F(TEpochModelNegBinTest, _calculateLL) {
    TBirpCoreBridge test;
    test._readFile("../tests/test_files/method_diggingInDirt_LLandJP_NB.csv");
    auto methods = test.methods();
    auto times = test.timepoints();
    auto locations = test.locations();
    std::vector<size_t> timesOfChange = {3};
    std::vector<double> NewGamma = {-1, 1};
    std::vector<double> NewMu = {1.0};
    std::vector<double> NewLogA = {0}; // a=1 -> log(a)=0
    initialize(timesOfChange, times, methods, locations, NewGamma, NewMu);
    // define expected result
    // double expectedLL = -361.9808; this is from old test -> delete
    double expectedLL = -108.7528;
    auto [LL, JP] = epochModel->_calculateLLandJP(methods);
    // call function
    EXPECT_FLOAT_EQ(expectedLL, coretools::containerSum(LL));
}
// TODO test JP in NB?

// testing _calculateLLandJP
TEST_F(TEpochModelNegBinTest, _getMethAndLocIndex) {
    TBirpCoreBridge test;
    std::vector<std::string> filenames = {"/home/liam/git/birp_cpp2/tests/test_files/method_butterflies1.csv", "/home/liam/git/birp_cpp2/tests/test_files/method_butterflies2.csv",
                                          "/home/liam/git/birp_cpp2/tests/test_files/method_butterflies3.csv"};
    for (const auto& name: filenames) {
        test._readFile(name);
    }
    auto methods = test.methods();
    auto times = test.timepoints();
    auto locations = test.locations();
    std::vector<size_t> timesOfChange = {3}; // is not needed but pero bueno
    initialize(timesOfChange, times, methods, locations);

    std::vector<size_t> method_id_truthVec = {0, 0, 0, 0, 1, 1, 1, 2};
    std::vector<size_t> location_id_truthVec = {0, 1, 2, 3, 0, 1, 2, 0};
    std::size_t counter = 0;
    for (auto m: methods) {
        for (auto loc: m) {
            auto [m_id, l_id] = epochModel->_getMethAndLocIndex(counter, methods);
            // print out stuff
            std::cout << "method: " << m.method_id() << " m_id: " << m_id << " loc: " << loc.location_id() << " l_id: " << l_id << std::endl;
            // call function
            EXPECT_EQ(m_id, method_id_truthVec[counter]);
            EXPECT_EQ(l_id, location_id_truthVec[counter]);
            counter++;
        }
    }
}


