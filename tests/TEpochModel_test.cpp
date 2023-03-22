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

	std::pair<std::vector<double>, double> _calculateLLandJP(const TMethods &Data) const {
		return TPoissonEpochModel::_calculateLLandJP(Data);
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
		epochModel->initializeInferred(Timepoints, {}, {});
	}
};

TEST_F(TEpochModelPoissonTest, _fillRho_easy1) {
	TBirpCoreBridge test;
	test._readFile("../tests/test_files/test_epochModel1.csv");
	test._sortTimepoints();
	auto times                        = test.timepoints();
	std::vector<size_t> timesOfChange = {2013, 2018};
	initialize(timesOfChange, times);

	// define expected 1D output vector (one row is one epoch)
	std::vector<std::vector<size_t>> outExp = {
	    {2010, 2012, 2013, 2013, 2013, 2013, 2013}, {0, 0, 1, 2, 4, 5, 5}, {0, 0, 0, 0, 0, 2, 11}};
	// test if result is as expected
	for (size_t i = 0; i < times.size(); i++) {
		for (size_t j = 0; j < timesOfChange.size() + 1; j++) { EXPECT_EQ(epochModel->rho(i, j), outExp.at(j).at(i)); }
	}
}

TEST_F(TEpochModelPoissonTest, _fillRho_easy2) {
	TBirpCoreBridge test;
	test._readFile("../tests/test_files/test_epochModel2.csv");
	auto times                        = test.timepoints();
	std::vector<size_t> timesOfChange = {2};
	initialize(timesOfChange, times);

	// define expected 1D output vector (one row is one epoch)
	std::vector<std::vector<size_t>> outExp = {{1, 2, 2, 2, 2}, {0, 0, 1, 2, 3}};
	// test if result is as expected
	for (size_t i = 0; i < times.size(); i++) {
		for (size_t j = 0; j < timesOfChange.size() + 1; j++) { EXPECT_EQ(epochModel->rho(i, j), outExp.at(j).at(i)); }
	}
}

TEST_F(TEpochModelPoissonTest, _fillRho_TimesOfChange_same_as_times) {
	TBirpCoreBridge test;
	test._readFile("../tests/test_files/test_epochModel2_TOC_same_as_times.csv");
	auto times                        = test.timepoints();
	std::vector<size_t> timesOfChange = {1, 2};
	initialize(timesOfChange, times);

	// define expected 1D output vector (one row is one epoch)
	std::vector<std::vector<size_t>> outExp = {{1, 1}, {0, 1}, {0, 0}};
	// test if result is as expected
	for (size_t i = 0; i < times.size(); i++) {
		for (size_t j = 0; j < timesOfChange.size(); j++) { EXPECT_EQ(epochModel->rho(i, j), outExp.at(j).at(i)); }
	}
}

TEST_F(TEpochModelPoissonTest, _fillRho_TimesOfChange_empty) {
	TBirpCoreBridge test;
	test._readFile("../tests/test_files/test_epochModel2_TOC_empty.csv");
	auto times                        = test.timepoints();
	std::vector<size_t> timesOfChange = {};
	initialize(timesOfChange, times);

	EXPECT_EQ(timesOfChange.size(), 0);
}

// testing _getEpochStarts
TEST_F(TEpochModelPoissonTest, _getEpochStarts_easy1) {
	TBirpCoreBridge test;
	test._readFile("../tests/test_files/test_epochStarts1.csv");
	auto times                        = test.timepoints();
	std::vector<size_t> timesOfChange = {3, 8};
	initialize(timesOfChange, times);

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
	test._readFile("../tests/test_files/test_calcLLandJP1.csv");
	auto data                         = test.methods();
	auto times                        = test.timepoints();
	std::vector<size_t> timesOfChange = {2};
	std::vector<double> NewGamma      = {-2, 1.5};

	initialize(timesOfChange, times, NewGamma);

	// define expected result
	double expectedLL = -52.768731;
	double expectedJP = 4.420733;
	auto [LL, JP]     = epochModel->_calculateLLandJP(data[0]);
	// call function
	EXPECT_FLOAT_EQ(expectedLL, coretools::containerSum(LL));
	EXPECT_FLOAT_EQ(expectedJP, JP);
}

// testing _calculateLLandJP
TEST_F(TEpochModelPoissonTest, _calculateLLandJP_2) {
	TBirpCoreBridge test;
	test._readFile("../tests/test_files/test_calcLLandJP2.csv");
	auto data                         = test.methods();
	auto times                        = test.timepoints();
	std::vector<size_t> timesOfChange = {4};
	std::vector<double> NewGamma      = {-1, 1};

	initialize(timesOfChange, times, NewGamma);

	// define expected result
	double expectedLL = -1422.29668;
	double expectedJP = 12.89151;
	auto [LL, JP]     = epochModel->_calculateLLandJP(data[0]);
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

	std::pair<std::vector<double>, double> _calculateLLandJP(const TMethods &Data) const {
		return TNegBinEpochModel::_calculateLLandJP(Data);
	}
};

class TEpochModelNegBinTest : public testing::Test {
	// class for google test suite: define variables needed for testing TBirpPrior
public:
	// define parameters for test
	std::unique_ptr<stattools::TParameterTyped<TypeGamma, 1>> gamma;
	std::unique_ptr<stattools::TParameterTyped<TypeTime, 1>> timesOfChange;
	std::unique_ptr<stattools::TParameterTyped<TypeMu, 2>> mu;
    std::unique_ptr<stattools::TParameterTyped<TypeLogA, 1>> logA; // TODO
	std::shared_ptr<TNegBinEpochModelBridge> epochModel;

	void initialize(const std::vector<size_t> &TimesOfChange, const TUniqueContainer<size_t> &Timepoints,
	                const std::vector<TMethods> &Methods, const TUniqueContainer<std::string> &Locations,
	                std::vector<double> NewGamma = {}, std::vector<double> NewMu = {},
	                std::vector<double> OldGamma = {}, std::vector<double> OldMu = {}) {
		using namespace stattools::prior;

		const size_t numEpochs = TimesOfChange.size() + 1;
		coretools::instances::parameters().addParameter("timesOfChange",
		                                                coretools::str::concatenateString(TimesOfChange, ","));

		// create a gamma
		if (OldGamma.empty()) { OldGamma.resize(numEpochs, 0.0); }
		if (NewGamma.empty()) { NewGamma.resize(numEpochs, 0.0); }
		testing::initializeStorageParameter(gamma, "gamma", {numEpochs}, OldGamma, NewGamma);

		// create times of change
		testing::initializeStorageParameter(timesOfChange, "timesOfChange", {numEpochs - 1}, TimesOfChange,
		                                    TimesOfChange);

		// create mu
		if (OldMu.empty()) { OldMu.resize(Methods.size()*Locations.size(), 1.0); }
		if (NewMu.empty()) { NewMu.resize(Methods.size()*Locations.size(), 1.0); }
		std::array<size_t, 2> dimensionsMu = {Methods.size(), Locations.size()};
		testing::initializeStorageParameter(mu, "mu", dimensionsMu, OldMu, NewMu);

		epochModel = std::make_shared<TNegBinEpochModelBridge>(gamma.get(), timesOfChange.get(),
		                                                       mu.get(), logA.get()); // get(): unique pointer gets raw pointer // TODO
		epochModel->initializeInferred(Timepoints, Locations, Methods);
	}
};

// testing _calculateLLandJP
TEST_F(TEpochModelNegBinTest, _calculateLLandJP) {
	TBirpCoreBridge test;
	test._readFile("../tests/test_files/NB_dirichletMultinomial_LL_test1a.csv");
	auto data                         = test.methods();
	auto times                        = test.timepoints();
	auto methods                      = test.methods();
	auto locations                    = test.locations();
	std::vector<size_t> timesOfChange = {3};
	std::vector<double> NewGamma      = {-1, 1};
	std::vector<double> NewMu         = {3,1,7,8};
	initialize(timesOfChange, times, methods, locations, NewGamma, NewMu);
	// define expected result
	double expectedLL = -361.9808;
	auto [LL, JP]     = epochModel->_calculateLLandJP(data[0]);
	// call function
	EXPECT_FLOAT_EQ(expectedLL, coretools::containerSum(LL));
}