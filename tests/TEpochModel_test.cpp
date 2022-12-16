//
// Created by liam on 12/12/22.
//

#include "../stattools/tests/commonTestingFunctions.h"
#include "TBirpPrior.h"
#include "TEpochModel.h"
#include "gtest/gtest.h"


class TBirpPriorBridge : public TBirpPrior {
	// Dummy class inheriting from TBirpPrior: allows for testing protected functions
public:
	void _readFile(const std::string &Filename){
		return TBirpPrior::_readFile(Filename);
	}
	auto gettimes(){
		return _timepoints;
	}
	void sortTimepoints(){
		return TBirpPrior::_sortTimepoints();
	}

};

class TEpochModelBridge : public TEpochModel {
	// Dummy class inheriting from TEpochModel: allows for testing protected functions
public:
	TEpochModelBridge(stattools::TParameterTyped<TypeGamma, 1> *Gamma, stattools::TParameterTyped<TypeTime, 1>* TimesOfChange, const std::string &Filename): TEpochModel(Gamma, TimesOfChange, Filename){}

	void _fillRho(const TUniqueContainer<size_t> &Timepoints) {
		return TEpochModel::_fillRho(Timepoints);
	}

    std::vector<size_t> _getEpochStarts(const TUniqueContainer<size_t> &Timepoints) {
        return TEpochModel::_getEpochStarts(Timepoints);
    }

};

class TEpochModelTest : public testing::Test {
	// class for google test suite: define variables needed for testing TBirpPrior
public:
	// define parameters for test
	std::unique_ptr<stattools::TParameterTyped<TypeGamma, 1>> gamma;
	std::shared_ptr<TEpochModelBridge> epochModel;
	std::unique_ptr<stattools::TParameterTyped<TypeTime,1>> timesOfChange;

	void initialize(const std::vector<size_t> &TimesOfChange,
	                std::vector<double> NewGamma = {}, std::vector<double> OldGamma = {}) {
		// this function is always called whenever a test is called TEST_F(TBirpPriorTest, ...)
		using namespace stattools::prior;

		const size_t numEpochs     = TimesOfChange.size() + 1;

		// create a gamma
		if (OldGamma.empty()) { OldGamma.resize(numEpochs, 0.0); }
		if (NewGamma.empty()) { NewGamma.resize(numEpochs, 0.0); }
		testing::initializeStorageParameter(gamma, "gamma", {numEpochs}, OldGamma, NewGamma);
		testing::initializeStorageParameter(timesOfChange, "timesOfChange", {numEpochs-1}, TimesOfChange, TimesOfChange);

		epochModel = std::make_shared<TEpochModelBridge>(gamma.get(), timesOfChange.get(), "out");  // get(): unique pointer gets raw pointer

		epochModel->initialize();
	}
};



TEST_F(TEpochModelTest, _fillRho_easy1) {
	TBirpPriorBridge test;
	test._readFile("../tests/test_files/test_epochModel1.csv");
	test.sortTimepoints();
	auto times = test.gettimes();
	std::vector<size_t> timesOfChange = {2013, 2018};
	initialize(timesOfChange);
	// call function
	epochModel->_fillRho(times);

	// define expected 1D output vector (one row is one epoch)
	std::vector<std::vector<size_t>> outExp = {{2010,2012,2013,2013,2013,2013,2013},
	                                           {0,0,1,2,4,5,5},
	                                           {0,0,0,0,0,2,11}};
	// test if result is as expected
	for (size_t i = 0; i < times.size(); i++) {
		for (size_t j = 0; j < timesOfChange.size()+1; j++) {
			EXPECT_EQ(epochModel->rho(i, j), outExp.at(j).at(i)); }
	}
}


TEST_F(TEpochModelTest, _fillRho_easy2) {
	TBirpPriorBridge test;
	test._readFile("../tests/test_files/test_epochModel2.csv");
	auto times = test.gettimes();
	std::vector<size_t> timesOfChange = {2};
	initialize(timesOfChange);
	// call function
	epochModel->_fillRho(times);

	// define expected 1D output vector (one row is one epoch)
	std::vector<std::vector<size_t>> outExp = {{1,2,2,2,2},
	                                           {0,0,1,2,3}};
	// test if result is as expected
	for (size_t i = 0; i < times.size(); i++) {
		for (size_t j = 0; j < timesOfChange.size()+1; j++) {
			EXPECT_EQ(epochModel->rho(i, j), outExp.at(j).at(i));
		}
	}
}

TEST_F(TEpochModelTest, _fillRho_TimesOfChange_same_as_times) {
    TBirpPriorBridge test;
    test._readFile("../tests/test_files/test_epochModel2_TOC_same_as_times.csv");
    auto times = test.gettimes();
    std::vector<size_t> timesOfChange = {1,2};
    initialize(timesOfChange);
    // call function
    epochModel->_fillRho(times);

    // define expected 1D output vector (one row is one epoch)
    std::vector<std::vector<size_t>> outExp = {{1,1},
                                               {0,1},
                                               {0,0}};
    // test if result is as expected
    for (size_t i = 0; i < times.size(); i++) {
        for (size_t j = 0; j < timesOfChange.size(); j++) { EXPECT_EQ(epochModel->rho(i, j), outExp.at(j).at(i)); }
    }
}


TEST_F(TEpochModelTest, _fillRho_TimesOfChange_empty) {
    TBirpPriorBridge test;
    test._readFile("../tests/test_files/test_epochModel2_TOC_empty.csv");
    auto times = test.gettimes();
    std::vector<size_t> timesOfChange = {};
    initialize(timesOfChange);
    epochModel->_fillRho(times);
    EXPECT_EQ(timesOfChange.size(), 0);
}


// testing _getEpochStarts
TEST_F(TEpochModelTest, _getEpochStarts_easy1) {
    TBirpPriorBridge test;
    test._readFile("../tests/test_files/test_epochStarts1.csv");
    auto times = test.gettimes();
    std::vector<size_t> timesOfChange = {3,8};
    initialize(timesOfChange);
    // call function
    auto epochStarts = epochModel->_getEpochStarts(times);

	// define expected 1D output vector (rows are columns)
	std::vector<size_t> outExp = {0, 3, 8, 10};

	// test if result is as expected
    EXPECT_EQ(epochStarts.size(), outExp.size());
	for (size_t i = 0; i < timesOfChange.size() + 2; i++) { EXPECT_EQ(epochStarts[i], outExp[i]); }
}


/*

// testing _calculateLL
TEST_F(TBirpPriorTest, _calculateLL_oneDataset) {
    // define values for parameters: used to test functions
    std::vector<double> oldGamma      = {0.1, 0.9};
    std::vector<double> newGamma      = {0.2, -0.1};
    std::vector<size_t> valuesEffort  = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
    std::vector<size_t> valuesCounts  = {15, 9, 17, 24, 30, 26, 25, 20, 26, 14};
    std::vector<size_t> timesOfChange = {5};
    std::vector<size_t> times         = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    initialize(1, times, timesOfChange, newGamma, oldGamma, valuesEffort, valuesCounts);

    // call function
    EXPECT_FLOAT_EQ(birpPrior->_calculateLL(), -467.4725); // not the same as in R due to normalization constant not calculated here
}


// testing _calculateLL
TEST_F(TBirpPriorTest, _calculateLL_twoDatasets) {
    // define values for parameters: used to test functions
    std::vector<double> oldGamma      = {0.1, 0.9};
    std::vector<double> newGamma      = {0.8, -0.3};
    std::vector<size_t> valuesEffort  = {1,2,3,4,5,6,7,8,9,10, 11, 12, 13, 14, 15, 16, 17, 18, 19,20};
    std::vector<size_t> valuesCounts  = {1, 9, 32, 95, 259, 236, 217, 179, 154, 137, 14, 66, 155, 349, 829, 675, 541, 335, 315, 265};
    std::vector<size_t> timesOfChange = {5};
    std::vector<size_t> times         = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};

    initialize(2, times, timesOfChange, newGamma, oldGamma, valuesEffort, valuesCounts);

    // call function
    EXPECT_FLOAT_EQ(birpPrior->_calculateLL(), -9867.187); // not the same as in R due to normalization constant not calculated here
}
*/