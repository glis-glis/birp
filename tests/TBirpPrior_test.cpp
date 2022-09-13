
#include "../stattools/tests/commonTestingFunctions.h"
#include "TBirpCore.h"
#include "TBirpPrior.h"
#include "gtest/gtest.h"

class TBirpPriorBridge : public TBirpPrior {
	// Dummy class inheriting from TBirpPrior: allows for testing protected functions
public:
	void _readFile(const std::string &Filename){
		return TBirpPrior::_readFile(Filename);
	}
};

TEST(TBirpPriorTest, readFile){
	// construct birp prior
	TBirpPriorBridge bPrior;

	bPrior._readFile("test_counts.txt");

	const auto locations = bPrior.locations();
	EXPECT_EQ(locations[0], "Im Bach");

}




























/*
TEST_F(TBirpPriorTest, _fillRho_easy1) {
	std::vector<size_t> timesOfChange = {2, 5, 9};
	std::vector<size_t> times         = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	initialize(1, times, timesOfChange);

	// call function
	birpPrior->_fillRho(timesOfChange, times);

	// define expected 1D output vector (rows are columns)
	std::vector<std::vector<size_t>> outExp = {{1, 2, 2, 2, 2, 2, 2, 2, 2, 2},
	                                           {0, 0, 1, 2, 3, 3, 3, 3, 3, 3},
	                                           {0, 0, 0, 0, 0, 1, 2, 3, 4, 4},
	                                           {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}};

	// test if result is as expected
	for (size_t i = 0; i < times.size(); i++) {
		for (size_t j = 0; j < timesOfChange.size(); j++) { EXPECT_EQ(birpPrior->rho(i, j), outExp.at(j).at(i)); }
	}
}

TEST_F(TBirpPriorTest, _fillRho_easy2) {
	std::vector<size_t> timesOfChange = {2, 3};
	std::vector<size_t> times         = {1, 2, 3, 4, 5, 6};
	initialize(1, times, timesOfChange);

	// call function
	birpPrior->_fillRho(timesOfChange, times);

	// define expected 1D output vector (rows are columns)
	std::vector<std::vector<size_t>> outExp = {{1, 2, 2, 2, 2, 2}, {0, 0, 1, 1, 1, 1}, {0, 0, 0, 1, 2, 3}};

	// test if result is as expected
	for (size_t i = 0; i < times.size(); i++) {
		for (size_t j = 0; j < timesOfChange.size(); j++) { EXPECT_EQ(birpPrior->rho(i, j), outExp.at(j).at(i)); }
	}
}

TEST_F(TBirpPriorTest, _fillRho_TimesOfChange_same_as_times) {
	std::vector<size_t> timesOfChange = {1, 2};
	std::vector<size_t> times         = {1, 2};
	initialize(1, times, timesOfChange);

	// call function
	birpPrior->_fillRho(timesOfChange, times);

	// define expected 1D output vector (rows are columns)
	std::vector<std::vector<size_t>> outExp = {{1, 1}, {0, 1}, {0, 0}};

	// test if result is as expected
	for (size_t i = 0; i < times.size(); i++) {
		for (size_t j = 0; j < timesOfChange.size(); j++) { EXPECT_EQ(birpPrior->rho(i, j), outExp.at(j).at(i)); }
	}
}

TEST_F(TBirpPriorTest, _fillRho_TimesOfChange_empty) {
	std::vector<size_t> timesOfChange = {};
	std::vector<size_t> times         = {1, 2, 3};
	initialize(1, times, timesOfChange);

	// call function
	birpPrior->_fillRho(timesOfChange, times);

	// define expected 1D output vector (rows are columns)
	std::vector<size_t> outExp = {1, 2, 3};

	// test if result is as expected
	for (size_t i = 0; i < times.size(); i++) { EXPECT_EQ(birpPrior->rho(i, 0), outExp.at(i)); }
}

// testing _getEpochStarts
TEST_F(TBirpPriorTest, _getEpochStarts_easy1) {
	std::vector<size_t> timesOfChange = {3, 8};
	std::vector<size_t> times         = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	initialize(1, times, timesOfChange);

	// call function
	auto result = birpPrior->_getEpochStarts(timesOfChange, times);

	// define expected 1D output vector (rows are columns)
	std::vector<size_t> outExp = {0, 3, 8, 10};

	// test if result is as expected
	for (size_t i = 0; i < timesOfChange.size() + 2; i++) { EXPECT_EQ(result.at(i), outExp.at(i)); }
}

TEST_F(TBirpPriorTest, _getEpochStarts_easy2) {
	std::vector<size_t> timesOfChange = {3, 10};
	std::vector<size_t> times         = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	initialize(1, times, timesOfChange);

	// call function
	auto result = birpPrior->_getEpochStarts(timesOfChange, times);

	// define expected 1D output vector (rows are columns)
	std::vector<size_t> outExp = {0, 3, 10, 11};

	// test if result is as expected
	for (size_t i = 0; i < timesOfChange.size() + 2; i++) { EXPECT_EQ(result.at(i), outExp.at(i)); }
}

TEST_F(TBirpPriorTest, _getEpochStarts_TimesOfChangeEmpty) {
	std::vector<size_t> timesOfChange = {};
	std::vector<size_t> times         = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	initialize(1, times, timesOfChange);

	// call function
	auto result = birpPrior->_getEpochStarts(timesOfChange, times);

	// define expected 1D output vector (rows are columns)
	std::vector<size_t> outExp = {0, 10};

	// test if result is as expected
	for (size_t i = 0; i < timesOfChange.size() + 2; i++) { EXPECT_EQ(result.at(i), outExp.at(i)); }
}

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