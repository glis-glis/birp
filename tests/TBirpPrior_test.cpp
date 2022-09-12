
#include "../stattools/tests/commonTestingFunctions.h"
#include "TBirpCore.h"
#include "TBirpPrior.h"
#include "gtest/gtest.h"

/*

class TBirpPriorBridge : public TBirpPrior {
	// Dummy class inheriting from TBirpPrior: allows for testing protected functions
public:
	TBirpPriorBridge(stattools::TParameterTyped<TypeGamma, 1> *Gamma,
	                 stattools::TObservationTyped<TypeEffort, 2> *Effort, const std::string &Filename)
	    : TBirpPrior(Gamma, Effort, Filename){};

	// make protected functions public
	void _fillRho(const std::vector<size_t> &TimesOfChange, const std::vector<size_t> &Times) {
		TBirpPrior::_fillRho(TimesOfChange, Times);
	}
	std::vector<size_t> _getEpochStarts(const std::vector<size_t> &TimesOfChange, const std::vector<size_t> &Times) {
		return TBirpPrior::_getEpochStarts(TimesOfChange, Times);
	}
	double _calculateLL() const { return TBirpPrior::_calculateLL(); }
};

class TBirpPriorTest : public testing::Test {
	// class for google test suite: define variables needed for testing TBirpPrior
public:
	// define parameters for test
	std::unique_ptr<stattools::TParameterTyped<TypeGamma, 1>> gamma;
	std::unique_ptr<stattools::TObservationTyped<TypeEffort, 2>> effort;
	std::unique_ptr<stattools::TObservationTyped<TypeCounts, 2>> counts;

	std::shared_ptr<TBirpPriorBridge> birpPrior;

	void initialize(size_t NumMethods, const std::vector<size_t> &TimePoints, const std::vector<size_t> &TimesOfChange,
	                std::vector<double> NewGamma = {}, std::vector<double> OldGamma = {},
	                std::vector<size_t> ValuesEffort = {}, std::vector<size_t> ValuesCount = {}) {
		// this function is always called whenever a test is called TEST_F(TBirpPriorTest, ...)
		using namespace stattools::prior;

		const size_t numTimePoints = TimePoints.size();
		const size_t numEpochs     = TimesOfChange.size() + 1;

		// add times of change to command line
		coretools::instances::parameters().addParameter("timesOfChange",
		                                                coretools::str::concatenateString(TimesOfChange, ","));

		// create a gamma
		if (OldGamma.empty()) { OldGamma.resize(numEpochs, 0.0); }
		if (NewGamma.empty()) { NewGamma.resize(numEpochs, 0.0); }
		testing::initializeStorageParameter(gamma, "gamma", {numEpochs}, OldGamma, NewGamma);

		// create dimension names
		auto timeNames = std::make_shared<coretools::TNamesStrings>();
		for (size_t t = 0; t < numTimePoints; t++) { timeNames->addName({coretools::str::toString(TimePoints[t])}); }
		auto methodNames = std::make_shared<coretools::TNamesIndices>(NumMethods);

		// create efforts
		if (ValuesEffort.empty()) { ValuesEffort.resize(NumMethods * numTimePoints); }
		auto priorEffort = std::make_shared<TUniformFixed<stattools::TObservationBase, TypeEffort, 2>>();
		testing::initializeStorageObservation(effort, priorEffort, "effort", {NumMethods, numTimePoints}, ValuesEffort,
		                                      {methodNames, timeNames});

		// create counts
		if (ValuesCount.empty()) { ValuesCount.resize(NumMethods * numTimePoints); }
		birpPrior = std::make_shared<TBirpPriorBridge>(gamma.get(), effort.get(), "birp");
		testing::initializeStorageObservation(counts, birpPrior, "counts", {NumMethods, numTimePoints}, ValuesCount,
		                                      {methodNames, timeNames});
	}
};

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
}*/

TEST(Dummy,dummy){
	TUniqueContainer<std::string> _locations;
	_locations.add("im Bach");
	_locations.add("Ufem Bärg");
	_locations.add("ir Höhli");
	std::cout << _locations.size() << std::endl;;
	_locations.add("im Bach");
	std::cout << _locations.size() << std::endl;
	std::cout << _locations.getIndex("im Bach") << std::endl;
}