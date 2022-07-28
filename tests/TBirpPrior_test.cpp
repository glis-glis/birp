#include "../stattools/tests/commonTestingFunctions.h"
#include "TBirpCore.h"
#include "TBirpPrior.h"
#include "gtest/gtest.h"

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
};

class TBirpPriorTest : public testing::Test {
	// class for google test suite: define variables needed for testing TBirpPrior
public:
	// define values for parameters: used to test functions
	std::vector<double> oldGamma     = {0.1, 0.9, 0.3, 0.7};
	std::vector<double> newGamma     = {0.2, 0.8, 0.4, 0.6};
	std::vector<double> valuesEffort = {0.09, 0.7,  0.21, 1.85, 0.14, 0.04, 1.27, 0.31, 0.83, 0,
	                                    1.31, 0.53, 0.19, 2.9,  2.48, 0.8,  0.16, 0.42, 0.65, 3.31};
	std::vector<size_t> valuesCounts = {26, 19, 6, 70, 83, 58, 96, 11, 0, 49, 26, 45, 96, 41, 67, 36, 11, 67, 43, 32};

	// define parameters for test
	std::unique_ptr<stattools::TParameterTyped<TypeGamma, 1>> gamma;
	std::unique_ptr<stattools::TObservationTyped<TypeEffort, 2>> effort;
	std::unique_ptr<stattools::TObservationTyped<TypeCounts, 2>> counts;

	std::shared_ptr<TBirpPriorBridge> birpPrior;

	void SetUp() override {
		// this function is always called whenever a test is called TEST_F(TBirpPriorTest, ...)
		using namespace stattools::prior;

		// create a gamma
		testing::initializeStorageParameter(gamma, "gamma", {oldGamma.size()}, oldGamma, newGamma);

		// create efforts
		auto priorEffort = std::make_shared<TUniformFixed<stattools::TObservationBase, TypeEffort, 2>>();
		testing::initializeStorageObservation(effort, priorEffort, "effort", {2, 10}, valuesEffort);

		// create counts
		birpPrior = std::make_shared<TBirpPriorBridge>(gamma.get(), effort.get(), "birp");
		testing::initializeStorageObservation(counts, birpPrior, "counts", {2, 10}, valuesCounts);
	}
};

// testing _fillRho
TEST_F(TBirpPriorTest, _fillRho_easy1) {
	// vector definieren mit times of change und dann vector mit times. --> werte geben (Sinnvoll). also test border
	// cases (empty vector), only one vector, ...
	std::vector<size_t> timesOfChange = {2, 5, 9};
	std::vector<size_t> times         = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
	coretools::instances::parameters().addParameter("timesOfChange",
	                                                coretools::str::concatenateString(timesOfChange, ","));

	// call function
	birpPrior->_fillRho(timesOfChange, times);


    // define expected 1D output vector (rows are columns)
    std::vector<std::vector<size_t>> outExp = {
            {1, 2, 2, 2, 2, 2, 2, 2, 2, 2},
            {0, 0, 1, 2, 3, 3, 3, 3, 3, 3},
            {0, 0, 0, 0, 0, 1, 2, 3, 4, 4},
            {0, 0, 0, 0, 0, 0, 0, 0, 0, 1}
    };

	// test if result is as expected
	for (size_t i = 0; i < times.size(); i++) {
		for (size_t j = 0; j < timesOfChange.size(); j++) {
			EXPECT_EQ(birpPrior->rho(i, j), outExp.at(j).at(i));
		}
	}

}

TEST_F(TBirpPriorTest, _fillRho_easy2) {
    std::vector<size_t> timesOfChange = {2, 3};
    std::vector<size_t> times         = {1, 2, 3, 4, 5, 6};
    coretools::instances::parameters().addParameter("timesOfChange",
                                                    coretools::str::concatenateString(timesOfChange, ","));

    // call function
    birpPrior->_fillRho(timesOfChange, times);

    // define expected 1D output vector (rows are columns)
    std::vector<std::vector<size_t>> outExp = {
            {1, 2, 2, 2, 2, 2},
            {0, 0, 1, 1, 1, 1},
            {0, 0, 0, 1, 2, 3}
    };

    // test if result is as expected
    for (size_t i = 0; i < times.size(); i++) {
        for (size_t j = 0; j < timesOfChange.size(); j++) {
            EXPECT_EQ(birpPrior->rho(i, j), outExp.at(j).at(i));
        }
    }

}

TEST_F(TBirpPriorTest, _fillRho_TimesOfChange_same_as_times) {
    std::vector<size_t> timesOfChange = {1,2};
    std::vector<size_t> times         = {1,2};
    coretools::instances::parameters().addParameter("timesOfChange",
                                                    coretools::str::concatenateString(timesOfChange, ","));

    // call function
    birpPrior->_fillRho(timesOfChange, times);

    // define expected 1D output vector (rows are columns)
    std::vector<std::vector<size_t>> outExp = {
            {1, 1},
            {0, 1},
            {0, 0}
    };

    // test if result is as expected
    for (size_t i = 0; i < times.size(); i++) {
        for (size_t j = 0; j < timesOfChange.size(); j++) {
            EXPECT_EQ(birpPrior->rho(i, j), outExp.at(j).at(i));
        }
    }

}

TEST_F(TBirpPriorTest, _fillRho_TimesOfChange_empty) {
    std::vector<size_t> timesOfChange = {};
    std::vector<size_t> times         = {1,2,3};
    coretools::instances::parameters().addParameter("timesOfChange",
                                                    coretools::str::concatenateString(timesOfChange, ","));

	// call function
	birpPrior->_fillRho(timesOfChange, times);

	// define expected 1D output vector (rows are columns)
	std::vector<size_t> outExp = {1,2,3};

	// test if result is as expected
	for (size_t i = 0; i < times.size(); i++) {
			EXPECT_EQ(birpPrior->rho(i,0), outExp.at(i));
	}
}

// testing _getEpochStarts
TEST_F(TBirpPriorTest, _getEpochStarts_easy1) {
    std::vector<size_t> timesOfChange = {3,8};
    std::vector<size_t> times         = {1,2,3,4,5,6,7,8,9,10};
    coretools::instances::parameters().addParameter("timesOfChange",
                                                    coretools::str::concatenateString(timesOfChange, ","));
    // call function
    auto result = birpPrior->_getEpochStarts(timesOfChange, times);

    // define expected 1D output vector (rows are columns)
    std::vector<size_t> outExp = {0,3,8,10};

    // test if result is as expected
    for (size_t i = 0; i < timesOfChange.size() + 2; i++) {
            EXPECT_EQ(result.at(i), outExp.at(i));
    }

}

TEST_F(TBirpPriorTest, _getEpochStarts_easy2) {
    std::vector<size_t> timesOfChange = {3,10};
    std::vector<size_t> times         = {1,2,3,4,5,6,7,8,9,10};
    coretools::instances::parameters().addParameter("timesOfChange",
                                                    coretools::str::concatenateString(timesOfChange, ","));
    // call function
	auto result = birpPrior->_getEpochStarts(timesOfChange, times);

    // define expected 1D output vector (rows are columns)
    std::vector<size_t> outExp = {0,3,10,11};

    // test if result is as expected
    for (size_t i = 0; i < timesOfChange.size() + 2; i++) {
        EXPECT_EQ(result.at(i), outExp.at(i));
    }

}

TEST_F(TBirpPriorTest, _getEpochStarts_TimesOfChangeEmpty) {
	std::vector<size_t> timesOfChange = {};
	std::vector<size_t> times         = {1,2,3,4,5,6,7,8,9,10};
	coretools::instances::parameters().addParameter("timesOfChange",
	                                                coretools::str::concatenateString(timesOfChange, ","));
	// call function
	auto result = birpPrior->_getEpochStarts(timesOfChange, times);

	// define expected 1D output vector (rows are columns)
	std::vector<size_t> outExp = {0,10};

	// test if result is as expected
	for (size_t i = 0; i < timesOfChange.size() + 2; i++) {
		EXPECT_EQ(result.at(i), outExp.at(i));
	}

}