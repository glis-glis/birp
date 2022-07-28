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
	double _calculateLL() const {
		return TBirpPrior::_calculateLL();
	}
};

class TBirpPriorTest : public testing::Test {
	// class for google test suite: define variables needed for testing TBirpPrior
public:
	// define parameters for test
	std::unique_ptr<stattools::TParameterTyped<TypeGamma, 1>> gamma;
	std::unique_ptr<stattools::TObservationTyped<TypeEffort, 2>> effort;
	std::unique_ptr<stattools::TObservationTyped<TypeCounts, 2>> counts;

	std::shared_ptr<TBirpPriorBridge> birpPrior;

	void initialize(const std::vector<double> &newGamma={}, const std::vector<double> &oldGamma={}, const std::vector<size_t> &Effort={}, const std::vector<size_t> &Counts={}) {
		// this function is always called whenever a test is called TEST_F(TBirpPriorTest, ...)
		using namespace stattools::prior;

		// create a gamma
		testing::initializeStorageParameter(gamma, "gamma", {oldGamma.size()}, oldGamma, newGamma);

		// create efforts
		auto priorEffort = std::make_shared<TUniformFixed<stattools::TObservationBase, TypeEffort, 2>>();
		testing::initializeStorageObservation(effort, priorEffort, "effort", {2, 10}, Effort);

		// create counts
		birpPrior = std::make_shared<TBirpPriorBridge>(gamma.get(), effort.get(), "birp");
		testing::initializeStorageObservation(counts, birpPrior, "counts", {2, 10}, Counts);
	}


};



// testing _fillRho
TEST_F(TBirpPriorTest, _fillRho_easy1) {
	initialize();
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
	initialize();
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
	initialize();
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
	initialize();
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
	initialize();
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
	initialize();
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
	initialize();
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


// testing _calculateLL
TEST_F(TBirpPriorTest, _calculateLL) {
	// define values for parameters: used to test functions
	std::vector<double> oldGamma     = {0.1, 0.9};
	std::vector<double> newGamma     = {0.2,-0.1};
	std::vector<size_t> valuesEffort = {10,10,10,10,10,10,10,10,10,10};
	std::vector<size_t> valuesCounts = {15,9,17,24,30,26,25,20,26,14};

	initialize(newGamma, oldGamma, valuesEffort, valuesCounts);

	// call function
	EXPECT_FLOAT_EQ(birpPrior->_calculateLL(), -24.04967);

}