//
// Created by madleina on 22.07.22.
//

#include "TData.h"
#include "gtest/gtest.h"

class TDataTest : public testing::Test {
public:
	coretools::TMultiDimensionalStorage<TypeCounts, 2> countsWriting;
	coretools::TMultiDimensionalStorage<TypeCounts, 2> countsReading;

	coretools::TMultiDimensionalStorage<TypeEffort, 2> effortsWriting;
	coretools::TMultiDimensionalStorage<TypeEffort, 2> effortsReading;

	void SetUp() override {
		coretools::instances::parameters().addParameter("counts", "birp_test_counts.txt");
		coretools::instances::parameters().addParameter("efforts", "birp_test_efforts.txt");
	}
};

TEST_F(TDataTest, singleMethod_singleTime) {
	// first write...
	countsWriting.resize(std::array<size_t, 2>({1, 1}));
	effortsWriting.resize(std::array<size_t, 2>({1, 1}));
	countsWriting(0, 0)  = 10;
	effortsWriting(0, 0) = 0.5;

	TData::write("birp_test", countsWriting, effortsWriting);

	// should throw: not sufficient data
	EXPECT_ANY_THROW(TData::read(countsReading, effortsReading));
}

TEST_F(TDataTest, singleMethod_twoTimes_zeroEffort_nonZeroCounts) {
	// first write...
	countsWriting.resize(std::array<size_t, 2>({1, 2}));
	effortsWriting.resize(std::array<size_t, 2>({1, 2}));

	countsWriting(0, 0)  = 10; // should also be zero since efforts is zero
	countsWriting(0, 1)  = 20;
	effortsWriting(0, 0) = 0.0;
	effortsWriting(0, 1) = 1.0;

	TData::write("birp_test", countsWriting, effortsWriting);

	// should throw: not sufficient data
	EXPECT_ANY_THROW(TData::read(countsReading, effortsReading));
}

TEST_F(TDataTest, singleMethod_twoTimes_oneHasMissingData) {
	// first write...
	countsWriting.resize(std::array<size_t, 2>({1, 2}));
	effortsWriting.resize(std::array<size_t, 2>({1, 2}));

	countsWriting(0, 0)  = 0;
	countsWriting(0, 1)  = 20;
	effortsWriting(0, 0) = 0.0;
	effortsWriting(0, 1) = 1.0;

	TData::write("birp_test", countsWriting, effortsWriting);

	// should throw: not sufficient data
	EXPECT_ANY_THROW(TData::read(countsReading, effortsReading));
}

TEST_F(TDataTest, differentNumMethods) {
	// first write...
	countsWriting.resize(std::array<size_t, 2>({1, 2}));
	effortsWriting.resize(std::array<size_t, 2>({2, 2}));

	countsWriting(0, 0)  = 10;
	countsWriting(0, 1)  = 20;
	effortsWriting(0, 0) = 0.5;
	effortsWriting(0, 1) = 1.0;
	effortsWriting(1, 0) = 1.5;
	effortsWriting(1, 1) = 2.0;

	TData::write("birp_test", countsWriting, effortsWriting);

	// should throw: different number of methods
	EXPECT_ANY_THROW(TData::read(countsReading, effortsReading));
}

TEST_F(TDataTest, differentNumTimePoints) {
	// first write...
	countsWriting.resize(std::array<size_t, 2>({1, 2}));
	effortsWriting.resize(std::array<size_t, 2>({1, 3}));

	countsWriting(0, 0)  = 10;
	countsWriting(0, 1)  = 20;
	effortsWriting(0, 0) = 0.5;
	effortsWriting(0, 1) = 1.0;
	effortsWriting(0, 2) = 1.5;

	TData::write("birp_test", countsWriting, effortsWriting);

	// should throw: different number of time points
	EXPECT_ANY_THROW(TData::read(countsReading, effortsReading));
}

TEST_F(TDataTest, methodNamesDontMatch) {
	// first write...
	countsWriting.resize(std::array<size_t, 2>({1, 2}));
	effortsWriting.resize(std::array<size_t, 2>({1, 2}));

	countsWriting.setDimensionName(std::make_shared<coretools::TNamesStrings>(), 0);
	(*countsWriting.getDimensionName(0)).addName({"method1"});

	countsWriting(0, 0)  = 10;
	countsWriting(0, 1)  = 20;
	effortsWriting(0, 0) = 0.5;
	effortsWriting(0, 1) = 1.0;

	TData::write("birp_test", countsWriting, effortsWriting);

	// should throw: different method names
	EXPECT_ANY_THROW(TData::read(countsReading, effortsReading));
}

TEST_F(TDataTest, timeNamesDontMatch) {
	// first write...
	countsWriting.resize(std::array<size_t, 2>({1, 2}));
	effortsWriting.resize(std::array<size_t, 2>({1, 2}));

	countsWriting.setDimensionName(std::make_shared<coretools::TNamesStrings>(), 1);
	(*countsWriting.getDimensionName(1)).addName({"1"});
	(*countsWriting.getDimensionName(1)).addName({"100"});

	countsWriting(0, 0)  = 10;
	countsWriting(0, 1)  = 20;
	effortsWriting(0, 0) = 0.5;
	effortsWriting(0, 1) = 1.0;

	TData::write("birp_test", countsWriting, effortsWriting);

	// should throw: different time point names
	EXPECT_ANY_THROW(TData::read(countsReading, effortsReading));
}

TEST_F(TDataTest, singleMethod_twoTimes) {
	// first write...
	countsWriting.resize(std::array<size_t, 2>({1, 2}));
	effortsWriting.resize(std::array<size_t, 2>({1, 2}));

	countsWriting(0, 0)  = 10;
	countsWriting(0, 1)  = 20;
	effortsWriting(0, 0) = 0.5;
	effortsWriting(0, 1) = 1.0;

	TData::write("birp_test", countsWriting, effortsWriting);

	// should work
	TData::read(countsReading, effortsReading);
	EXPECT_EQ(countsWriting.dimensions()[0], countsReading.dimensions()[0]);
	EXPECT_EQ(countsWriting.dimensions()[1], countsReading.dimensions()[1]);
	EXPECT_EQ(effortsWriting.dimensions()[0], effortsReading.dimensions()[0]);
	EXPECT_EQ(effortsWriting.dimensions()[1], effortsReading.dimensions()[1]);

	for (size_t i = 0; i < countsWriting.dimensions()[0]; i++) {
		for (size_t k = 0; k < countsWriting.dimensions()[1]; k++) {
			EXPECT_EQ(countsWriting(i, k), countsReading(i, k));
			EXPECT_EQ(effortsWriting(i, k), effortsReading(i, k));
		}
	}
}

TEST_F(TDataTest, threeMethods_threeTimes) {
	// first write...
	countsWriting.resize(std::array<size_t, 2>({3, 3}));
	effortsWriting.resize(std::array<size_t, 2>({3, 3}));

	for (size_t i = 0; i < countsWriting.dimensions()[0]; i++) {
		for (size_t k = 0; k < countsWriting.dimensions()[1]; k++) {
			countsWriting(i, k)  = i + k;
			effortsWriting(i, k) = (double)(i + k) / 10.;
		}
	}
	TData::write("birp_test", countsWriting, effortsWriting);

	// should work
	TData::read(countsReading, effortsReading);
	EXPECT_EQ(countsWriting.dimensions()[0], countsReading.dimensions()[0]);
	EXPECT_EQ(countsWriting.dimensions()[1], countsReading.dimensions()[1]);
	EXPECT_EQ(effortsWriting.dimensions()[0], effortsReading.dimensions()[0]);
	EXPECT_EQ(effortsWriting.dimensions()[1], effortsReading.dimensions()[1]);

	for (size_t i = 0; i < countsWriting.dimensions()[0]; i++) {
		for (size_t k = 0; k < countsWriting.dimensions()[1]; k++) {
			EXPECT_EQ(countsWriting(i, k), countsReading(i, k));
			EXPECT_EQ(effortsWriting(i, k), effortsReading(i, k));
		}
	}
}

TEST_F(TDataTest, threeMethods_threeTimes_filterOutMethod) {
	// first write...
	countsWriting.resize(std::array<size_t, 2>({3, 3}));
	effortsWriting.resize(std::array<size_t, 2>({3, 3}));

	// second method: missing data
	for (size_t i = 0; i < countsWriting.dimensions()[0]; i++) {
		for (size_t k = 0; k < countsWriting.dimensions()[1]; k++) {
			if (i == 1 && k > 0) {
				countsWriting(i, k)  = 0;
				effortsWriting(i, k) = 0.0;
			} else {
				countsWriting(i, k)  = i + k;
				effortsWriting(i, k) = (double)(i + k) / 10.;
			}
		}
	}
	TData::write("birp_test", countsWriting, effortsWriting);

	// should work, but filter out method 2!
	TData::read(countsReading, effortsReading);
	EXPECT_EQ(countsReading.dimensions()[0], 2);
	EXPECT_EQ(countsWriting.dimensions()[1], countsReading.dimensions()[1]);
	EXPECT_EQ(effortsReading.dimensions()[0], 2);
	EXPECT_EQ(effortsWriting.dimensions()[1], effortsReading.dimensions()[1]);

	// first method
	for (size_t k = 0; k < countsWriting.dimensions()[1]; k++) {
		EXPECT_EQ(countsWriting(0, k), countsReading(0, k));
		EXPECT_EQ(effortsWriting(0, k), effortsReading(0, k));
	}

	// third method
	for (size_t k = 0; k < countsWriting.dimensions()[1]; k++) {
		EXPECT_EQ(countsWriting(2, k), countsReading(1, k));
		EXPECT_EQ(effortsWriting(2, k), effortsReading(1, k));
	}

	// dimension names
	for (size_t k = 0; k < countsWriting.dimensions()[1]; k++) {
		EXPECT_EQ((*countsWriting.getDimensionName(1))[k], (*countsReading.getDimensionName(1))[k]);
		EXPECT_EQ((*effortsWriting.getDimensionName(1))[k], (*effortsReading.getDimensionName(1))[k]);
	}
	EXPECT_EQ((*effortsWriting.getDimensionName(0))[0], (*effortsReading.getDimensionName(0))[0]);
	EXPECT_EQ((*countsWriting.getDimensionName(0))[0], (*countsReading.getDimensionName(0))[0]);
	EXPECT_EQ((*effortsWriting.getDimensionName(0))[2], (*effortsReading.getDimensionName(0))[1]);
	EXPECT_EQ((*countsWriting.getDimensionName(0))[2], (*countsReading.getDimensionName(0))[1]);
}

TEST_F(TDataTest, threeMethods_expandEffort) {
	// first write...
	countsWriting.resize(std::array<size_t, 2>({3, 3}));
	effortsWriting.resize(std::array<size_t, 2>({3, 1}));

	for (size_t i = 0; i < countsWriting.dimensions()[0]; i++) {
		for (size_t k = 0; k < countsWriting.dimensions()[1]; k++) { countsWriting(i, k) = i + k; }
	}
	effortsWriting(0, 0) = 1.;
	effortsWriting(1, 0) = 2.;
	effortsWriting(2, 0) = 3.;
	TData::write("birp_test", countsWriting, effortsWriting);

	// should work: repeat efforts over all time points
	TData::read(countsReading, effortsReading);
	EXPECT_EQ(countsWriting.dimensions()[0], countsReading.dimensions()[0]);
	EXPECT_EQ(countsWriting.dimensions()[1], countsReading.dimensions()[1]);
	EXPECT_EQ(countsWriting.dimensions()[0], effortsReading.dimensions()[0]);
	EXPECT_EQ(countsWriting.dimensions()[1], effortsReading.dimensions()[1]);

	for (size_t i = 0; i < countsWriting.dimensions()[0]; i++) {
		for (size_t k = 0; k < countsWriting.dimensions()[1]; k++) {
			EXPECT_EQ(countsWriting(i, k), countsReading(i, k));
			EXPECT_EQ(effortsReading(i, k), i+1.0);
		}
	}
}