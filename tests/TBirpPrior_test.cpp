
#include "TBirpCore.h"
#include "TEpochModel.h"
#include "gtest/gtest.h"
#include "TBirpTests.h"

/* these are the old tests for the outdated TSpecies data structure; only keep them until I copied all relevant code to the new tests

TEST(TBirpPriorTest, readFile1) {
	// construct birp
	TBirpCoreBridge birp;
	birp._readFile("../tests/test_files/test_counts.csv");
	// test locations
	const auto locations = birp.locations();
	EXPECT_EQ(locations[0], "beach");
	EXPECT_EQ(locations[1], "river");
	EXPECT_EQ(locations[2], "mountain");
	// test species
	const auto species = birp.species();
	EXPECT_EQ(species[0].name(), "monkey");
	EXPECT_EQ(species[1].name(), "wolf");
	EXPECT_EQ(species[2].name(), "snake");
	// test method
	const auto method = birp.methods();
	EXPECT_EQ(method[0], "loud_yelling");
	EXPECT_EQ(method[1], "cameraTrap");
	EXPECT_EQ(method[2], "tracks");
	// test timepoints
	const auto times = birp.timepoints();
	EXPECT_EQ(times[0], 2010);
	EXPECT_EQ(times[1], 2015);
	EXPECT_EQ(times[2], 1017);
	EXPECT_EQ(times[3], 2014);
	EXPECT_EQ(times[4], 2012);
	EXPECT_EQ(times[5], 2029);
	EXPECT_EQ(times[6], 2020);
}

TEST(TBirpPriorTest, readFile2) {
	// construct birp
	TBirpCoreBridge birp;
	birp._readFile("../tests/test_files/test_counts2.csv");
	// test locations
	const auto locations = birp.locations();
	EXPECT_EQ(locations[0], "basecamp");
	EXPECT_EQ(locations[1], "jungle");
	EXPECT_EQ(locations[2], "Thierrys_office");
	EXPECT_EQ(locations.size(), 3);
	// test species
	const auto species = birp.species();
	EXPECT_EQ(species[0].name(), "lion");
	EXPECT_EQ(species[1].name(), "giraffe");
	EXPECT_EQ(species[2].name(), "zebra");
	EXPECT_EQ(species.size(), 3);
	// test method
	const auto method = birp.methods();
	EXPECT_EQ(method[0], "looking");
	EXPECT_EQ(method[1], "listening");
	EXPECT_EQ(method[2], "dancing");
	EXPECT_EQ(method.size(), 3);
	// test timepoints
	const auto times = birp.timepoints();
	EXPECT_EQ(times[0], 2020);
	EXPECT_EQ(times.size(), 1);
}

TEST(TBirpPriorTest, sortTimes1) {
	// construct birp
	TBirpCoreBridge birp;
	birp._readFile("../tests/test_files/sortingTest.csv");
	birp._sortTimepoints();
	// test timepoints
	const auto times = birp.timepoints();
	EXPECT_EQ(times[0], 2001);
	EXPECT_EQ(times[1], 2005);
	EXPECT_EQ(times[2], 2020);
	EXPECT_EQ(times[3], 2022);
}

TEST(TBirpPriorTest, sortTimes2) {
	TBirpCoreBridge birp;
	birp._readFile("../tests/test_files/sortingTest.csv");
	birp._sortTimepoints();
	std::vector<size_t> truthVec{2001, 2005, 2020, 2022, 2005, 2020, 2001, 2005, 2020};
	std::vector<size_t> sorted_tp{2001, 2005, 2020, 2022};
	size_t counter = 0;
	auto species   = birp.species();
	for (auto sp : species) {
		for (auto ds : sp) {
			for (auto dp : ds) {
				EXPECT_EQ(sorted_tp[dp.timepoint_id], truthVec[counter]);
				counter = counter + 1;
			}
		}
	}
}

 */

// These are the new tests for the TData files

TEST(TBirpPriorTest, readFile1) {
    // construct birp
    TBirpCoreBridge birp;
    birp._readFile("../tests/test_files/method_examineSlime_threeCovariates.csv");
    // test locations (not sorted yet?)
    const auto locations = birp.locations();
    EXPECT_EQ(locations[0], "insideSock");
    EXPECT_EQ(locations[1], "tinyRock");
    EXPECT_EQ(locations[2], "branch");
    // test speciesName TODO ??
    // const auto species = birp);
    //EXPECT_EQ(species[0].name(), "snail");
    // test method name
    // const auto method = birp.methods();
    // EXPECT_EQ(method[0].methodName(), "examineSlime_threeCovariates");
    // test timepoints
    const auto times = birp.timepoints();
    EXPECT_EQ(times[0], 3);
    EXPECT_EQ(times[1], 2);
    EXPECT_EQ(times[2], 1);
    EXPECT_EQ(times[3], 4);
}