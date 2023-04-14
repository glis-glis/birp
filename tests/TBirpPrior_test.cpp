
#include "TBirpCore.h"
#include "TEpochModel.h"
#include "gtest/gtest.h"
#include "TBirpTests.h"

TEST(TBirpPriorTest, readFile1) {
    // construct birp
    TBirpCoreBridge birp;
    birp._readFile("../tests/test_files/method_examineSlime_threeCovariates.csv");
    // test locations (not sorted yet?)
    const auto locations = birp.locations();
    EXPECT_EQ(locations[0], "insideSock");
    EXPECT_EQ(locations[1], "tinyRock");
    EXPECT_EQ(locations[2], "branch");
    // test method name
    const auto method = birp.methods();
    // EXPECT_EQ(method[0].methodName(), "examineSlime_threeCovariates");
    // test species name
    EXPECT_EQ(birp.speciesNames()[method[0].speciesIDsinUniqueContainer()[0]], "snail");
    // test timepoints
    const auto times = birp.timepoints();
    EXPECT_EQ(times[0], 3);
    EXPECT_EQ(times[1], 2);
    EXPECT_EQ(times[2], 1);
    EXPECT_EQ(times[3], 4);
}

TEST(TBirpPriorTest, readFile2) {
    // construct birp
    TBirpCoreBridge birp;
    birp._readFile("../tests/test_files/method_loudYelling.csv");
    // test locations
    const auto locations = birp.locations();
    EXPECT_EQ(locations[0], "mountain");
    EXPECT_EQ(locations[1], "river");
    EXPECT_EQ(locations.size(), 2);
    // test method
    const auto method = birp.methods();
    // EXPECT_EQ(method[0].methodName(), "loudYelling");
    // test species name
    EXPECT_EQ(birp.speciesNames()[method[0].speciesIDsinUniqueContainer()[0]], "lion");
// test timepoints
    const auto times = birp.timepoints();
    EXPECT_EQ(times[0], 1);
    EXPECT_EQ(times[1], 4);
    EXPECT_EQ(times[2], 3);
    EXPECT_EQ(times[3], 2);
    EXPECT_EQ(times.size(), 4);
}

TEST(TBirpPriorTest, sortTimes1) {
    // construct birp
    TBirpCoreBridge birp;
    birp._readFile("../tests/test_files/method_examineSlime_threeCovariates.csv");
    birp._sortTimepoints();
    // test timepoints
    const auto times = birp.timepoints();
    EXPECT_EQ(times[0], 1);
    EXPECT_EQ(times[1], 2);
    EXPECT_EQ(times[2], 3);
    EXPECT_EQ(times[3], 4);
    EXPECT_EQ(times[4], 5);
}

TEST(TBirpPriorTest, sortTimes2) {
    TBirpCoreBridge birp;
    birp._readFile("../tests/test_files/method_cameraTrap.csv");
    birp.timepoints().print();
    birp._sortTimepoints();
    birp.timepoints().print();
    std::vector<size_t> truthVec{2001, 2005, 2020, 2022, 2001, 2005, 2007, 2020};
    std::vector<size_t> sorted_tp{2001, 2005, 2007, 2020, 2022};
    size_t counter = 0;
    const auto& methods   = birp.methods();
    for (auto m : methods) {
        for (auto loc : m) {
            for (const auto& tp : loc) {
                EXPECT_EQ(sorted_tp[tp.timepoint_id()], truthVec[counter]);
                //std::cout << tp.timepoint_id() << ": " << sorted_tp[tp.timepoint_id()] << " vs " << truthVec[counter] << std::endl;
                counter = counter + 1;
            }
        }
    }
}