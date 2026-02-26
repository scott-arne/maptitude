#include <gtest/gtest.h>
#include "maptitude/DensityScoreResult.h"

using namespace Maptitude;

TEST(DensityScoreResultTest, DefaultValues) {
    DensityScoreResult result;
    EXPECT_DOUBLE_EQ(result.overall, 0.0);
    EXPECT_TRUE(result.by_residue.empty());
    EXPECT_TRUE(result.by_atom.empty());
}

TEST(DensityScoreResultTest, ToStringEmpty) {
    DensityScoreResult result;
    result.overall = 0.95;
    std::string s = result.ToString();
    EXPECT_NE(s.find("0.950"), std::string::npos);
    EXPECT_NE(s.find("residues=0"), std::string::npos);
    EXPECT_NE(s.find("atoms=0"), std::string::npos);
}

TEST(DensityScoreResultTest, ToStringWithData) {
    DensityScoreResult result;
    result.overall = 0.875;

    Residue r1("ALA", 1, "A");
    Residue r2("GLY", 2, "A");
    result.by_residue[r1] = 0.9;
    result.by_residue[r2] = 0.85;

    result.by_atom[0] = 0.92;
    result.by_atom[1] = 0.88;
    result.by_atom[2] = 0.84;

    std::string s = result.ToString();
    EXPECT_NE(s.find("0.875"), std::string::npos);
    EXPECT_NE(s.find("residues=2"), std::string::npos);
    EXPECT_NE(s.find("atoms=3"), std::string::npos);
}

TEST(DensityScoreResultTest, ToStringNegativeOverall) {
    DensityScoreResult result;
    result.overall = -0.123;
    std::string s = result.ToString();
    EXPECT_NE(s.find("-0.123"), std::string::npos);
}
