#include <gtest/gtest.h>
#include "maptitude/SymOp.h"
#include "maptitude/Error.h"

#include <cmath>

using namespace Maptitude;

TEST(SymOpTest, ParseIdentity) {
    auto op = SymOp::Parse("x,y,z");

    // Rotation should be identity
    EXPECT_NEAR(op.R[0], 1.0, 1e-10);  // R[0][0]
    EXPECT_NEAR(op.R[4], 1.0, 1e-10);  // R[1][1]
    EXPECT_NEAR(op.R[8], 1.0, 1e-10);  // R[2][2]

    // Off-diagonal should be zero
    EXPECT_NEAR(op.R[1], 0.0, 1e-10);
    EXPECT_NEAR(op.R[2], 0.0, 1e-10);
    EXPECT_NEAR(op.R[3], 0.0, 1e-10);

    // Translation should be zero
    EXPECT_NEAR(op.t[0], 0.0, 1e-10);
    EXPECT_NEAR(op.t[1], 0.0, 1e-10);
    EXPECT_NEAR(op.t[2], 0.0, 1e-10);
}

TEST(SymOpTest, ParseNegation) {
    auto op = SymOp::Parse("-x,-y,-z");

    EXPECT_NEAR(op.R[0], -1.0, 1e-10);
    EXPECT_NEAR(op.R[4], -1.0, 1e-10);
    EXPECT_NEAR(op.R[8], -1.0, 1e-10);
}

TEST(SymOpTest, ParseWithFraction) {
    auto op = SymOp::Parse("-x,y+1/2,-z");

    EXPECT_NEAR(op.R[0], -1.0, 1e-10);
    EXPECT_NEAR(op.R[4], 1.0, 1e-10);
    EXPECT_NEAR(op.R[8], -1.0, 1e-10);
    EXPECT_NEAR(op.t[1], 0.5, 1e-10);
}

TEST(SymOpTest, ParseScrewAxis) {
    auto op = SymOp::Parse("x,y,z+1/4");

    EXPECT_NEAR(op.R[0], 1.0, 1e-10);
    EXPECT_NEAR(op.R[4], 1.0, 1e-10);
    EXPECT_NEAR(op.R[8], 1.0, 1e-10);
    EXPECT_NEAR(op.t[2], 0.25, 1e-10);
}

TEST(SymOpTest, ApplyIdentity) {
    auto op = SymOp::Parse("x,y,z");
    auto result = op.Apply(0.25, 0.5, 0.75);

    EXPECT_NEAR(result[0], 0.25, 1e-10);
    EXPECT_NEAR(result[1], 0.5, 1e-10);
    EXPECT_NEAR(result[2], 0.75, 1e-10);
}

TEST(SymOpTest, ApplyWithTranslation) {
    auto op = SymOp::Parse("-x,y+1/2,-z+1/4");
    auto result = op.Apply(0.25, 0.3, 0.1);

    EXPECT_NEAR(result[0], -0.25, 1e-10);
    EXPECT_NEAR(result[1], 0.8, 1e-10);
    EXPECT_NEAR(result[2], 0.15, 1e-10);
}

TEST(SymOpTest, ParseAll) {
    auto ops = SymOp::ParseAll("x,y,z\n-x,y+1/2,-z+1/2\n-x,-y,z+1/2");
    EXPECT_EQ(ops.size(), 3u);
}

TEST(SymOpTest, ParseAllSemicolonSeparated) {
    auto ops = SymOp::ParseAll("x,y,z;-x,y+1/2,-z");
    EXPECT_EQ(ops.size(), 2u);
}

TEST(SymOpTest, ParseInvalidThrows) {
    EXPECT_THROW(SymOp::Parse("x,y"), SymOpError);
    EXPECT_THROW(SymOp::Parse("x,y,z,w"), SymOpError);
}

TEST(SymOpTest, ToString) {
    auto op = SymOp::Parse("x,y,z");
    std::string s = op.ToString();
    EXPECT_FALSE(s.empty());
}
