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

TEST(SymOpTest, RejectsHighByteCharactersWithoutUndefinedBehavior) {
    // std::isdigit on a negative char is UB. These inputs must produce a clean
    // SymOpError, not a crash or a garbage parse.
    EXPECT_THROW(SymOp::Parse("\xc3\xa9,y,z"), SymOpError);
    EXPECT_THROW(SymOp::Parse("x,\xff,z"), SymOpError);
    EXPECT_THROW(SymOp::Parse("x,y,\x80"), SymOpError);
}

TEST(SymOpTest, RejectsRepeatedAxisInOneComponent) {
    EXPECT_THROW(SymOp::Parse("x+x,y,z"), SymOpError);
    EXPECT_THROW(SymOp::Parse("x-x,y,z"), SymOpError);
    EXPECT_THROW(SymOp::Parse("x,y+y,z"), SymOpError);
}

TEST(SymOpTest, StillAcceptsRealSymmetryOperators) {
    EXPECT_NO_THROW(SymOp::Parse("x,y,z"));
    EXPECT_NO_THROW(SymOp::Parse("-x,-y,z"));
    EXPECT_NO_THROW(SymOp::Parse("-x,y+1/2,-z+1/2"));
    EXPECT_NO_THROW(SymOp::Parse("y,x,-z"));
    EXPECT_NO_THROW(SymOp::Parse("1/2-x,1/2+y,-z"));
}

TEST(SymOpTest, EveryAcceptedOperatorRoundTripsThroughToString) {
    // ToString can only emit coefficients in {-1, 0, +1}. Once the parser
    // rejects anything it cannot serialize, parse -> ToString -> parse is
    // stable for every accepted input.
    const char* operators[] = {"x,y,z", "-x,-y,z", "y,x,-z", "-x,y+1/2,-z+1/2",
                               "1/2-x,1/2+y,-z"};
    for (const char* text : operators) {
        SymOp first = SymOp::Parse(text);
        SymOp second = SymOp::Parse(first.ToString());
        EXPECT_EQ(first.ToString(), second.ToString()) << "input: " << text;

        auto a = first.Apply(0.3, 0.4, 0.5);
        auto b = second.Apply(0.3, 0.4, 0.5);
        EXPECT_NEAR(a[0], b[0], 1e-12) << "input: " << text;
        EXPECT_NEAR(a[1], b[1], 1e-12) << "input: " << text;
        EXPECT_NEAR(a[2], b[2], 1e-12) << "input: " << text;
    }
}
