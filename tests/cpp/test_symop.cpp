#include <gtest/gtest.h>
#include "maptitude/SymOp.h"
#include "maptitude/Error.h"

#include <cmath>
#include <string>

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

TEST(SymOpTest, StillAcceptsTrigonalAndHexagonalOperators) {
    // Space groups 143-194 have rows with two non-zero entries, so "at most one
    // non-zero per row" is NOT the rule the repeated-axis guard enforces -- it
    // rejects a repeated axis, which is a different and much narrower claim.
    // Tightening it to one-per-row would silently lose every trigonal and
    // hexagonal setting, so pin the distinction rather than leave it to a comment.
    EXPECT_NO_THROW(SymOp::Parse("x-y,x,z"));      // P3, 3+ about c
    EXPECT_NO_THROW(SymOp::Parse("-y,x-y,z"));     // P3, 3- about c
    EXPECT_NO_THROW(SymOp::Parse("y-x,-x,z"));     // P6, 6- about c
    EXPECT_NO_THROW(SymOp::Parse("-x+y,-x,z+1/3"));  // P3(1), screw component

    // The rotation really is populated, not silently dropped: "x-y,x,z" sends
    // (u, v, w) to (u - v, u, w).
    const SymOp op = SymOp::Parse("x-y,x,z");
    const auto image = op.Apply(0.3, 0.4, 0.5);
    EXPECT_NEAR(image[0], -0.1, 1e-12);
    EXPECT_NEAR(image[1], 0.3, 1e-12);
    EXPECT_NEAR(image[2], 0.5, 1e-12);
}

TEST(SymOpTest, RepresentativeOperatorsRoundTripThroughToString) {
    // A representative sample, not a proof over the whole acceptance set -- the
    // universal claim belongs to UnusualTranslationsRoundTripExactly, which covers
    // the translations this list does not reach. Compare the operators themselves,
    // not two ToString() outputs: string equality survives information lost in the
    // FIRST serialization, because the second serialization loses it identically.
    const char* operators[] = {"x,y,z", "-x,-y,z", "y,x,-z", "-x,y+1/2,-z+1/2",
                               "1/2-x,1/2+y,-z", "x-y,x,z", "-y,x-y,z+1/3"};
    for (const char* text : operators) {
        SymOp first = SymOp::Parse(text);
        SymOp second = SymOp::Parse(first.ToString());
        EXPECT_TRUE(first == second) << "input: " << text
                                     << " serialized as " << first.ToString();
    }
}

TEST(SymOpTest, MalformedNumbersRaiseSymOpErrorNotStdExceptions) {
    // std::stod signals failure with std::invalid_argument and std::out_of_range.
    // Neither is a SymOpError, so both crossed the Python boundary through the
    // generic std::exception arm as RuntimeError -- invisible to a caller catching
    // SymOpError around a malformed CCP4 or CIF operator.
    EXPECT_THROW(SymOp::Parse(".,y,z"), SymOpError);
    EXPECT_THROW(SymOp::Parse("x+1/,y,z"), SymOpError);
    EXPECT_THROW(SymOp::Parse("x+1e309,y,z"), SymOpError);
}

TEST(SymOpTest, RejectsZeroDenominator) {
    // "x+1/0" did not fail: it parsed to t[0] = inf and serialized as "x+inf,y,z",
    // carrying a non-finite translation into density expansion.
    //
    // Assert on the message, not just the type. The trailing !std::isfinite(trans)
    // check would reject this input on its own, so a type-only assertion stays green
    // when the zero-denominator branch is deleted and the diagnostic silently
    // degrades to "Non-finite translation". Pinning the message is what makes the
    // specific branch load-bearing.
    try {
        SymOp::Parse("x+1/0,y,z");
        FAIL() << "expected SymOpError for a zero denominator";
    } catch (const SymOpError& error) {
        EXPECT_NE(std::string(error.what()).find("Zero denominator"), std::string::npos)
            << "rejected by the wrong branch: " << error.what();
    }
}

TEST(SymOpTest, RejectsMalformedSignSequences) {
    // The loop carried 'sign' as mutable state with no record of whether a term was
    // owed, so a sign with nothing to apply to was silently absorbed.
    EXPECT_THROW(SymOp::Parse("--x,y,z"), SymOpError);  // parsed as -x
    EXPECT_THROW(SymOp::Parse("x+,y,z"), SymOpError);   // trailing '+' discarded
    EXPECT_THROW(SymOp::Parse(",y,z"), SymOpError);     // empty component, all-zero row
}

TEST(SymOpTest, StillAcceptsAPureTranslationComponent) {
    // An all-zero rotation row is explicitly still legal: "1/2" is a translation, not
    // a malformed component. The empty-component rejection must not reach it.
    const SymOp op = SymOp::Parse("1/2,y,z");
    EXPECT_NEAR(op.t[0], 0.5, 1e-12);
    EXPECT_NEAR(op.R[0], 0.0, 1e-12);
    EXPECT_NEAR(op.R[1], 0.0, 1e-12);
    EXPECT_NEAR(op.R[2], 0.0, 1e-12);
}

TEST(SymOpTest, UnusualTranslationsRoundTripExactly) {
    // ToString only recognizes denominators 2..12. Everything else fell back to the
    // stream's default six significant digits, so 1/13 serialized as 0.0769231 and
    // reparsed to a different double -- a lossy round trip the seven-case sample
    // below never reached.
    const SymOp first = SymOp::Parse("x+1/13,y,z");
    const SymOp second = SymOp::Parse(first.ToString());
    EXPECT_TRUE(first == second) << "serialized as " << first.ToString();
    EXPECT_DOUBLE_EQ(first.t[0], second.t[0]);
}

TEST(SymOpTest, RejectsNonFiniteAndNonDecimalNumberTokens) {
    // std::stod accepts "inf", "infinity", "nan", and C99 hex floats. The named
    // literals reach the parser only in the denominator, because 'i' and 'n' never
    // enter the digit branch -- and 1/inf is a finite 0.0, so every downstream
    // finiteness check passes and the component silently becomes the identity.
    EXPECT_THROW(SymOp::Parse("x+1/inf,y,z"), SymOpError);
    EXPECT_THROW(SymOp::Parse("x+1/infinity,y,z"), SymOpError);
    EXPECT_THROW(SymOp::Parse("x+1/-inf,y,z"), SymOpError);
    // A hex float starts with a digit, so only the consumed-token scan rejects it.
    EXPECT_THROW(SymOp::Parse("x+0x10,y,z"), SymOpError);
    EXPECT_THROW(SymOp::Parse("x+1/0x10,y,z"), SymOpError);
}

TEST(SymOpTest, RejectsASignedNumberToken) {
    // The consumed-token scan allows '+' and '-' because an exponent needs them, so
    // it cannot reject a leading sign; only the leading-character test can. Signs
    // belong to the sign branch, which is what tracks pending_sign -- a sign that
    // reaches ParseNumber has bypassed that state entirely. Without this test the
    // leading-character test is unpinned and "x+1/+2,y,z" quietly parses as x+1/2.
    EXPECT_THROW(SymOp::Parse("x+1/+2,y,z"), SymOpError);
    EXPECT_THROW(SymOp::Parse("x+1/-2,y,z"), SymOpError);
}

TEST(SymOpTest, RejectsATranslationThatOverflowsWhileSumming) {
    // Every operand here is finite and every token is well formed, so no per-token
    // check fires: the sum is what overflows. This is the only input the trailing
    // !std::isfinite(trans) guard rejects on its own, and the round-2 implementer
    // correctly reported that nothing exercised it.
    EXPECT_THROW(SymOp::Parse("x+1e308+1e308,y,z"), SymOpError);
}

TEST(SymOpTest, StillAcceptsEveryDecimalFormTheGrammarIntends) {
    // The token validation must not narrow the numbers real operators use.
    EXPECT_NO_THROW(SymOp::Parse("x+1/2,y,z"));
    EXPECT_NO_THROW(SymOp::Parse("x+0.5,y,z"));
    EXPECT_NO_THROW(SymOp::Parse("x+.5,y,z"));
    EXPECT_NO_THROW(SymOp::Parse("x+5e-1,y,z"));
    const SymOp op = SymOp::Parse("x+5e-1,y,z");
    EXPECT_NEAR(op.t[0], 0.5, 1e-12);
}
