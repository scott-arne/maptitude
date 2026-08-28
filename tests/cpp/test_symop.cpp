#include <gtest/gtest.h>
#include "maptitude/SymOp.h"
#include "maptitude/Error.h"

#include <cmath>
#include <limits>
#include <sstream>
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

TEST(SymOpTest, SerializesAZeroComponentAsAParseableZero) {
    // ParseComponent accepts "0" as a pure translation of zero, and ToString suppressed
    // both the all-zero rotation row and the sub-threshold translation -- so
    // Parse("0,y,z") serialized to ",y,z", which the empty-component guard added in
    // 339b3a5 then rejects. The serializer was emitting a string its own parser refuses.
    const SymOp op = SymOp::Parse("0,y,z");
    EXPECT_EQ(op.ToString(), "0,y,z");
    const SymOp round_tripped = SymOp::Parse(op.ToString());
    EXPECT_TRUE(op == round_tripped);
}

TEST(SymOpTest, EverySerializedOperatorReparses) {
    // The contract: every operator Parse() produces serializes to a string Parse()
    // accepts. Equality is deliberately NOT claimed here -- "1e-11,y,z" serializes to
    // "0,y,z" because of ToString's 1e-10 suppression threshold, so it reparses
    // successfully but not to an equal operator. Exact round-tripping is pinned
    // separately by UnusualTranslationsRoundTripExactly.
    const char* operators[] = {"0,y,z",  "0,0,0",      "-0,y,z",      "1e-11,y,z",
                               "1/2,y,z", "3,y,z",     "x,y,z",       "-x,-y,z",
                               "x-y,x,z", "-y,x-y,z+1/3", "x+1/13,y,z", "x-1/2,y,z"};
    for (const char* text : operators) {
        const SymOp op = SymOp::Parse(text);
        const std::string serialized = op.ToString();
        EXPECT_NO_THROW(SymOp::Parse(serialized))
            << "input: " << text << " serialized as " << serialized;
    }
}

TEST(SymOpTest, PureTranslationComponentsRoundTripExactly) {
    // A row with an all-zero rotation is the only place the translation block's
    // "nothing emitted yet" bookkeeping is observable: if the block forgets to clear the
    // flag, the row emits its translation AND the trailing zero, so "1/2,y,z" serializes
    // as "1/20,y,z" and "3,y,z" as "6/20,y,z". Both are grammatically valid, so every
    // reparses-without-throwing check stays green while the translation silently becomes
    // 0.05 and 0.3. Only exact equality on a translation-only row detects that, and no
    // other test serializes one -- StillAcceptsAPureTranslationComponent parses without
    // serializing, and every operator in the two round-trip samples has a rotation
    // coefficient in each row carrying a translation.
    const char* operators[] = {"1/2,y,z", "3,y,z", "-1/2,y,z"};
    for (const char* text : operators) {
        const SymOp first = SymOp::Parse(text);
        const SymOp second = SymOp::Parse(first.ToString());
        EXPECT_TRUE(first == second) << "input: " << text
                                     << " serialized as " << first.ToString();
    }
}

TEST(SymOpTest, LargeTranslationsRoundTripWithoutIntegerOverflow) {
    // ToString's fraction search cast std::round(frac * denom) to int. For a translation
    // at or above 2^30 that double is outside int's range, which is undefined behavior;
    // here it saturated, so "x+1073741824,y,z" serialized as "x+2147483647/2,y,z" and
    // reparsed as 1073741823.5. The corrupted output is grammatically valid, so
    // EverySerializedOperatorReparses stays green -- only exact equality detects it.
    // The two values below the boundary are controls: they must keep using the fraction
    // path rather than being pushed onto the decimal fallback by an over-broad guard.
    const char* operators[] = {"x+1073741822,y,z", "x+1073741823,y,z",
                               "x+1073741824,y,z", "x+1073741825,y,z",
                               "x+2147483647,y,z", "x+1e20,y,z",
                               "x-1e20,y,z",       "1e20,y,z"};
    for (const char* text : operators) {
        const SymOp first = SymOp::Parse(text);
        const SymOp second = SymOp::Parse(first.ToString());
        EXPECT_TRUE(first == second) << "input: " << text
                                     << " serialized as " << first.ToString();
    }
}

TEST(SymOpTest, NegativeNearZeroTranslationsSerializeWithASeparator) {
    // ToString chose the '+' separator from the sign of t[row] but chose what to print from
    // the rounded numerator. For 1e-10 < |t| < 5e-9 the numerator rounds to zero, and a
    // negative translation then emitted neither a '+' (t > 0 is false) nor a '-' (the sign
    // does not survive the cast of -0.0 to int). "x-1e-9,y,z" served as "x0/2,y,z", which
    // the term-boundary rule added at f8bedfc rejects -- the serializer emitting a string
    // its own parser refuses.
    const char* operators[] = {"x-1e-9,y,z", "x-2e-10,y,z", "x-1e-9,y-1e-9,z-1e-9",
                               "y-1e-9,x,z", "x-y-1e-9,x,z"};
    for (const char* text : operators) {
        const SymOp op = SymOp::Parse(text);
        const std::string serialized = op.ToString();
        EXPECT_NO_THROW(SymOp::Parse(serialized))
            << "input: " << text << " serialized as " << serialized;
    }
    // The band's two edges must keep behaving as documented rather than being swept up.
    EXPECT_EQ(SymOp::Parse("x-1e-10,y,z").ToString(), "x,y,z");
    EXPECT_EQ(SymOp::Parse("x-1e-8,y,z").ToString(), "x-1e-08,y,z");
}

TEST(SymOpTest, EveryTranslationMagnitudeAndSignRoundTrips) {
    // A property test over the translation space, not another point pin. Rounds 4, 5, and 6
    // each found a distinct ToString defect that every existing test missed, because each
    // test named the handful of values its own finding used. This sweeps the magnitudes that
    // select different branches -- suppressed, fraction, near-fraction, decimal, int-range
    // boundary, non-finite-product -- in both signs and in three row shapes, because the
    // round-6 defect only appears when a translation follows an already-emitted term.
    const double magnitudes[] = {
        1e-11, 1e-10, 2e-10, 1e-9, 5e-9, 1e-8, 1e-7, 1e-3,
        1.0 / 12, 1.0 / 6, 1.0 / 4, 1.0 / 3, 0.5, 2.0 / 3, 0.75, 1.0 / 13,
        1.0, 1.5, 3.0, 1e6, 1073741823.0, 1073741824.0, 2147483647.0, 1e20,
        std::numeric_limits<double>::max()};
    const char* prefixes[] = {"", "x", "x-y"};

    for (const double magnitude : magnitudes) {
        for (const double sign : {1.0, -1.0}) {
            const double value = sign * magnitude;
            for (const char* prefix : prefixes) {
                std::ostringstream component;
                component.precision(std::numeric_limits<double>::max_digits10);
                component << prefix;
                // A leading '+' on an empty prefix is legal but redundant; a negative value
                // supplies its own sign either way.
                if (*prefix != '\0' && value > 0) component << "+";
                component << value;

                const std::string text = component.str() + ",y,z";
                const SymOp op = SymOp::Parse(text);
                const std::string serialized = op.ToString();
                ASSERT_NO_THROW(SymOp::Parse(serialized))
                    << "input: " << text << " serialized as " << serialized;

                // Above both suppression thresholds -- the 1e-10 magnitude test and the
                // 1e-8 fraction tolerance, whose reach is 5e-9 at denominator 2 -- the
                // round trip must be exact, not merely parseable. That is the assertion the
                // round-5 defect needed: its corrupted output was grammatically valid.
                if (magnitude >= 1e-7) {
                    const SymOp reparsed = SymOp::Parse(serialized);
                    EXPECT_TRUE(op == reparsed)
                        << "input: " << text << " serialized as " << serialized;
                }
            }
        }
    }
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

TEST(SymOpTest, RejectsTermsWithNoSeparator) {
    // Juxtaposition is not addition. Each of these parsed silently into a DIFFERENT,
    // valid-looking operator: "2x" became x+2, "xy" became x+y with an extra rotation
    // column, "x.5" became x+1/2. A malformed CCP4 or CIF symmetry record shifted or
    // rotated atoms with no error reaching the caller.
    EXPECT_THROW(SymOp::Parse("2x,y,z"), SymOpError);
    EXPECT_THROW(SymOp::Parse("x2,y,z"), SymOpError);
    EXPECT_THROW(SymOp::Parse("2.5x,y,z"), SymOpError);
    EXPECT_THROW(SymOp::Parse(".5x,y,z"), SymOpError);
    EXPECT_THROW(SymOp::Parse("x.5,y,z"), SymOpError);
    EXPECT_THROW(SymOp::Parse("x1/2,y,z"), SymOpError);
    EXPECT_THROW(SymOp::Parse("xy,y,z"), SymOpError);
    EXPECT_THROW(SymOp::Parse("xyz,y,z"), SymOpError);
    // A coefficient on an axis: "x+2y" gained both a y column and a translation of 2.
    EXPECT_THROW(SymOp::Parse("x+2y,y,z"), SymOpError);
}

TEST(SymOpTest, StillAcceptsTermsSeparatedByASign) {
    // The boundary rule must not reject operators that do separate their terms. A leading
    // sign starts a term rather than following one, so it must not trip the check either.
    EXPECT_NO_THROW(SymOp::Parse("x+1/2,y,z"));
    EXPECT_NO_THROW(SymOp::Parse("1/2-x,y,z"));
    EXPECT_NO_THROW(SymOp::Parse("+x,y,z"));
    EXPECT_NO_THROW(SymOp::Parse("-x+y-z,y,z"));
    const SymOp op = SymOp::Parse("x-y+1/3,y,z");
    EXPECT_NEAR(op.R[0], 1.0, 1e-12);
    EXPECT_NEAR(op.R[1], -1.0, 1e-12);
    EXPECT_NEAR(op.R[2], 0.0, 1e-12);
    EXPECT_NEAR(op.t[0], 1.0 / 3.0, 1e-12);
}
