#include <gtest/gtest.h>
#include "maptitude/ScatteringFactors.h"

#include <cmath>
#include <set>

using namespace Maptitude;

TEST(ScatteringFactorsTest, CarbonExists) {
    auto* coeffs = GetScatteringFactors(6, 0);
    ASSERT_NE(coeffs, nullptr);
}

TEST(ScatteringFactorsTest, HydrogenExists) {
    auto* coeffs = GetScatteringFactors(1, 0);
    ASSERT_NE(coeffs, nullptr);
}

TEST(ScatteringFactorsTest, OxygenExists) {
    auto* coeffs = GetScatteringFactors(8, 0);
    ASSERT_NE(coeffs, nullptr);
}

TEST(ScatteringFactorsTest, NitrogenExists) {
    auto* coeffs = GetScatteringFactors(7, 0);
    ASSERT_NE(coeffs, nullptr);
}

TEST(ScatteringFactorsTest, IronChargeStates) {
    auto* fe0 = GetScatteringFactors(26, 0);
    auto* fe2 = GetScatteringFactors(26, 2);
    auto* fe3 = GetScatteringFactors(26, 3);

    ASSERT_NE(fe0, nullptr);
    ASSERT_NE(fe2, nullptr);
    ASSERT_NE(fe3, nullptr);
}

TEST(ScatteringFactorsTest, FallbackToNeutral) {
    // Request a charge state that doesn't exist for carbon
    auto* coeffs = GetScatteringFactors(6, 3);
    auto* neutral = GetScatteringFactors(6, 0);

    // Should fall back to neutral
    ASSERT_NE(coeffs, nullptr);
    EXPECT_EQ(coeffs, neutral);
}

TEST(ScatteringFactorsTest, NonexistentElement) {
    // Element 200 doesn't exist
    auto* coeffs = GetScatteringFactors(200, 0);
    EXPECT_EQ(coeffs, nullptr);
}

TEST(ScatteringFactorsTest, EvaluateAtZero) {
    auto* coeffs = GetScatteringFactors(6, 0);  // Carbon
    ASSERT_NE(coeffs, nullptr);

    // At s^2 = 0, f0 = c + sum(a_i) = total number of electrons approximately
    double f0 = coeffs->Evaluate(0.0);
    EXPECT_NEAR(f0, 6.0, 0.2);  // Carbon has 6 electrons
}

TEST(ScatteringFactorsTest, EvaluateDecaysWithS) {
    auto* coeffs = GetScatteringFactors(6, 0);  // Carbon
    ASSERT_NE(coeffs, nullptr);

    double f_low = coeffs->Evaluate(0.01);
    double f_high = coeffs->Evaluate(1.0);

    // Scattering factor should decrease with increasing s^2
    EXPECT_GT(f_low, f_high);
}

TEST(ScatteringFactorsTest, TableNotEmpty) {
    size_t count = 0;
    auto* table = GetScatteringFactorTable(count);
    ASSERT_NE(table, nullptr);
    EXPECT_GT(count, 0u);
}

// ---- Full table coverage tests ----

TEST(ScatteringFactorsTest, FullTableHas209Entries) {
    size_t count = 0;
    GetScatteringFactorTable(count);
    EXPECT_EQ(count, 209u);
}

TEST(ScatteringFactorsTest, AllCommonElementsExist) {
    // Elements commonly found in macromolecular structures
    unsigned int elements[] = {1, 6, 7, 8, 15, 16, 11, 12, 17, 19, 20, 26, 30};
    for (unsigned int z : elements) {
        auto* coeffs = GetScatteringFactors(z, 0);
        EXPECT_NE(coeffs, nullptr) << "Missing neutral atom Z=" << z;
    }
}

TEST(ScatteringFactorsTest, AllElementsZ1to98HaveNeutral) {
    // Elements Z=1..98 should all have neutral entries
    size_t count = 0;
    auto* table = GetScatteringFactorTable(count);

    std::set<unsigned int> z_values;
    for (size_t i = 0; i < count; ++i) {
        if (table[i].formal_charge == 0) {
            z_values.insert(table[i].atomic_number);
        }
    }

    // Verify we cover Z=1 through at least Z=98
    EXPECT_GE(z_values.size(), 98u);
    EXPECT_NE(z_values.find(1), z_values.end());   // H
    EXPECT_NE(z_values.find(98), z_values.end());   // Cf
}

TEST(ScatteringFactorsTest, HydrogenAnionExists) {
    auto* coeffs = GetScatteringFactors(1, -1);
    ASSERT_NE(coeffs, nullptr);
    // H- has 2 electrons
    double f0 = coeffs->Evaluate(0.0);
    EXPECT_NEAR(f0, 2.0, 0.2);
}

TEST(ScatteringFactorsTest, OxygenAnionExists) {
    auto* coeffs = GetScatteringFactors(8, -1);
    ASSERT_NE(coeffs, nullptr);
    // O- has 9 electrons
    double f0 = coeffs->Evaluate(0.0);
    EXPECT_NEAR(f0, 9.0, 0.3);
}

TEST(ScatteringFactorsTest, ChlorineAnionExists) {
    auto* coeffs = GetScatteringFactors(17, -1);
    ASSERT_NE(coeffs, nullptr);
    // Cl- has 18 electrons
    double f0 = coeffs->Evaluate(0.0);
    EXPECT_NEAR(f0, 18.0, 0.5);
}

TEST(ScatteringFactorsTest, ElectronCountAtZero) {
    // At s^2=0, f0 should approximately equal Z (for neutral atoms)
    struct TestCase {
        unsigned int z;
        double expected_electrons;
        double tolerance;
    };

    TestCase cases[] = {
        {1, 1.0, 0.2},     // H
        {6, 6.0, 0.3},     // C
        {7, 7.0, 0.3},     // N
        {8, 8.0, 0.3},     // O
        {16, 16.0, 0.5},   // S
        {26, 26.0, 1.0},   // Fe
        {30, 30.0, 1.0},   // Zn
    };

    for (const auto& tc : cases) {
        auto* coeffs = GetScatteringFactors(tc.z, 0);
        ASSERT_NE(coeffs, nullptr) << "Z=" << tc.z;
        double f0 = coeffs->Evaluate(0.0);
        EXPECT_NEAR(f0, tc.expected_electrons, tc.tolerance)
            << "Z=" << tc.z << " f0=" << f0;
    }
}

TEST(ScatteringFactorsTest, EvaluateAlwaysPositiveAtLowS) {
    size_t count = 0;
    auto* table = GetScatteringFactorTable(count);
    for (size_t i = 0; i < count; ++i) {
        double f0 = table[i].coeffs.Evaluate(0.01);
        EXPECT_GT(f0, 0.0)
            << "Z=" << static_cast<int>(table[i].atomic_number)
            << " charge=" << static_cast<int>(table[i].formal_charge);
    }
}
