#include <algorithm>
#include <set>

#include <gtest/gtest.h>

#include "fixtures.h"

using namespace MaptitudeTest;

TEST(FixturesTest, AtomMolIsFreshEachCall) {
    OEChem::OEGraphMol first = MakeAtomMol(6, 0.0, 0.0, 0.0);
    OEChem::OEGraphMol second = MakeAtomMol(6, 0.0, 0.0, 0.0);
    OEChem::OEAtomBase* atom = first.GetAtom(OEChem::OEHasAtomicNum(6));
    ASSERT_NE(atom, nullptr);
    atom->SetRadius(9.0);

    OEChem::OEAtomBase* other = second.GetAtom(OEChem::OEHasAtomicNum(6));
    ASSERT_NE(other, nullptr);
    EXPECT_NE(other->GetRadius(), 9.0);
}

TEST(FixturesTest, GaussianPeaksAtCentre) {
    OESystem::OEScalarGrid grid = MakeGaussianGrid(0.0, 0.0, 0.0, 0.6, 4.0, 0.5);
    float peak = 0.0f;
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        peak = std::max(peak, grid[i]);
    }
    EXPECT_NEAR(peak, 1.0f, 1e-6);
}

TEST(FixturesTest, UniformGridHasZeroSpread) {
    OESystem::OEScalarGrid grid = MakeUniformGrid(2.5f, 4.0, 0.5);
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        EXPECT_FLOAT_EQ(grid[i], 2.5f);
    }
}

TEST(FixturesTest, RampGridHasNoDuplicateValues) {
    OESystem::OEScalarGrid grid = MakeRampGrid(2.0, 1.0);
    std::set<float> seen;
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        EXPECT_TRUE(seen.insert(grid[i]).second) << "duplicate value at element " << i;
    }
}

TEST(FixturesTest, NegatedPairIsExactlyOpposite) {
    auto pair = MakeNegatedPair(MakeRampGrid(2.0, 1.0));
    ASSERT_EQ(pair.first.GetSize(), pair.second.GetSize());
    for (unsigned int i = 0; i < pair.first.GetSize(); ++i) {
        EXPECT_FLOAT_EQ(pair.first[i], -pair.second[i]);
    }
}
