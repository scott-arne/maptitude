#include <algorithm>
#include <set>
#include <stdexcept>

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

TEST(FixturesTest, GaussianMatchesClosedForm) {
    const double cx = 1.0, cy = -1.0, cz = 2.0, sigma = 0.6;
    OESystem::OEScalarGrid grid = MakeGaussianGrid(cx, cy, cz, sigma, 4.0, 0.5);

    float peak = 0.0f;
    float peak_x = 0.0f, peak_y = 0.0f, peak_z = 0.0f;
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        float gx, gy, gz;
        grid.ElementToSpatialCoord(i, gx, gy, gz);
        const double dx = gx - cx, dy = gy - cy, dz = gz - cz;
        const double expected = std::exp(-(dx * dx + dy * dy + dz * dz) / (2.0 * sigma * sigma));
        EXPECT_NEAR(grid[i], static_cast<float>(expected), 1e-5);

        if (grid[i] > peak) {
            peak = grid[i];
            peak_x = gx;
            peak_y = gy;
            peak_z = gz;
        }
    }

    EXPECT_NEAR(peak, 1.0f, 1e-5);
    EXPECT_NEAR(peak_x, cx, 1e-5);
    EXPECT_NEAR(peak_y, cy, 1e-5);
    EXPECT_NEAR(peak_z, cz, 1e-5);
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

TEST(FixturesTest, RampGridAcceptsLargestValidGeometry) {
    OESystem::OEScalarGrid grid = MakeRampGrid(4.5, 1.0);
    EXPECT_EQ(grid.GetXDim(), 10u);
    EXPECT_EQ(grid.GetYDim(), 10u);
    EXPECT_EQ(grid.GetZDim(), 10u);
    EXPECT_EQ(grid.GetSize(), 1000u);

    std::set<float> seen;
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        EXPECT_TRUE(seen.insert(grid[i]).second) << "duplicate value at element " << i;
    }
}

TEST(FixturesTest, RampGridRejectsOversizedGeometry) {
    EXPECT_THROW(MakeRampGrid(5.0, 1.0), std::invalid_argument);
}

TEST(FixturesTest, EmptyGridIsZeroFilled) {
    OESystem::OEScalarGrid grid = MakeEmptyGrid(4.0, 0.5);
    EXPECT_EQ(grid.GetXDim(), 17u);
    EXPECT_EQ(grid.GetYDim(), 17u);
    EXPECT_EQ(grid.GetZDim(), 17u);
    EXPECT_EQ(grid.GetSize(), 4913u);
    EXPECT_FLOAT_EQ(grid.GetSpacing(), 0.5f);

    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        EXPECT_FLOAT_EQ(grid[i], 0.0f);
    }
}

TEST(FixturesTest, AtomMolHasResidueMetadata) {
    OEChem::OEGraphMol mol = MakeAtomMol(6, 1.0, 2.0, 3.0);
    EXPECT_EQ(mol.NumAtoms(), 1u);

    OEChem::OEAtomBase* atom = mol.GetAtom(OEChem::OEHasAtomicNum(6));
    ASSERT_NE(atom, nullptr);

    EXPECT_EQ(atom->GetAtomicNum(), 6);

    double coords[3];
    mol.GetCoords(atom, coords);
    EXPECT_NEAR(coords[0], 1.0, 1e-6);
    EXPECT_NEAR(coords[1], 2.0, 1e-6);
    EXPECT_NEAR(coords[2], 3.0, 1e-6);

    const OEChem::OEResidue &residue = OEChem::OEAtomGetResidue(atom);
    EXPECT_STREQ(residue.GetName(), "LIG");
    EXPECT_EQ(residue.GetResidueNumber(), 1);
    EXPECT_EQ(residue.GetChainID(), 'A');
    EXPECT_FLOAT_EQ(residue.GetBFactor(), 0.0f);
}
