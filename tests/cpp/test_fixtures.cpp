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
    OESystem::OESkewGrid grid = MakeGaussianGrid(cx, cy, cz, sigma, 4.0, 0.5);
    const float* values = grid.GetValues();

    float peak = 0.0f;
    float peak_x = 0.0f, peak_y = 0.0f, peak_z = 0.0f;
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        float gx, gy, gz;
        grid.ElementToSpatialCoord(i, gx, gy, gz);
        const double dx = gx - cx, dy = gy - cy, dz = gz - cz;
        const double expected = std::exp(-(dx * dx + dy * dy + dz * dz) / (2.0 * sigma * sigma));
        EXPECT_NEAR(values[i], static_cast<float>(expected), 1e-5);

        if (values[i] > peak) {
            peak = values[i];
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
    OESystem::OESkewGrid grid = MakeUniformGrid(2.5f, 4.0, 0.5);
    const float* values = grid.GetValues();
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        EXPECT_FLOAT_EQ(values[i], 2.5f);
    }
}

TEST(FixturesTest, RampGridHasNoDuplicateValues) {
    OESystem::OESkewGrid grid = MakeRampGrid(2.0, 1.0);
    const float* values = grid.GetValues();
    std::set<float> seen;
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        EXPECT_TRUE(seen.insert(values[i]).second) << "duplicate value at element " << i;
    }
}

TEST(FixturesTest, NegatedPairIsExactlyOpposite) {
    auto pair = MakeNegatedPair(MakeRampGrid(2.0, 1.0));
    ASSERT_EQ(pair.first.GetSize(), pair.second.GetSize());
    const float* first_values = pair.first.GetValues();
    const float* second_values = pair.second.GetValues();
    for (unsigned int i = 0; i < pair.first.GetSize(); ++i) {
        EXPECT_FLOAT_EQ(first_values[i], -second_values[i]);
    }
}

TEST(FixturesTest, RampGridAcceptsLargestValidGeometry) {
    OESystem::OESkewGrid grid = MakeRampGrid(4.5, 1.0);
    EXPECT_EQ(grid.GetXDim(), 10u);
    EXPECT_EQ(grid.GetYDim(), 10u);
    EXPECT_EQ(grid.GetZDim(), 10u);
    EXPECT_EQ(grid.GetSize(), 1000u);

    const float* values = grid.GetValues();
    std::set<float> seen;
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        EXPECT_TRUE(seen.insert(values[i]).second) << "duplicate value at element " << i;
    }
}

TEST(FixturesTest, RampGridRejectsOversizedGeometry) {
    EXPECT_THROW(MakeRampGrid(5.0, 1.0), std::invalid_argument);
}

TEST(FixturesTest, EmptyGridIsZeroFilled) {
    OESystem::OESkewGrid grid = MakeEmptyGrid(4.0, 0.5);
    EXPECT_EQ(grid.GetXDim(), 17u);
    EXPECT_EQ(grid.GetYDim(), 17u);
    EXPECT_EQ(grid.GetZDim(), 17u);
    EXPECT_EQ(grid.GetSize(), 4913u);
    EXPECT_FLOAT_EQ(grid.GetSpacing(), 0.5f);

    // A change in SetUnitCell/SetMid semantics would move the nodes, and the
    // node coordinates are what the metric pins read: pin the two extreme
    // corners so such a change surfaces here rather than as an unattributable
    // drift in a pin. The cell is not incidental either -- same_grid_geometry
    // compares its six parameters, and generate_pins.cpp feeds fixture grids to
    // combine_maps (src/GridOps.cpp:108) and diff_to_calc (:140), both of which
    // call it -- so those six are pinned below as well.
    float x0, y0, z0, x1, y1, z1;
    grid.ElementToSpatialCoord(0u, x0, y0, z0);
    grid.ElementToSpatialCoord(grid.GetSize() - 1u, x1, y1, z1);
    EXPECT_NEAR(x0, -4.0f, 1e-5);
    EXPECT_NEAR(y0, -4.0f, 1e-5);
    EXPECT_NEAR(z0, -4.0f, 1e-5);
    EXPECT_NEAR(x1, 4.0f, 1e-5);
    EXPECT_NEAR(y1, 4.0f, 1e-5);
    EXPECT_NEAR(z1, 4.0f, 1e-5);

    // 17 nodes at 0.5 A give an 8.5 A cell edge, one spacing more than the
    // 8.0 A the nodes themselves span.
    float ca, cb, cc, alpha, beta, gamma;
    ASSERT_TRUE(grid.GetUnitCell(ca, cb, cc, alpha, beta, gamma));
    EXPECT_NEAR(ca, 8.5f, 1e-5);
    EXPECT_NEAR(cb, 8.5f, 1e-5);
    EXPECT_NEAR(cc, 8.5f, 1e-5);
    EXPECT_NEAR(alpha, 90.0f, 1e-5);
    EXPECT_NEAR(beta, 90.0f, 1e-5);
    EXPECT_NEAR(gamma, 90.0f, 1e-5);

    const float* values = grid.GetValues();
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        EXPECT_FLOAT_EQ(values[i], 0.0f);
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
