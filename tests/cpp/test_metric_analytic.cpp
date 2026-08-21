#include <gtest/gtest.h>

#include "maptitude/CoverageOptions.h"
#include "maptitude/Metric.h"
#include "maptitude/QScoreOptions.h"
#include "maptitude/RsccOptions.h"
#include "maptitude/RsrOptions.h"

#include "fixtures.h"

using namespace Maptitude;
using namespace MaptitudeTest;

namespace {
constexpr double HALF_WIDTH = 4.0;
constexpr double SPACING = 0.5;
constexpr double RESOLUTION = 2.0;
}  // namespace

TEST(MetricAnalyticTest, RsccOfMapWithItselfIsOne) {
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    OESystem::OEScalarGrid grid = MakeGaussianGrid(0.0, 0.0, 0.0, 1.0, HALF_WIDTH, SPACING);
    DensityScoreResult result = rscc(mol, grid, RESOLUTION, nullptr, &grid);
    EXPECT_NEAR(result.overall, 1.0, 1e-9);
}

TEST(MetricAnalyticTest, RsccOfExactNegationIsMinusOne) {
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    auto pair = MakeNegatedPair(MakeGaussianGrid(0.0, 0.0, 0.0, 1.0, HALF_WIDTH, SPACING));
    DensityScoreResult result = rscc(mol, pair.first, RESOLUTION, nullptr, &pair.second);
    EXPECT_NEAR(result.overall, -1.0, 1e-9);
}

TEST(MetricAnalyticTest, RsccIsInvariantUnderPositiveAffineRescaling) {
    OESystem::OEScalarGrid grid = MakeGaussianGrid(0.0, 0.0, 0.0, 1.0, HALF_WIDTH, SPACING);
    OESystem::OEScalarGrid calc = MakeGaussianGrid(0.3, 0.0, 0.0, 1.0, HALF_WIDTH, SPACING);

    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    const double baseline = rscc(mol, grid, RESOLUTION, nullptr, &calc).overall;

    // rho -> 3*rho + 7 must not move a correlation coefficient.
    OESystem::OEScalarGrid rescaled(grid);
    for (unsigned int i = 0; i < rescaled.GetSize(); ++i) {
        rescaled[i] = 3.0f * rescaled[i] + 7.0f;
    }
    OEChem::OEGraphMol mol2 = MakeAtomMol(6, 0.0, 0.0, 0.0);
    const double rescaled_score = rscc(mol2, rescaled, RESOLUTION, nullptr, &calc).overall;

    // Tolerance accounts for accumulated floating-point errors in the correlation
    // computation over the 4913-element grid (grid rescaling, trilinear sampling,
    // computed-map generation, mean/variance calculation, and the final correlation).
    EXPECT_NEAR(rescaled_score, baseline, 1e-8);
}

TEST(MetricAnalyticTest, RsrOfIdenticalMapsIsZero) {
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    OESystem::OEScalarGrid grid = MakeGaussianGrid(0.0, 0.0, 0.0, 1.0, HALF_WIDTH, SPACING);
    DensityScoreResult result = rsr(mol, grid, RESOLUTION, nullptr, &grid);
    EXPECT_NEAR(result.overall, 0.0, 1e-9);
}

TEST(MetricAnalyticTest, RsrIsSymmetricInItsTwoMaps) {
    OESystem::OEScalarGrid a = MakeGaussianGrid(0.0, 0.0, 0.0, 1.0, HALF_WIDTH, SPACING);
    OESystem::OEScalarGrid b = MakeGaussianGrid(0.5, 0.0, 0.0, 1.2, HALF_WIDTH, SPACING);

    OEChem::OEGraphMol mol_ab = MakeAtomMol(6, 0.0, 0.0, 0.0);
    OEChem::OEGraphMol mol_ba = MakeAtomMol(6, 0.0, 0.0, 0.0);
    const double forward = rsr(mol_ab, a, RESOLUTION, nullptr, &b).overall;
    const double reverse = rsr(mol_ba, b, RESOLUTION, nullptr, &a).overall;

    EXPECT_NEAR(forward, reverse, 1e-9);
}

TEST(MetricAnalyticTest, CoverageOnUniformMapIsOne) {
    // A uniform map has stddev == 0, so threshold == mean for every sigma, and
    // the scoring comparison is >=. Every sampled point is a hit. This is a
    // property of the definition, not a bug: do not "fix" it to 0.
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    OESystem::OEScalarGrid grid = MakeUniformGrid(1.0f, HALF_WIDTH, SPACING);
    DensityScoreResult result = coverage(mol, grid);
    EXPECT_NEAR(result.overall, 1.0, 1e-12);
}

TEST(MetricAnalyticTest, CoverageIsMonotonicNonIncreasingInSigma) {
    OESystem::OEScalarGrid grid = MakeGaussianGrid(0.0, 0.0, 0.0, 1.5, HALF_WIDTH, SPACING);

    double previous = 1.0;
    for (double sigma : {0.0, 0.5, 1.0, 2.0, 3.0}) {
        CoverageOptions options;
        options.SetSigma(sigma);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        const double score = coverage(mol, grid, nullptr, options).overall;
        EXPECT_LE(score, previous + 1e-12) << "coverage rose when sigma rose to " << sigma;
        previous = score;
    }
}

TEST(MetricAnalyticTest, QScoreIsNearOneForTheReferenceGaussian) {
    QScoreOptions options;
    OESystem::OEScalarGrid grid =
        MakeGaussianGrid(0.0, 0.0, 0.0, options.GetSigma(), HALF_WIDTH, SPACING);
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);

    // Trilinear interpolation off the reference-Gaussian nodes keeps this shy of
    // an exact 1.0; the correlation is invariant under the map normalization.
    DensityScoreResult result = qscore(mol, grid, RESOLUTION, nullptr, options);
    EXPECT_GT(result.overall, 0.99);
    EXPECT_LE(result.overall, 1.0 + 1e-9);
}

TEST(MetricAnalyticTest, EdiamIsMonotonicNonDecreasingInUniformLevel) {
    double previous = -1e30;
    for (float level : {0.1f, 0.5f, 1.0f, 2.0f}) {
        OESystem::OEScalarGrid grid = MakeUniformGrid(level, HALF_WIDTH, SPACING);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        const double score = ediam(mol, grid, RESOLUTION).overall;
        EXPECT_GE(score, previous - 1e-12) << "EDIAm fell when the level rose to " << level;
        EXPECT_GE(score, -1.0);
        EXPECT_LE(score, 1.0);
        previous = score;
    }
}
