#include <gtest/gtest.h>
#include <limits>
#include <vector>

#include "maptitude/CoverageOptions.h"
#include "maptitude/Error.h"
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

    // The grid stores float. Offsetting rho from [0,1] to [7,8] quantizes values at
    // the ~4e-7 float spacing around 7, then the correlation's mean-subtraction
    // surfaces that loss as a ~1e-9 perturbation (measured: 1.16e-9). The same
    // pipeline driven with an unmodified grid reproduces the baseline bit-exactly,
    // so this looser tolerance is specific to the offset and must not be copied.
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

    // The two Gaussians differ across the sampled support, so a positive residual follows from RSR's definition.
    EXPECT_GT(forward, 0.0);
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

    std::vector<double> scores;
    double previous = 1.0;
    for (double sigma : {0.0, 1.0, 2.0, 4.0, 8.0}) {
        CoverageOptions options;
        options.SetSigma(sigma);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        const double score = coverage(mol, grid, nullptr, options).overall;
        scores.push_back(score);
        EXPECT_LE(score, previous + 1e-12) << "coverage rose when sigma rose to " << sigma;
        previous = score;
    }

    // The atom sits on the grid maximum. At sigma 0 the threshold is the mean, and
    // the maximum is always at least the mean, so the first score is 1.0 for any map.
    // At high sigma the threshold exceeds the map's maximum, so no atom anywhere can
    // be covered. Coverage scores each atom as a binary hit and averages, so a
    // single-atom molecule can only return 1.0 or 0.0 — the property shows as one step.
    EXPECT_DOUBLE_EQ(scores.front(), 1.0);
    EXPECT_DOUBLE_EQ(scores.back(), 0.0);
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
        // EDIAm returns values in [0, 1] per Metric.h:106.
        EXPECT_GE(score, 0.0);
        EXPECT_LE(score, 1.0);
        previous = score;
    }
}

TEST(MetricAnalyticTest, QScoreRejectsAnAdaptiveSweepWithAnUnusableStep) {
    // ADAPTIVE computes step = min(grid_spacing, resolution / 7) and never consults
    // the validated options, so it bypasses every setter guard. A subnormal
    // resolution passes the `resolution <= 0.0` check and then underflows that
    // division to exactly 0.0, leaving a step of zero and a sweep that never
    // advances. Before this guard existed the call hung instead of returning.
    QScoreOptions options;
    options.SetRadialSampling(RadialSampling::ADAPTIVE);
    OESystem::OEScalarGrid grid =
        MakeGaussianGrid(0.0, 0.0, 0.0, options.GetSigma(), HALF_WIDTH, SPACING);
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);

    EXPECT_THROW(qscore(mol, grid, std::numeric_limits<double>::denorm_min(), nullptr, options),
                 GridError);
}

TEST(MetricAnalyticTest, QScoreRejectsAFixedSweepWithNoShells) {
    // Each setter's value is legal on its own; the combination yields zero shells,
    // a constant reference vector, and a nan correlation.
    QScoreOptions options;
    options.SetRadialStep(1.0);
    options.SetMaxRadius(0.1);
    OESystem::OEScalarGrid grid =
        MakeGaussianGrid(0.0, 0.0, 0.0, options.GetSigma(), HALF_WIDTH, SPACING);
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);

    EXPECT_THROW(qscore(mol, grid, RESOLUTION, nullptr, options), GridError);
}

TEST(MetricAnalyticTest, QScoreRejectsAFixedSweepWithTooManyShells) {
    // The smallest legal step against the default radius is the other end of the
    // same cross-field failure: (2.0 + 0.01) / 1e-6 is just over two million
    // shells, each allocating a full sample set per atom. Both values pass their
    // own setter, so only the consumption-point check can catch the combination.
    QScoreOptions options;
    options.SetRadialStep(MIN_RADIAL_STEP);
    OESystem::OEScalarGrid grid =
        MakeGaussianGrid(0.0, 0.0, 0.0, options.GetSigma(), HALF_WIDTH, SPACING);
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);

    EXPECT_THROW(qscore(mol, grid, RESOLUTION, nullptr, options), GridError);
}

TEST(MetricAnalyticTest, QScoreStillAcceptsTheDefaultSweep) {
    // The guard must not narrow the working configuration. This is the neutrality
    // half of the two tests above.
    QScoreOptions options;
    OESystem::OEScalarGrid grid =
        MakeGaussianGrid(0.0, 0.0, 0.0, options.GetSigma(), HALF_WIDTH, SPACING);
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);

    EXPECT_NO_THROW(qscore(mol, grid, RESOLUTION, nullptr, options));
}
