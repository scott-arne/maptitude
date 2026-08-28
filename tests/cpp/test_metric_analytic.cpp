#include <gtest/gtest.h>
#include <limits>
#include <string>
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

TEST(MetricAnalyticTest, QScoreAdaptiveFailsOnlyTheAtomWithNoRadius) {
    // ADAPTIVE derives its maximum radius as twice the atom's own, so an atom that
    // carries no radius makes the sweep unusable. That is a fact about one atom and not
    // about the request: this guard used to throw from inside the per-atom loop, so a
    // single such atom discarded the scores of every other atom in the molecule. The
    // failure now scopes to the atom, matching how an out-of-grid atom is handled.
    QScoreOptions options;
    options.SetRadialSampling(RadialSampling::ADAPTIVE);
    OESystem::OEScalarGrid grid =
        MakeGaussianGrid(0.0, 0.0, 0.0, options.GetSigma(), HALF_WIDTH, SPACING);

    // OEAssignBondiVdWRadii leaves an atomic number of zero at radius 0.0.
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    OEChem::OEAtomBase* second = mol.NewAtom(6u);
    const double second_coords[3] = {1.0, 0.0, 0.0};
    mol.SetCoords(second, second_coords);
    OEChem::OEAtomSetResidue(second, OEChem::OEAtomGetResidue(mol.GetAtom(OEChem::OEHasAtomIdx(0))));
    OEChem::OEAtomBase* radiusless = mol.NewAtom(0u);
    const double radiusless_coords[3] = {-1.0, 0.0, 0.0};
    mol.SetCoords(radiusless, radiusless_coords);
    OEChem::OEAtomSetResidue(radiusless, OEChem::OEAtomGetResidue(mol.GetAtom(OEChem::OEHasAtomIdx(0))));

    DensityScoreResult result;
    ASSERT_NO_THROW(result = qscore(mol, grid, RESOLUTION, nullptr, options));
    ASSERT_EQ(result.by_atom.count(2u), 1u) << "the radiusless atom must still be reported";
    EXPECT_TRUE(std::isnan(result.by_atom.at(2u)));
    ASSERT_EQ(result.by_atom.count(0u), 1u);
    ASSERT_EQ(result.by_atom.count(1u), 1u);
    EXPECT_FALSE(std::isnan(result.by_atom.at(0u))) << "a scorable atom was discarded with it";
    EXPECT_FALSE(std::isnan(result.by_atom.at(1u))) << "a scorable atom was discarded with it";
    EXPECT_FALSE(std::isnan(result.overall));
}

TEST(MetricAnalyticTest, QScoreRejectsAnAdaptiveStepNoAtomCanSweep) {
    // The adaptive step is min(grid_spacing, resolution / 7): a fact about the call, not
    // about any atom. A step that no atom in the molecule can sweep with is therefore a
    // caller error however the sweep check happens to be reached, and must be as loud as
    // the step floor above. Both geometries here returned a silent overall = nan once the
    // sweep check moved into the per-atom loop.
    struct Case { double resolution; double spacing; };
    for (const Case& c : {Case{30.0, 4.0}, Case{25.0, 3.5}}) {
        QScoreOptions options;
        options.SetRadialSampling(RadialSampling::ADAPTIVE);
        OESystem::OEScalarGrid grid =
            MakeGaussianGrid(0.0, 0.0, 0.0, options.GetSigma(), HALF_WIDTH, c.spacing);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);

        EXPECT_THROW(qscore(mol, grid, c.resolution, nullptr, options), GridError)
            << "resolution " << c.resolution << " A at spacing " << c.spacing << " A";
    }
}

TEST(MetricAnalyticTest, QScoreStillAcceptsAnAdaptiveSweepAtWorkingSpacing) {
    // Positive control for the rejection above: the same single carbon at a spacing and
    // resolution a real map would use must still score. The value is not pinned here --
    // that is pin_values.h's job -- only that it is a real score rather than NaN.
    QScoreOptions options;
    options.SetRadialSampling(RadialSampling::ADAPTIVE);
    OESystem::OEScalarGrid grid =
        MakeGaussianGrid(0.0, 0.0, 0.0, options.GetSigma(), HALF_WIDTH, SPACING);
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);

    DensityScoreResult result;
    ASSERT_NO_THROW(result = qscore(mol, grid, RESOLUTION, nullptr, options));
    EXPECT_FALSE(std::isnan(result.overall));
    EXPECT_GT(result.overall, 0.9);
}

TEST(MetricAnalyticTest, QScoreAdaptiveFailsOnlyTheAtomTooSmallForTheStep) {
    // The shell-count conditions read the step and the maximum radius jointly, so they
    // are call-level only when no atom can satisfy them. When one atom can, a second that
    // cannot is a fact about that atom's radius and has to scope to it -- otherwise the
    // rejection above swings too far and discards a scorable molecule.
    //
    // Radii are assigned explicitly: PrepareStructure leaves a molecule alone once any
    // heavy atom carries one, so the Bondi table never runs and the two radii here are
    // the test's own. At step 2 / 7 = 0.286 A the 0.1 A atom's sweep spans 0.2 A and
    // produces no shell, while the 1.7 A atom's spans 3.4 A and sweeps normally.
    QScoreOptions options;
    options.SetRadialSampling(RadialSampling::ADAPTIVE);
    OESystem::OEScalarGrid grid =
        MakeGaussianGrid(0.0, 0.0, 0.0, options.GetSigma(), HALF_WIDTH, SPACING);

    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    OEChem::OEAtomBase* scorable = mol.GetAtom(OEChem::OEHasAtomIdx(0));
    scorable->SetRadius(1.7);
    OEChem::OEAtomBase* too_small = mol.NewAtom(6u);
    const double too_small_coords[3] = {3.5, 0.0, 0.0};
    mol.SetCoords(too_small, too_small_coords);
    too_small->SetRadius(0.1);
    OEChem::OEAtomSetResidue(too_small, OEChem::OEAtomGetResidue(scorable));

    DensityScoreResult result;
    ASSERT_NO_THROW(result = qscore(mol, grid, RESOLUTION, nullptr, options));
    ASSERT_EQ(result.by_atom.count(1u), 1u) << "the too-small atom must still be reported";
    EXPECT_TRUE(std::isnan(result.by_atom.at(1u)));
    ASSERT_EQ(result.by_atom.count(0u), 1u);
    EXPECT_FALSE(std::isnan(result.by_atom.at(0u))) << "a scorable atom was discarded with it";
    EXPECT_FALSE(std::isnan(result.overall));
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

TEST(MetricAnalyticTest, QScoreRejectsASweepWhoseShellPointProductIsUnbounded) {
    // Bounding shells and points separately is not enough: 510000 shells is under
    // MAX_SHELLS and 10000 points is exactly MAX_NUM_POINTS, yet the fixed precompute
    // would store one 10000-point sphere per shell -- about 122 GB -- before any
    // scoring runs. Only the product check rejects this.
    QScoreOptions options;
    options.SetRadialStep(MIN_RADIAL_STEP);
    options.SetMaxRadius(0.5);
    options.SetNumPoints(MAX_NUM_POINTS);
    OESystem::OEScalarGrid grid =
        MakeGaussianGrid(0.0, 0.0, 0.0, options.GetSigma(), HALF_WIDTH, SPACING);
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);

    // Assert on the message, not just the type: 510000 shells is under MAX_SHELLS and
    // every setter accepted, so any other branch firing here would mean the arithmetic
    // moved and this test had stopped covering the product bound. Deleting the branch
    // cannot be used to check that, because the unguarded path allocates rather than
    // failing.
    try {
        qscore(mol, grid, RESOLUTION, nullptr, options);
        FAIL() << "expected GridError for a 5.1e9-sample sweep";
    } catch (const GridError& error) {
        EXPECT_NE(std::string(error.what()).find("samples per atom"), std::string::npos)
            << "rejected by the wrong branch: " << error.what();
    }
}

TEST(MetricAnalyticTest, QScoreStillAcceptsAMidBandSampleCount) {
    // Pins MAX_TOTAL_SAMPLES from below at 1608 samples: 201 shells x 8 points. Other
    // tests pin it from below too, but far lower -- QScoreCarbonSigma08 is the highest
    // of them at 4.02 shells x 16 points = 64.32 samples, and the default sweep is
    // 32.16 -- so without this case a narrowing anywhere in (64.32, 1608] leaves the
    // whole suite green. That band is where ordinary fine sweeps live, and 1000 or 200
    // is a plausible transcription slip. The product test above cannot cover it: it
    // uses a 5.1e9-sample sweep, so every value below that throws from the same branch
    // with the same message.
    QScoreOptions options;
    options.SetRadialStep(0.01);
    options.SetMaxRadius(2.0);
    OESystem::OEScalarGrid grid =
        MakeGaussianGrid(0.0, 0.0, 0.0, options.GetSigma(), HALF_WIDTH, SPACING);
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);

    EXPECT_NO_THROW(qscore(mol, grid, RESOLUTION, nullptr, options));
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
