/// Tier 2 characterization ("pin") tests.
///
/// These assert only that today's numbers equal yesterday's. They encode no
/// claim of correctness. A failure here during Phase 1 means a fix changed a
/// score -- reconcile it against the documented exceptions in the spec before
/// regenerating pin_values.h.
#include <algorithm>
#include <cmath>

#include <gtest/gtest.h>

#include "maptitude/CoverageOptions.h"
#include "maptitude/Metric.h"
#include "maptitude/QScoreOptions.h"
#include "maptitude/RsccOptions.h"
#include "maptitude/RsrOptions.h"

#include "fixtures.h"
#include "pin_values.h"

using namespace Maptitude;
using namespace MaptitudeTest;

namespace {

constexpr double HALF_WIDTH = 6.0;
constexpr double SPACING = 0.5;
constexpr double RESOLUTION = 2.0;
constexpr double PIN_RELATIVE_TOLERANCE = 1e-12;

OESystem::OEScalarGrid ObsGrid() {
    return MakeGaussianGrid(0.0, 0.0, 0.0, 1.0, HALF_WIDTH, SPACING);
}

OESystem::OEScalarGrid CalcGrid() {
    return MakeGaussianGrid(0.25, -0.1, 0.15, 1.1, HALF_WIDTH, SPACING);
}

/// Compare against a pinned value with a tolerance scaled by the larger of 1.0
/// or |pinned|. All metrics here are bounded in [-1, 1], so all current pins use
/// the absolute branch (scale = 1.0). The relative branch exists for Task 5's
/// DensityCalculator pins, which can exceed 1.0.
void ExpectPinned(double actual, double pinned) {
    const double scale = std::max(1.0, std::abs(pinned));
    EXPECT_NEAR(actual, pinned, PIN_RELATIVE_TOLERANCE * scale);
}

}  // namespace

TEST(MetricCharacterizationTest, RsccCarbonBinned) {
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    const OESystem::OEScalarGrid obs = ObsGrid();
    const OESystem::OEScalarGrid calc = CalcGrid();
    ExpectPinned(rscc(mol, obs, RESOLUTION, nullptr, &calc).overall, MaptitudePins::RSCC_CARBON_BINNED);
}

TEST(MetricCharacterizationTest, RsccCarbonFixed) {
    RsccOptions options;
    options.SetAtomRadiusMethod(AtomRadius::FIXED);
    options.SetFixedAtomRadius(1.5);
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    const OESystem::OEScalarGrid obs = ObsGrid();
    const OESystem::OEScalarGrid calc = CalcGrid();
    ExpectPinned(rscc(mol, obs, RESOLUTION, nullptr, &calc, options).overall,
                 MaptitudePins::RSCC_CARBON_FIXED);
}

TEST(MetricCharacterizationTest, RsccCarbonScaled) {
    // Carbon's 1.7 A Bondi radius times the plan's original 1.2 scaling gave
    // 2.04 A, indistinguishable from the 2.00 A binned radius at 0.5 A spacing.
    // Scaling 1.5 reaches a strictly larger point set.
    RsccOptions options;
    options.SetAtomRadiusMethod(AtomRadius::SCALED);
    options.SetAtomRadiusScaling(1.5);
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    const OESystem::OEScalarGrid obs = ObsGrid();
    const OESystem::OEScalarGrid calc = CalcGrid();
    ExpectPinned(rscc(mol, obs, RESOLUTION, nullptr, &calc, options).overall,
                 MaptitudePins::RSCC_CARBON_SCALED);
}

TEST(MetricCharacterizationTest, RsccOxygenOffset) {
    OEChem::OEGraphMol mol = MakeAtomMol(8, 0.3, 0.0, 0.0);
    const OESystem::OEScalarGrid obs = ObsGrid();
    const OESystem::OEScalarGrid calc = CalcGrid();
    ExpectPinned(rscc(mol, obs, RESOLUTION, nullptr, &calc).overall, MaptitudePins::RSCC_OXYGEN_OFFSET);
}

TEST(MetricCharacterizationTest, RsrCarbonBinned) {
    // RsrOptions defaults to ADAPTIVE, so the binned path must be selected
    // explicitly -- otherwise this pin silently duplicates RSR_CARBON_ADAPTIVE.
    RsrOptions options;
    options.SetAtomRadiusMethod(AtomRadius::BINNED);
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    const OESystem::OEScalarGrid obs = ObsGrid();
    const OESystem::OEScalarGrid calc = CalcGrid();
    ExpectPinned(rsr(mol, obs, RESOLUTION, nullptr, &calc, options).overall, MaptitudePins::RSR_CARBON_BINNED);
}

TEST(MetricCharacterizationTest, RsrCarbonAdaptive) {
    RsrOptions options;
    options.SetAtomRadiusMethod(AtomRadius::ADAPTIVE);
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    const OESystem::OEScalarGrid obs = ObsGrid();
    const OESystem::OEScalarGrid calc = CalcGrid();
    ExpectPinned(rsr(mol, obs, RESOLUTION, nullptr, &calc, options).overall,
                 MaptitudePins::RSR_CARBON_ADAPTIVE);
}

TEST(MetricCharacterizationTest, QScoreCarbonDefault) {
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    ExpectPinned(qscore(mol, ObsGrid(), RESOLUTION).overall, MaptitudePins::QSCORE_CARBON_DEFAULT);
}

TEST(MetricCharacterizationTest, QScoreCarbonSigma08) {
    QScoreOptions options;
    options.SetSigma(0.8);
    options.SetNumPoints(16);
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    ExpectPinned(qscore(mol, ObsGrid(), RESOLUTION, nullptr, options).overall,
                 MaptitudePins::QSCORE_CARBON_SIGMA08);
}

TEST(MetricCharacterizationTest, EdiamCarbonDefault) {
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    ExpectPinned(ediam(mol, ObsGrid(), RESOLUTION).overall, MaptitudePins::EDIAM_CARBON_DEFAULT);
}

TEST(MetricCharacterizationTest, CoverageCarbonDefault) {
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    ExpectPinned(coverage(mol, ObsGrid()).overall, MaptitudePins::COVERAGE_CARBON_DEFAULT);
}

TEST(MetricCharacterizationTest, CoverageCarbonSigma24) {
    // This grid's threshold passes the map maximum at sigma 18.79, so the only
    // sigma that pins the uncovered branch is one well above it. Paired with
    // COVERAGE_CARBON_DEFAULT, the two pins bracket the transition.
    CoverageOptions options;
    options.SetSigma(24.0);
    OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
    ExpectPinned(coverage(mol, ObsGrid(), nullptr, options).overall,
                 MaptitudePins::COVERAGE_CARBON_SIGMA24);
}
