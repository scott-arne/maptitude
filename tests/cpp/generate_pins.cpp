/// Prints a complete pin_values.h to stdout.
///
/// Run before any behavior-changing fix; commit the output. The characterization
/// tests then fail loudly on any numeric drift. Regenerating this file is only
/// correct when a drift has been reviewed and accepted as one of the documented
/// exceptions in the Phase 1 spec.
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "maptitude/CoverageOptions.h"
#include "maptitude/DensityCalculator.h"
#include "maptitude/GridOps.h"
#include "maptitude/Metric.h"
#include "maptitude/QScoreOptions.h"
#include "maptitude/RsccOptions.h"
#include "maptitude/RsrOptions.h"
#include "maptitude/SymOp.h"
#include "maptitude/UnitCell.h"

#include "fixtures.h"

using namespace Maptitude;
using namespace MaptitudeTest;

namespace {

constexpr double HALF_WIDTH = 6.0;
constexpr double SPACING = 0.5;
constexpr double RESOLUTION = 2.0;

void Emit(std::ostream& os, const std::string& name, double value) {
    os << "constexpr double " << name << " = " << std::setprecision(17) << value << ";\n";
}

OESystem::OEScalarGrid ObsGrid() {
    return MakeGaussianGrid(0.0, 0.0, 0.0, 1.0, HALF_WIDTH, SPACING);
}

OESystem::OEScalarGrid CalcGrid() {
    return MakeGaussianGrid(0.25, -0.1, 0.15, 1.1, HALF_WIDTH, SPACING);
}

/// Two carbons: one on the observed grid's maximum, one 3.0 A out along x.
OEChem::OEGraphMol MakeTwoCarbonMol() {
    OEChem::OEGraphMol mol;
    const double coords[2][3] = {{0.0, 0.0, 0.0}, {3.0, 0.0, 0.0}};
    for (const auto& xyz : coords) {
        OEChem::OEAtomBase* atom = mol.NewAtom(6);
        mol.SetCoords(atom, xyz);
        OEChem::OEResidue residue;
        residue.SetName("LIG");
        residue.SetResidueNumber(1);
        residue.SetChainID('A');
        residue.SetBFactor(0.0);
        OEChem::OEAtomSetResidue(atom, residue);
    }
    return mol;
}

void EmitMetricPins(std::ostream& os) {
    const OESystem::OEScalarGrid obs = ObsGrid();
    const OESystem::OEScalarGrid calc = CalcGrid();

    {
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit(os, "RSCC_CARBON_BINNED", rscc(mol, obs, RESOLUTION, nullptr, &calc).overall);
    }
    {
        RsccOptions options;
        options.SetAtomRadiusMethod(AtomRadius::FIXED);
        options.SetFixedAtomRadius(1.5);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit(os, "RSCC_CARBON_FIXED", rscc(mol, obs, RESOLUTION, nullptr, &calc, options).overall);
    }
    {
        // Carbon's 1.7 A Bondi radius times the plan's original 1.2 scaling gave
        // 2.04 A, indistinguishable from the 2.00 A binned radius at 0.5 A spacing.
        // Scaling 1.5 reaches a strictly larger point set.
        RsccOptions options;
        options.SetAtomRadiusMethod(AtomRadius::SCALED);
        options.SetAtomRadiusScaling(1.5);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit(os, "RSCC_CARBON_SCALED", rscc(mol, obs, RESOLUTION, nullptr, &calc, options).overall);
    }
    {
        OEChem::OEGraphMol mol = MakeAtomMol(8, 0.3, 0.0, 0.0);
        Emit(os, "RSCC_OXYGEN_OFFSET", rscc(mol, obs, RESOLUTION, nullptr, &calc).overall);
    }
    {
        // RSCC's radius switch has no ADAPTIVE case, so ADAPTIVE falls through to the
        // binned default. Equal to RSCC_CARBON_BINNED by construction -- the equality
        // is the pin. Giving RSCC a real adaptive branch would move this and not that.
        RsccOptions options;
        options.SetAtomRadiusMethod(AtomRadius::ADAPTIVE);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit(os, "RSCC_CARBON_ADAPTIVE", rscc(mol, obs, RESOLUTION, nullptr, &calc, options).overall);
    }
    {
        // RsrOptions defaults to ADAPTIVE, so the binned path must be selected
        // explicitly -- otherwise this pin silently duplicates RSR_CARBON_ADAPTIVE.
        RsrOptions options;
        options.SetAtomRadiusMethod(AtomRadius::BINNED);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit(os, "RSR_CARBON_BINNED", rsr(mol, obs, RESOLUTION, nullptr, &calc, options).overall);
    }
    {
        RsrOptions options;
        options.SetAtomRadiusMethod(AtomRadius::ADAPTIVE);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit(os, "RSR_CARBON_ADAPTIVE", rsr(mol, obs, RESOLUTION, nullptr, &calc, options).overall);
    }
    {
        // Pins that RsrOptions still defaults to ADAPTIVE. Equal to
        // RSR_CARBON_ADAPTIVE by construction -- that equality IS the assertion, and
        // changing the default moves this pin while leaving the explicit one alone.
        RsrOptions options;
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit(os, "RSR_CARBON_DEFAULT", rsr(mol, obs, RESOLUTION, nullptr, &calc, options).overall);
    }
    {
        RsrOptions options;
        options.SetAtomRadiusMethod(AtomRadius::FIXED);
        options.SetFixedAtomRadius(1.5);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit(os, "RSR_CARBON_FIXED", rsr(mol, obs, RESOLUTION, nullptr, &calc, options).overall);
    }
    {
        RsrOptions options;
        options.SetAtomRadiusMethod(AtomRadius::SCALED);
        options.SetAtomRadiusScaling(1.5);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit(os, "RSR_CARBON_SCALED", rsr(mol, obs, RESOLUTION, nullptr, &calc, options).overall);
    }
    {
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit(os, "QSCORE_CARBON_DEFAULT", qscore(mol, obs, RESOLUTION).overall);
    }
    {
        QScoreOptions options;
        options.SetSigma(0.8);
        options.SetNumPoints(16);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit(os, "QSCORE_CARBON_SIGMA08", qscore(mol, obs, RESOLUTION, nullptr, options).overall);
    }
    {
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit(os, "EDIAM_CARBON_DEFAULT", ediam(mol, obs, RESOLUTION).overall);
    }
    {
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit(os, "COVERAGE_CARBON_DEFAULT", coverage(mol, obs).overall);
    }
    {
        // This grid's threshold passes the map maximum at sigma 18.79, so the only
        // sigma that pins the uncovered branch is one well above it. Paired with
        // COVERAGE_CARBON_DEFAULT, the two pins bracket the transition.
        CoverageOptions options;
        options.SetSigma(24.0);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit(os, "COVERAGE_CARBON_SIGMA24", coverage(mol, obs, nullptr, options).overall);
    }
    {
        // The only pin whose value no single-atom molecule can produce: one atom above
        // the threshold and one below, so 0.5 is the aggregation mean itself. Any
        // implementation that drops an atom or fails to average returns 1.0 or 0.0.
        CoverageOptions options;
        options.SetSigma(4.0);
        OEChem::OEGraphMol mol = MakeTwoCarbonMol();
        Emit(os, "COVERAGE_TWO_ATOM_SPLIT", coverage(mol, obs, nullptr, options).overall);
    }
}

/// Reduce a grid to five scalars: sum, sum-of-squares, min, max, and a
/// mean-centred order-sensitive index moment.
///
/// The index moment is sum((i+1) * (v[i] - mean)). Mean-centring removes the
/// permutation-invariant component: the mean contributes uniformly to every
/// index, so including it only inflates the tolerance without adding signal.
/// It catches deterministic axis-order and stride permutations (e.g., a
/// transposed grid with the same value multiset).
void EmitGridSummary(std::ostream& os, const std::string& prefix, const OESystem::OEScalarGrid& grid) {
    double sum = 0.0, sum_sq = 0.0, index_moment = 0.0;
    double lo = grid[0], hi = grid[0];

    // First pass: sum, sum_sq, min, max
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        const double v = grid[i];
        sum += v;
        sum_sq += v * v;
        lo = std::min(lo, v);
        hi = std::max(hi, v);
    }

    // Second pass: mean-centred index moment
    const double mean = sum / static_cast<double>(grid.GetSize());
    for (unsigned int i = 0; i < grid.GetSize(); ++i) {
        index_moment += static_cast<double>(i + 1) * (grid[i] - mean);
    }

    Emit(os, prefix + "_SUM", sum);
    Emit(os, prefix + "_SUM_SQ", sum_sq);
    Emit(os, prefix + "_MIN", lo);
    Emit(os, prefix + "_MAX", hi);
    Emit(os, prefix + "_INDEX_MOMENT", index_moment);
}

void EmitFcPins(std::ostream& os) {
    // P1 in an orthorhombic cell: the only geometry Phase 1 keeps supporting.
    UnitCell ortho(20.0, 25.0, 30.0, 90.0, 90.0, 90.0);
    std::vector<SymOp> symops = SymOp::ParseAll("x,y,z");

    OEChem::OEGraphMol mol = MakeAtomMol(6, 5.0, 5.0, 5.0);
    DensityCalculator calc(ortho, symops);
    OESystem::OEScalarGrid obs = MakeGaussianGrid(5.0, 5.0, 5.0, 1.0, 6.0, 0.5);
    std::unique_ptr<OESystem::OEScalarGrid> fc(calc.Calculate(mol, obs, 2.0));
    EmitGridSummary(os, "FC_ORTHORHOMBIC", *fc);

    // Monoclinic: produces numbers today. Task 9 makes this throw CellError.
    // The pin exists so that regression is a reviewed change, not a silent one.
    UnitCell mono(20.0, 25.0, 30.0, 90.0, 105.0, 90.0);
    OEChem::OEGraphMol mol2 = MakeAtomMol(6, 5.0, 5.0, 5.0);
    DensityCalculator calc2(mono, symops);
    std::unique_ptr<OESystem::OEScalarGrid> fc2(calc2.Calculate(mol2, obs, 2.0));
    EmitGridSummary(os, "FC_MONOCLINIC", *fc2);

    // Orthorhombic with n_scale_shells = 4, to cover the per-shell FFT scaling branch.
    OEChem::OEGraphMol mol3 = MakeAtomMol(6, 5.0, 5.0, 5.0);
    DensityCalculator calc3(ortho, symops);
    std::unique_ptr<OESystem::OEScalarGrid> fc3(
        calc3.Calculate(mol3, obs, RESOLUTION, nullptr, 0.35, 46.0, false, 4));
    EmitGridSummary(os, "FC_ORTHORHOMBIC_SHELLS4", *fc3);

    // Asymmetric atom position. The three cases above all place the atom at
    // (5,5,5) inside a cubic grid with an isotropic Gaussian, a configuration an
    // axis permutation maps almost onto itself -- their index moments move by
    // less than their own tolerance under an x/y swap and cannot detect a
    // stride or axis-order regression in the final interpolation loop. Distinct
    // coordinates break that symmetry: the weakest permutation moves this
    // moment by ~10^5 times its tolerance. Orthorhombic on purpose, so Task 9
    // leaves it in place as the surviving permutation guard.
    OEChem::OEGraphMol mol4 = MakeAtomMol(6, 5.0, 2.0, -1.0);
    DensityCalculator calc4(ortho, symops);
    OESystem::OEScalarGrid obs4 = MakeGaussianGrid(5.0, 2.0, -1.0, 1.0, HALF_WIDTH, SPACING);
    std::unique_ptr<OESystem::OEScalarGrid> fc4(calc4.Calculate(mol4, obs4, RESOLUTION));
    EmitGridSummary(os, "FC_ORTHORHOMBIC_ASYM", *fc4);
}

void EmitGridOpsPins(std::ostream& os) {
    // scale_map
    {
        OESystem::OEScalarGrid grid = MakeRampGrid(2.0, 0.5);
        scale_map(grid, 1.5);
        EmitGridSummary(os, "GRIDOPS_SCALE", grid);
    }

    // combine_maps: ADD
    {
        const OESystem::OEScalarGrid lhs = MakeRampGrid(2.0, 0.5);
        const OESystem::OEScalarGrid rhs = MakeUniformGrid(10.0f, 2.0, 0.5);
        std::unique_ptr<OESystem::OEScalarGrid> result(combine_maps(lhs, rhs, MapOp::ADD));
        EmitGridSummary(os, "GRIDOPS_ADD", *result);
    }

    // combine_maps: SUBTRACT
    {
        const OESystem::OEScalarGrid lhs = MakeRampGrid(2.0, 0.5);
        const OESystem::OEScalarGrid rhs = MakeUniformGrid(10.0f, 2.0, 0.5);
        std::unique_ptr<OESystem::OEScalarGrid> result(combine_maps(lhs, rhs, MapOp::SUBTRACT));
        EmitGridSummary(os, "GRIDOPS_SUBTRACT", *result);
    }

    // combine_maps: MIN
    {
        const OESystem::OEScalarGrid lhs = MakeRampGrid(2.0, 0.5);
        const OESystem::OEScalarGrid rhs = MakeUniformGrid(50.0f, 2.0, 0.5);
        std::unique_ptr<OESystem::OEScalarGrid> result(combine_maps(lhs, rhs, MapOp::MIN));
        EmitGridSummary(os, "GRIDOPS_MIN", *result);
    }

    // combine_maps: MAX
    {
        const OESystem::OEScalarGrid lhs = MakeRampGrid(2.0, 0.5);
        const OESystem::OEScalarGrid rhs = MakeUniformGrid(-50.0f, 2.0, 0.5);
        std::unique_ptr<OESystem::OEScalarGrid> result(combine_maps(lhs, rhs, MapOp::MAX));
        EmitGridSummary(os, "GRIDOPS_MAX", *result);
    }

    // diff_to_calc
    {
        const OESystem::OEScalarGrid obs = MakeRampGrid(2.0, 0.5);
        const OESystem::OEScalarGrid diff = MakeUniformGrid(3.0f, 2.0, 0.5);
        std::unique_ptr<OESystem::OEScalarGrid> result(diff_to_calc(obs, diff));
        EmitGridSummary(os, "GRIDOPS_DIFF_TO_CALC", *result);
    }
}

}  // namespace

int main() {
    // All-or-nothing generation: build the complete header in memory, then write
    // to stdout only if every emitter succeeds. On exception, the shell redirect
    // still empties the file, but the result is unambiguously empty and breaks
    // the build immediately, rather than a syntactically plausible truncation.
    // git checkout restores it.
    try {
        std::ostringstream buffer;
        buffer << "// GENERATED by tests/cpp/generate_pins.cpp. Do not hand-edit.\n";
        buffer << "//\n";
        buffer << "// Regenerate only when a numeric change has been reviewed and accepted:\n";
        buffer << "//   cmake --build build-debug --target maptitude_pin_generator\n";
        buffer << "//   ./build-debug/tests/cpp/maptitude_pin_generator > tests/cpp/pin_values.h\n";
        buffer << "#ifndef MAPTITUDE_TEST_PIN_VALUES_H\n";
        buffer << "#define MAPTITUDE_TEST_PIN_VALUES_H\n\n";
        buffer << "namespace MaptitudePins {\n\n";
        EmitMetricPins(buffer);
        EmitFcPins(buffer);
        EmitGridOpsPins(buffer);
        buffer << "\n}  // namespace MaptitudePins\n\n";
        buffer << "#endif  // MAPTITUDE_TEST_PIN_VALUES_H\n";

        std::cout << buffer.str();
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error generating pins: " << e.what() << "\n";
        return 1;
    }
}
