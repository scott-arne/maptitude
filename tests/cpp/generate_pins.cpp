/// Prints a complete pin_values.h to stdout.
///
/// Run before any behavior-changing fix; commit the output. The characterization
/// tests then fail loudly on any numeric drift. Regenerating this file is only
/// correct when a drift has been reviewed and accepted as one of the documented
/// exceptions in the Phase 1 spec.
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "maptitude/CoverageOptions.h"
#include "maptitude/Metric.h"
#include "maptitude/QScoreOptions.h"
#include "maptitude/RsccOptions.h"
#include "maptitude/RsrOptions.h"

#include "fixtures.h"

using namespace Maptitude;
using namespace MaptitudeTest;

namespace {

constexpr double HALF_WIDTH = 6.0;
constexpr double SPACING = 0.5;
constexpr double RESOLUTION = 2.0;

void Emit(const std::string& name, double value) {
    std::cout << "constexpr double " << name << " = " << std::setprecision(17) << value << ";\n";
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

void EmitMetricPins() {
    const OESystem::OEScalarGrid obs = ObsGrid();
    const OESystem::OEScalarGrid calc = CalcGrid();

    {
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit("RSCC_CARBON_BINNED", rscc(mol, obs, RESOLUTION, nullptr, &calc).overall);
    }
    {
        RsccOptions options;
        options.SetAtomRadiusMethod(AtomRadius::FIXED);
        options.SetFixedAtomRadius(1.5);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit("RSCC_CARBON_FIXED", rscc(mol, obs, RESOLUTION, nullptr, &calc, options).overall);
    }
    {
        // Carbon's 1.7 A Bondi radius times the plan's original 1.2 scaling gave
        // 2.04 A, indistinguishable from the 2.00 A binned radius at 0.5 A spacing.
        // Scaling 1.5 reaches a strictly larger point set.
        RsccOptions options;
        options.SetAtomRadiusMethod(AtomRadius::SCALED);
        options.SetAtomRadiusScaling(1.5);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit("RSCC_CARBON_SCALED", rscc(mol, obs, RESOLUTION, nullptr, &calc, options).overall);
    }
    {
        OEChem::OEGraphMol mol = MakeAtomMol(8, 0.3, 0.0, 0.0);
        Emit("RSCC_OXYGEN_OFFSET", rscc(mol, obs, RESOLUTION, nullptr, &calc).overall);
    }
    {
        // RSCC's radius switch has no ADAPTIVE case, so ADAPTIVE falls through to the
        // binned default. Equal to RSCC_CARBON_BINNED by construction -- the equality
        // is the pin. Giving RSCC a real adaptive branch would move this and not that.
        RsccOptions options;
        options.SetAtomRadiusMethod(AtomRadius::ADAPTIVE);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit("RSCC_CARBON_ADAPTIVE", rscc(mol, obs, RESOLUTION, nullptr, &calc, options).overall);
    }
    {
        // RsrOptions defaults to ADAPTIVE, so the binned path must be selected
        // explicitly -- otherwise this pin silently duplicates RSR_CARBON_ADAPTIVE.
        RsrOptions options;
        options.SetAtomRadiusMethod(AtomRadius::BINNED);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit("RSR_CARBON_BINNED", rsr(mol, obs, RESOLUTION, nullptr, &calc, options).overall);
    }
    {
        RsrOptions options;
        options.SetAtomRadiusMethod(AtomRadius::ADAPTIVE);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit("RSR_CARBON_ADAPTIVE", rsr(mol, obs, RESOLUTION, nullptr, &calc, options).overall);
    }
    {
        // Pins that RsrOptions still defaults to ADAPTIVE. Equal to
        // RSR_CARBON_ADAPTIVE by construction -- that equality IS the assertion, and
        // changing the default moves this pin while leaving the explicit one alone.
        RsrOptions options;
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit("RSR_CARBON_DEFAULT", rsr(mol, obs, RESOLUTION, nullptr, &calc, options).overall);
    }
    {
        RsrOptions options;
        options.SetAtomRadiusMethod(AtomRadius::FIXED);
        options.SetFixedAtomRadius(1.5);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit("RSR_CARBON_FIXED", rsr(mol, obs, RESOLUTION, nullptr, &calc, options).overall);
    }
    {
        RsrOptions options;
        options.SetAtomRadiusMethod(AtomRadius::SCALED);
        options.SetAtomRadiusScaling(1.5);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit("RSR_CARBON_SCALED", rsr(mol, obs, RESOLUTION, nullptr, &calc, options).overall);
    }
    {
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit("QSCORE_CARBON_DEFAULT", qscore(mol, obs, RESOLUTION).overall);
    }
    {
        QScoreOptions options;
        options.SetSigma(0.8);
        options.SetNumPoints(16);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit("QSCORE_CARBON_SIGMA08", qscore(mol, obs, RESOLUTION, nullptr, options).overall);
    }
    {
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit("EDIAM_CARBON_DEFAULT", ediam(mol, obs, RESOLUTION).overall);
    }
    {
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit("COVERAGE_CARBON_DEFAULT", coverage(mol, obs).overall);
    }
    {
        // This grid's threshold passes the map maximum at sigma 18.79, so the only
        // sigma that pins the uncovered branch is one well above it. Paired with
        // COVERAGE_CARBON_DEFAULT, the two pins bracket the transition.
        CoverageOptions options;
        options.SetSigma(24.0);
        OEChem::OEGraphMol mol = MakeAtomMol(6, 0.0, 0.0, 0.0);
        Emit("COVERAGE_CARBON_SIGMA24", coverage(mol, obs, nullptr, options).overall);
    }
    {
        // The only pin whose value no single-atom molecule can produce: one atom above
        // the threshold and one below, so 0.5 is the aggregation mean itself. Any
        // implementation that drops an atom or fails to average returns 1.0 or 0.0.
        CoverageOptions options;
        options.SetSigma(4.0);
        OEChem::OEGraphMol mol = MakeTwoCarbonMol();
        Emit("COVERAGE_TWO_ATOM_SPLIT", coverage(mol, obs, nullptr, options).overall);
    }
}

}  // namespace

int main() {
    std::cout << "// GENERATED by tests/cpp/generate_pins.cpp. Do not hand-edit.\n";
    std::cout << "//\n";
    std::cout << "// Regenerate only when a numeric change has been reviewed and accepted:\n";
    std::cout << "//   cmake --build build-debug --target maptitude_pin_generator\n";
    std::cout << "//   ./build-debug/tests/cpp/maptitude_pin_generator > tests/cpp/pin_values.h\n";
    std::cout << "#ifndef MAPTITUDE_TEST_PIN_VALUES_H\n";
    std::cout << "#define MAPTITUDE_TEST_PIN_VALUES_H\n\n";
    std::cout << "namespace MaptitudePins {\n\n";
    EmitMetricPins();
    std::cout << "\n}  // namespace MaptitudePins\n\n";
    std::cout << "#endif  // MAPTITUDE_TEST_PIN_VALUES_H\n";
    return 0;
}
