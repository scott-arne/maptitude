/// Tier 1: trilinear interpolation (spec §3.2, §3.4, §3.5).
///
/// Domain and blend tests use a 4x4x4 unit grid whose node origin is the
/// Cartesian origin, so a coordinate is its fractional index. The equivalence
/// guard compares against OpenEye on a non-multilinear field. The timing harness
/// measures the before/after on the 1d26 asset.
#include <gtest/gtest.h>

#include <oechem.h>
#include <oegrid.h>

#include <chrono>
#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "maptitude/Grid.h"

using namespace Maptitude;

namespace {

constexpr unsigned int N = 4u;
constexpr double OUTSIDE = -99.0;

GridParams UnitGridParams() {
    return GridParams{0.0, 0.0, 0.0, N, N, N, 1.0, 1.0, 1.0};
}

/// values[el] == el, with el = iz*N*N + iy*N + ix.
std::vector<float> ElementNumberValues() {
    std::vector<float> values(N * N * N);
    for (unsigned int i = 0; i < values.size(); ++i) {
        values[i] = static_cast<float>(i);
    }
    return values;
}

float GuardFieldValue(const double x, const double y, const double z) {
    return static_cast<float>(x * x + 2.0 * y * y * y - 0.5 * z * z +
                              3.0 * x * y - y * z);
}

// The guard's whole purpose is to compare against OpenEye's scalar-grid
// interpolation, so this arm cannot move to the skew carrier.
OESystem::OEScalarGrid MakeGuardScalarGrid() {  // OE-SCALARGRID-OK: the guard's OpenEye arm
    double minmax[6] = {0.0, 0.0, 0.0, 3.0, 3.0, 3.0};
    OESystem::OEScalarGrid grid(minmax, 1.0);  // OE-SCALARGRID-OK: the guard's OpenEye arm
    for (unsigned int iz = 0; iz < 4u; ++iz) {
        for (unsigned int iy = 0; iy < 4u; ++iy) {
            for (unsigned int ix = 0; ix < 4u; ++ix) {
                grid[iz * 16u + iy * 4u + ix] = GuardFieldValue(ix, iy, iz);
            }
        }
    }
    return grid;
}

// Fractional index, then OEFloatGridLinearInterpolate's return at it.
// Captured in Task 2; the guard's reference after Task 3 removes the call.
struct GuardSample { double f[3]; double expected; };
static constexpr GuardSample GUARD_SAMPLES[] = {
    // Every entry is strictly interior:
    // 0 < f_i < 3 on all three axes, so the domain split does not apply.
    {{0.13, 0.13, 0.13}, 0.35879999399185181},
    {{0.13, 0.13, 0.5}, 0.12569998204708099},
    {{0.13, 0.13, 1}, -0.18930001556873322},
    {{0.13, 0.13, 1.47}, -0.95540004968643188},
    {{0.13, 0.13, 2}, -1.8193000555038452},
    {{0.13, 0.13, 2.8599999999999999}, -4.0810995101928711},
    {{0.13, 0.5, 0.13}, 1.1949999332427979},
    {{0.13, 0.5, 0.5}, 0.82499998807907104},
    {{0.13, 0.5, 1}, 0.32499998807907104},
    {{0.13, 0.5, 1.47}, -0.61500006914138794},
    {{0.13, 0.5, 2}, -1.6749999523162842},
    {{0.13, 0.5, 2.8599999999999999}, -4.2549996376037598},
    {{0.13, 1, 0.13}, 2.3250000476837158},
    {{0.13, 1, 0.5}, 1.7699999809265137},
    {{0.13, 1, 1}, 1.0199999809265137},
    {{0.13, 1, 1.47}, -0.15500009059906006},
    {{0.13, 1, 2}, -1.4800000190734863},
    {{0.13, 1, 2.8599999999999999}, -4.4899997711181641},
    {{0.13, 1.47, 0.13}, 9.0272006988525391},
    {{0.13, 1.47, 0.5}, 8.2983007431030273},
    {{0.13, 1.47, 1}, 7.3133001327514648},
    {{0.13, 1.47, 1.47}, 5.9174003601074219},
    {{0.13, 1.47, 2}, 4.3433003425598145},
    {{0.13, 1.47, 2.8599999999999999}, 0.92910069227218628},
    {{0.13, 2, 0.13}, 16.584999084472656},
    {{0.13, 2, 0.5}, 15.659999847412109},
    {{0.13, 2, 1}, 14.409999847412109},
    {{0.13, 2, 1.47}, 12.764999389648438},
    {{0.13, 2, 2}, 10.909999847412109},
    {{0.13, 2, 2.8599999999999999}, 7.0400004386901855},
    {{0.13, 2.8599999999999999, 0.13}, 49.488594055175781},
    {{0.13, 2.8599999999999999, 0.5}, 48.245395660400391},
    {{0.13, 2.8599999999999999, 1}, 46.565395355224609},
    {{0.13, 2.8599999999999999, 1.47}, 44.516197204589844},
    {{0.13, 2.8599999999999999, 2}, 42.205394744873047},
    {{0.13, 2.8599999999999999, 2.8599999999999999}, 37.595798492431641},
    {{0.5, 0.13, 0.13}, 0.87309998273849487},
    {{0.5, 0.13, 0.5}, 0.63999998569488525},
    {{0.5, 0.13, 1}, 0.32499998807907104},
    {{0.5, 0.13, 1.47}, -0.44110006093978882},
    {{0.5, 0.13, 2}, -1.3050000667572021},
    {{0.5, 0.13, 2.8599999999999999}, -3.5667996406555176},
    {{0.5, 0.5, 0.13}, 2.119999885559082},
    {{0.5, 0.5, 0.5}, 1.75},
    {{0.5, 0.5, 1}, 1.25},
    {{0.5, 0.5, 1.47}, 0.30999994277954102},
    {{0.5, 0.5, 2}, -0.75},
    {{0.5, 0.5, 2.8599999999999999}, -3.3299996852874756},
    {{0.5, 1, 0.13}, 3.8050000667572021},
    {{0.5, 1, 0.5}, 3.25},
    {{0.5, 1, 1}, 2.5},
    {{0.5, 1, 1.47}, 1.3249999284744263},
    {{0.5, 1, 2}, 0},
    {{0.5, 1, 2.8599999999999999}, -3.0099997520446777},
    {{0.5, 1.47, 0.13}, 11.028900146484375},
    {{0.5, 1.47, 0.5}, 10.300000190734863},
    {{0.5, 1.47, 1}, 9.3150005340576172},
    {{0.5, 1.47, 1.47}, 7.919100284576416},
    {{0.5, 1.47, 2}, 6.3450002670288086},
    {{0.5, 1.47, 2.8599999999999999}, 2.9308006763458252},
    {{0.5, 2, 0.13}, 19.174999237060547},
    {{0.5, 2, 0.5}, 18.25},
    {{0.5, 2, 1}, 17},
    {{0.5, 2, 1.47}, 15.354999542236328},
    {{0.5, 2, 2}, 13.5},
    {{0.5, 2, 2.8599999999999999}, 9.630000114440918},
    {{0.5, 2.8599999999999999, 0.13}, 53.033195495605469},
    {{0.5, 2.8599999999999999, 0.5}, 51.789997100830078},
    {{0.5, 2.8599999999999999, 1}, 50.109996795654297},
    {{0.5, 2.8599999999999999, 1.47}, 48.060794830322266},
    {{0.5, 2.8599999999999999, 2}, 45.749996185302734},
    {{0.5, 2.8599999999999999, 2.8599999999999999}, 41.140396118164062},
    {{1, 0.13, 0.13}, 1.5680999755859375},
    {{1, 0.13, 0.5}, 1.3350000381469727},
    {{1, 0.13, 1}, 1.0199999809265137},
    {{1, 0.13, 1.47}, 0.25389993190765381},
    {{1, 0.13, 2}, -0.61000001430511475},
    {{1, 0.13, 2.8599999999999999}, -2.8717997074127197},
    {{1, 0.5, 0.13}, 3.369999885559082},
    {{1, 0.5, 0.5}, 3},
    {{1, 0.5, 1}, 2.5},
    {{1, 0.5, 1.47}, 1.559999942779541},
    {{1, 0.5, 2}, 0.5},
    {{1, 0.5, 2.8599999999999999}, -2.0799996852874756},
    {{1, 1, 0.13}, 5.804999828338623},
    {{1, 1, 0.5}, 5.25},
    {{1, 1, 1}, 4.5},
    {{1, 1, 1.47}, 3.3249998092651367},
    {{1, 1, 2}, 2},
    {{1, 1, 2.8599999999999999}, -1.0099996328353882},
    {{1, 1.47, 0.13}, 13.73390007019043},
    {{1, 1.47, 0.5}, 13.005000114440918},
    {{1, 1.47, 1}, 12.020000457763672},
    {{1, 1.47, 1.47}, 10.624100685119629},
    {{1, 1.47, 2}, 9.0500001907348633},
    {{1, 1.47, 2.8599999999999999}, 5.635800838470459},
    {{1, 2, 0.13}, 22.674999237060547},
    {{1, 2, 0.5}, 21.75},
    {{1, 2, 1}, 20.5},
    {{1, 2, 1.47}, 18.854999542236328},
    {{1, 2, 2}, 17},
    {{1, 2, 2.8599999999999999}, 13.130000114440918},
    {{1, 2.8599999999999999, 0.13}, 57.823196411132812},
    {{1, 2.8599999999999999, 0.5}, 56.579994201660156},
    {{1, 2.8599999999999999, 1}, 54.899993896484375},
    {{1, 2.8599999999999999, 1.47}, 52.850795745849609},
    {{1, 2.8599999999999999, 2}, 50.539997100830078},
    {{1, 2.8599999999999999, 2.8599999999999999}, 45.930397033691406},
    {{1.47, 0.13, 0.13}, 3.1614000797271729},
    {{1.47, 0.13, 0.5}, 2.928300142288208},
    {{1.47, 0.13, 1}, 2.613300085067749},
    {{1.47, 0.13, 1.47}, 1.8472000360488892},
    {{1.47, 0.13, 2}, 0.98330008983612061},
    {{1.47, 0.13, 2.8599999999999999}, -1.2784996032714844},
    {{1.47, 0.5, 0.13}, 5.4850001335144043},
    {{1.47, 0.5, 0.5}, 5.1150002479553223},
    {{1.47, 0.5, 1}, 4.6150002479553223},
    {{1.47, 0.5, 1.47}, 3.6750001907348633},
    {{1.47, 0.5, 2}, 2.6150002479553223},
    {{1.47, 0.5, 2.8599999999999999}, 0.035000443458557129},
    {{1.47, 1, 0.13}, 8.625},
    {{1.47, 1, 0.5}, 8.0699996948242188},
    {{1.47, 1, 1}, 7.320000171661377},
    {{1.47, 1, 1.47}, 6.1449999809265137},
    {{1.47, 1, 2}, 4.820000171661377},
    {{1.47, 1, 2.8599999999999999}, 1.8100005388259888},
    {{1.47, 1.47, 0.13}, 17.21660041809082},
    {{1.47, 1.47, 0.5}, 16.487701416015625},
    {{1.47, 1.47, 1}, 15.502700805664062},
    {{1.47, 1.47, 1.47}, 14.10680103302002},
    {{1.47, 1.47, 2}, 12.532700538635254},
    {{1.47, 1.47, 2.8599999999999999}, 9.1185007095336914},
    {{1.47, 2, 0.13}, 26.905000686645508},
    {{1.47, 2, 0.5}, 25.979999542236328},
    {{1.47, 2, 1}, 24.729999542236328},
    {{1.47, 2, 1.47}, 23.085000991821289},
    {{1.47, 2, 2}, 21.229999542236328},
    {{1.47, 2, 2.8599999999999999}, 17.360000610351562},
    {{1.47, 2.8599999999999999, 0.13}, 63.265796661376953},
    {{1.47, 2.8599999999999999, 0.5}, 62.022594451904297},
    {{1.47, 2.8599999999999999, 1}, 60.342594146728516},
    {{1.47, 2.8599999999999999, 1.47}, 58.29339599609375},
    {{1.47, 2.8599999999999999, 2}, 55.982597351074219},
    {{1.47, 2.8599999999999999, 2.8599999999999999}, 51.372997283935547},
    {{2, 0.13, 0.13}, 4.9580998420715332},
    {{2, 0.13, 0.5}, 4.7249999046325684},
    {{2, 0.13, 1}, 4.4099998474121094},
    {{2, 0.13, 1.47}, 3.6438999176025391},
    {{2, 0.13, 2}, 2.7799999713897705},
    {{2, 0.13, 2.8599999999999999}, 0.51820027828216553},
    {{2, 0.5, 0.13}, 7.869999885559082},
    {{2, 0.5, 0.5}, 7.5},
    {{2, 0.5, 1}, 7},
    {{2, 0.5, 1.47}, 6.059999942779541},
    {{2, 0.5, 2}, 5},
    {{2, 0.5, 2.8599999999999999}, 2.4200003147125244},
    {{2, 1, 0.13}, 11.805000305175781},
    {{2, 1, 0.5}, 11.25},
    {{2, 1, 1}, 10.5},
    {{2, 1, 1.47}, 9.3249998092651367},
    {{2, 1, 2}, 8},
    {{2, 1, 2.8599999999999999}, 4.9900002479553223},
    {{2, 1.47, 0.13}, 21.143899917602539},
    {{2, 1.47, 0.5}, 20.415000915527344},
    {{2, 1.47, 1}, 19.430000305175781},
    {{2, 1.47, 1.47}, 18.034099578857422},
    {{2, 1.47, 2}, 16.460000991821289},
    {{2, 1.47, 2.8599999999999999}, 13.045801162719727},
    {{2, 2, 0.13}, 31.674999237060547},
    {{2, 2, 0.5}, 30.75},
    {{2, 2, 1}, 29.5},
    {{2, 2, 1.47}, 27.854999542236328},
    {{2, 2, 2}, 26},
    {{2, 2, 2.8599999999999999}, 22.130001068115234},
    {{2, 2.8599999999999999, 0.13}, 69.4031982421875},
    {{2, 2.8599999999999999, 0.5}, 68.159996032714844},
    {{2, 2.8599999999999999, 1}, 66.479995727539062},
    {{2, 2.8599999999999999, 1.47}, 64.430793762207031},
    {{2, 2.8599999999999999, 2}, 62.1199951171875},
    {{2, 2.8599999999999999, 2.8599999999999999}, 57.510395050048828},
    {{2.8599999999999999, 0.13, 0.13}, 9.5934991836547852},
    {{2.8599999999999999, 0.13, 0.5}, 9.3603992462158203},
    {{2.8599999999999999, 0.13, 1}, 9.0453996658325195},
    {{2.8599999999999999, 0.13, 1.47}, 8.2792997360229492},
    {{2.8599999999999999, 0.13, 2}, 7.4153995513916016},
    {{2.8599999999999999, 0.13, 2.8599999999999999}, 5.153599739074707},
    {{2.8599999999999999, 0.5, 0.13}, 13.459999084472656},
    {{2.8599999999999999, 0.5, 0.5}, 13.089999198913574},
    {{2.8599999999999999, 0.5, 1}, 12.589999198913574},
    {{2.8599999999999999, 0.5, 1.47}, 11.649999618530273},
    {{2.8599999999999999, 0.5, 2}, 10.589999198913574},
    {{2.8599999999999999, 0.5, 2.8599999999999999}, 8.0099992752075195},
    {{2.8599999999999999, 1, 0.13}, 18.684999465942383},
    {{2.8599999999999999, 1, 0.5}, 18.129999160766602},
    {{2.8599999999999999, 1, 1}, 17.379999160766602},
    {{2.8599999999999999, 1, 1.47}, 16.204999923706055},
    {{2.8599999999999999, 1, 2}, 14.879999160766602},
    {{2.8599999999999999, 1, 2.8599999999999999}, 11.869999885559082},
    {{2.8599999999999999, 1.47, 0.13}, 29.236499786376953},
    {{2.8599999999999999, 1.47, 0.5}, 28.507598876953125},
    {{2.8599999999999999, 1.47, 1}, 27.522600173950195},
    {{2.8599999999999999, 1.47, 1.47}, 26.126699447631836},
    {{2.8599999999999999, 1.47, 2}, 24.55259895324707},
    {{2.8599999999999999, 1.47, 2.8599999999999999}, 21.138399124145508},
    {{2.8599999999999999, 2, 0.13}, 41.134998321533203},
    {{2.8599999999999999, 2, 0.5}, 40.209999084472656},
    {{2.8599999999999999, 2, 1}, 38.959999084472656},
    {{2.8599999999999999, 2, 1.47}, 37.314998626708984},
    {{2.8599999999999999, 2, 2}, 35.459999084472656},
    {{2.8599999999999999, 2, 2.8599999999999999}, 31.590000152587891},
    {{2.8599999999999999, 2.8599999999999999, 0.13}, 81.081993103027344},
    {{2.8599999999999999, 2.8599999999999999, 0.5}, 79.838790893554688},
    {{2.8599999999999999, 2.8599999999999999, 1}, 78.158790588378906},
    {{2.8599999999999999, 2.8599999999999999, 1.47}, 76.109596252441406},
    {{2.8599999999999999, 2.8599999999999999, 2}, 73.798797607421875},
    {{2.8599999999999999, 2.8599999999999999, 2.8599999999999999}, 69.189193725585938},
};

}  // namespace

TEST(InterpolateDensityAt, ReturnsTheNodeValueAtANode) {
    const GridParams gp = UnitGridParams();
    const std::vector<float> v = ElementNumberValues();
    // (ix, iy, iz) = (1, 2, 3) -> el = 3*16 + 2*4 + 1 = 57.
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, 2.0, 3.0, OUTSIDE), 57.0);
}

TEST(InterpolateDensityAt, BlendsAlongEachAxisWithTheRightStride) {
    const GridParams gp = UnitGridParams();
    const std::vector<float> v = ElementNumberValues();
    // Element number is linear in the indices, so a half step along x adds 1/2,
    // along y adds 4/2, along z adds 16/2. A transposed stride moves all three.
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.5, 0.0, 0.0, OUTSIDE), 0.5);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, 0.5, 0.0, OUTSIDE), 2.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, 0.0, 0.5, OUTSIDE), 8.0);
}

TEST(InterpolateDensityAt, InterpolatesOnTheClosedFarFace) {
    const GridParams gp = UnitGridParams();
    const std::vector<float> v = ElementNumberValues();
    // f == n - 1 exactly. OEFloatGridLinearInterpolate returned the default
    // here; the closed interval is a deliberate departure that moves metrics.
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 3.0, 0.0, 0.0, OUTSIDE), 3.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 3.0, 3.0, 3.0, OUTSIDE), 63.0);
    // The face is closed, so a face-interior point bilinearly blends.
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 3.0, 0.5, 0.5, OUTSIDE), 13.0);
}

TEST(InterpolateDensityAt, ReturnsTheDefaultOneUlpOutsideEachFace) {
    const GridParams gp = UnitGridParams();
    const std::vector<float> v = ElementNumberValues();
    const double just_over = std::nextafter(3.0, 4.0);
    const double just_under = std::nextafter(0.0, -1.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), just_over, 0.0, 0.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, just_over, 0.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, 0.0, just_over, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), just_under, 0.0, 0.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, just_under, 0.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, 0.0, just_under, OUTSIDE), OUTSIDE);
}

TEST(InterpolateDensityAt, ReturnsTheDefaultInTheOuterHalfSpacingShell) {
    const GridParams gp = UnitGridParams();
    const std::vector<float> v = ElementNumberValues();
    // The scalar carrier's IsInGrid admitted this shell -- the old box ran to
    // 3.5 on each face. The node span does not, so the domain narrows here.
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 3.25, 1.0, 1.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), -0.25, 1.0, 1.0, OUTSIDE), OUTSIDE);
}

TEST(InterpolateDensityAt, ReturnsTheDefaultForANonFiniteCoordinate) {
    const GridParams gp = UnitGridParams();
    const std::vector<float> v = ElementNumberValues();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double inf = std::numeric_limits<double>::infinity();
    // A NaN fractional index fails no range comparison, so ContainsFractionalIndex
    // needs its isfinite conjuncts as well as its range ones. The infinities would
    // be caught by the range test alone; the three NaN cases would not.
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), nan, 1.0, 1.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, nan, 1.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, 1.0, nan, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), inf, 1.0, 1.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), -inf, 1.0, 1.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, inf, 1.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, -inf, 1.0, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, 1.0, inf, OUTSIDE), OUTSIDE);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, 1.0, -inf, OUTSIDE), OUTSIDE);
}

TEST(InterpolateDensityAt, MeasuresFromTheGridOriginNotFromZero) {
    // Every other case here sits on a zero origin, which cannot distinguish
    // (x - x_origin) / spacing from x / spacing. This one can: the whole grid
    // is translated, so world (0,0,0) falls outside it entirely.
    const GridParams gp{2.0, -3.0, 10.0, N, N, N, 1.0, 1.0, 1.0};
    const std::vector<float> v = ElementNumberValues();
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 2.0, -3.0, 10.0, OUTSIDE), 0.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 3.0, -3.0, 10.0, OUTSIDE), 1.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 2.0, -2.0, 10.0, OUTSIDE), 4.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 2.0, -3.0, 11.0, OUTSIDE), 16.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 2.5, -3.0, 10.0, OUTSIDE), 0.5);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, 0.0, 0.0, OUTSIDE), OUTSIDE);
}

TEST(InterpolateDensityAt, UsesThePerAxisSpacingIndependently) {
    const GridParams gp{0.0, 0.0, 0.0, N, N, N, 0.5, 1.0, 2.0};
    const std::vector<float> v = ElementNumberValues();
    // One full node step on each axis: 0.5 A in x, 1.0 A in y, 2.0 A in z.
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.5, 0.0, 0.0, OUTSIDE), 1.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, 1.0, 0.0, OUTSIDE), 4.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, 0.0, 2.0, OUTSIDE), 16.0);
    // A grid that collapsed the three spacings onto one would place these
    // elsewhere; 1.0 A along x is two node steps, not one.
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, 0.0, 0.0, OUTSIDE), 2.0);
}

TEST(InterpolateDensityAt, BlendsAcrossASingleCellOnATwoNodeAxis) {
    // n == 2 is the minimum get_grid_params permits, and the case the clamp's
    // `n_i - 2u` safety argument rests on: every point falls in the one cell.
    const GridParams gp{0.0, 0.0, 0.0, 2u, 2u, 2u, 1.0, 1.0, 1.0};
    std::vector<float> v(8);
    for (unsigned int i = 0; i < v.size(); ++i) {
        v[i] = static_cast<float>(i);
    }

    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.0, 0.0, 0.0, OUTSIDE), 0.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, 1.0, 1.0, OUTSIDE), 7.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.0, 0.0, 0.0, OUTSIDE), 1.0);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 0.5, 0.5, 0.5, OUTSIDE), 3.5);
    EXPECT_DOUBLE_EQ(interpolate_density_at(gp, v.data(), 1.5, 0.0, 0.0, OUTSIDE), OUTSIDE);
}

// Disabled by default: it reads a 222 KB asset and runs millions of samples.
// Run with:
//   ./build-release/tests/cpp/maptitude_tests \
//     --gtest_also_run_disabled_tests --gtest_filter=Interpolation*Timing*
TEST(InterpolationTiming, DISABLED_OpenEyeVersusMaptitudeOn1d26) {
#ifndef NDEBUG
    GTEST_SKIP() << "timing is only meaningful against an optimized build; "
                    "configure with `cmake --preset local-release` and run "
                    "./build-release/tests/cpp/maptitude_tests";
#endif

    const std::string path = std::string(MAPTITUDE_TEST_ASSET_DIR) + "/1d26_2fofc.ccp4";

    // The OpenEye arm of the timing comparison needs the scalar overload of
    // OEFloatGridLinearInterpolate, which is what this test measures against.
    OESystem::OEScalarGrid scalar;  // OE-SCALARGRID-OK: the OpenEye arm of the comparison
    ASSERT_TRUE(OESystem::OEReadGrid(path, scalar)) << "failed to read " << path;
    const OESystem::OESkewGrid skew(scalar);

    const GridParams gp = get_grid_params(skew);
    const float* values = skew.GetValues();
    ASSERT_NE(values, nullptr);

    // A deterministic sweep of the interior, avoiding the faces so both arms
    // sample the same domain.
    std::vector<double> pts;
    for (unsigned int iz = 1; iz + 1 < gp.z_dim; ++iz) {
        for (unsigned int iy = 1; iy + 1 < gp.y_dim; ++iy) {
            for (unsigned int ix = 1; ix + 1 < gp.x_dim; ++ix) {
                pts.push_back(gp.x_origin + (ix + 0.37) * gp.x_spacing);
                pts.push_back(gp.y_origin + (iy + 0.61) * gp.y_spacing);
                pts.push_back(gp.z_origin + (iz + 0.19) * gp.z_spacing);
            }
        }
    }
    const size_t n = pts.size() / 3;
    constexpr int REPEATS = 20;

    double sink = 0.0;
    auto t0 = std::chrono::steady_clock::now();
    for (int r = 0; r < REPEATS; ++r) {
        for (size_t i = 0; i < n; ++i) {
            sink += OESystem::OEFloatGridLinearInterpolate(
                scalar, static_cast<float>(pts[i * 3]),
                static_cast<float>(pts[i * 3 + 1]),
                static_cast<float>(pts[i * 3 + 2]), 0.0f);
        }
    }
    auto t1 = std::chrono::steady_clock::now();
    for (int r = 0; r < REPEATS; ++r) {
        for (size_t i = 0; i < n; ++i) {
            sink += interpolate_density_at(gp, values, pts[i * 3], pts[i * 3 + 1],
                                           pts[i * 3 + 2], 0.0);
        }
    }
    auto t2 = std::chrono::steady_clock::now();

    const double oe_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    const double mt_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
    std::cout << "points=" << n << " repeats=" << REPEATS
              << " OEFloatGridLinearInterpolate=" << oe_ms << " ms"
              << " interpolate_density_at=" << mt_ms << " ms"
              << " ratio=" << (oe_ms / mt_ms) << "\n";
    EXPECT_NE(sink, 0.0) << "sink keeps both loops from being optimized away";
}

// ---- Interpolation equivalence guard (spec §3.4) ----
//
// Fixture: a 4x4x4 isotropic grid with unit spacing whose node origin is the
// Cartesian origin, so a coordinate is its own fractional index.
//
// The value field is deliberately NOT multilinear. A trilinear blend
// reproduces a multilinear field exactly, so a guard built on one would agree
// with any index or weight error that happens to be self-consistent.

TEST(InterpolationEquivalence, MatchesOpenEyeOnTheSharedInteriorDomain) {
    const OESystem::OEScalarGrid scalar = MakeGuardScalarGrid();  // OE-SCALARGRID-OK: the guard's OpenEye arm
    const OESystem::OESkewGrid skew(scalar);
    const GridParams gp = get_grid_params(skew);
    const float* values = skew.GetValues();
    ASSERT_NE(values, nullptr);

    for (const GuardSample& s : GUARD_SAMPLES) {
        const float oe = OESystem::OEFloatGridLinearInterpolate(
            scalar, static_cast<float>(s.f[0]), static_cast<float>(s.f[1]),
            static_cast<float>(s.f[2]), 0.0f);
        EXPECT_NEAR(oe, s.expected, 1e-5)
            << "pinned reference drifted at (" << s.f[0] << ", " << s.f[1]
            << ", " << s.f[2] << ")";
        // OpenEye's double→float→double round-trip incurs representation error;
        // use relative tolerance for the maptitude arm.
        EXPECT_NEAR(interpolate_density_at(gp, values, s.f[0], s.f[1], s.f[2], 0.0),
                    s.expected, 1e-5 * std::max(1.0, std::abs(s.expected)))
            << "maptitude disagrees at (" << s.f[0] << ", " << s.f[1]
            << ", " << s.f[2] << ")";
    }
}

// ---- Far-face divergence (spec §3.4) ----
//
// OEFloatGridLinearInterpolate treats the far face as outside the grid and
// returns the caller's default. maptitude closes the interval and blends from
// the nodes that are there. Each test pins both halves; Task 3 keeps only the
// maptitude half, because by then there is no OpenEye call left to make.
//
// -99.0f is the sentinel default: no value in this field is near it, so a
// "returned the default" result cannot be confused with a blend.

TEST(InterpolationFarFace, DivergesFromOpenEyeAtTheFarCorner) {
    const OESystem::OEScalarGrid scalar = MakeGuardScalarGrid();  // OE-SCALARGRID-OK: the divergence test's OpenEye arm
    const OESystem::OESkewGrid skew(scalar);
    const GridParams gp = get_grid_params(skew);
    const float* values = skew.GetValues();
    ASSERT_NE(values, nullptr);

    EXPECT_FLOAT_EQ(
        OESystem::OEFloatGridLinearInterpolate(scalar, 3.0f, 3.0f, 3.0f, -99.0f),
        -99.0f);
    // All three fractional indices at n - 1: no blend, the corner node itself.
    EXPECT_NEAR(interpolate_density_at(gp, values, 3.0, 3.0, 3.0, -99.0),
                GuardFieldValue(3.0, 3.0, 3.0), 1e-4);
}

TEST(InterpolationFarFace, DivergesFromOpenEyeOnAFarFace) {
    const OESystem::OEScalarGrid scalar = MakeGuardScalarGrid();  // OE-SCALARGRID-OK: the divergence test's OpenEye arm
    const OESystem::OESkewGrid skew(scalar);
    const GridParams gp = get_grid_params(skew);
    const float* values = skew.GetValues();
    ASSERT_NE(values, nullptr);

    EXPECT_FLOAT_EQ(
        OESystem::OEFloatGridLinearInterpolate(scalar, 3.0f, 1.5f, 1.5f, -99.0f),
        -99.0f);
    // One index pinned at n - 1, two interior: a bilinear blend of four nodes.
    const double expected =
        0.25 * (GuardFieldValue(3.0, 1.0, 1.0) + GuardFieldValue(3.0, 2.0, 1.0) +
                GuardFieldValue(3.0, 1.0, 2.0) + GuardFieldValue(3.0, 2.0, 2.0));
    EXPECT_NEAR(interpolate_density_at(gp, values, 3.0, 1.5, 1.5, -99.0),
                expected, 1e-4);
}

TEST(InterpolationFarFace, DivergesFromOpenEyeOnAFarEdge) {
    const OESystem::OEScalarGrid scalar = MakeGuardScalarGrid();  // OE-SCALARGRID-OK: the divergence test's OpenEye arm
    const OESystem::OESkewGrid skew(scalar);
    const GridParams gp = get_grid_params(skew);
    const float* values = skew.GetValues();
    ASSERT_NE(values, nullptr);

    EXPECT_FLOAT_EQ(
        OESystem::OEFloatGridLinearInterpolate(scalar, 3.0f, 3.0f, 1.5f, -99.0f),
        -99.0f);
    // Two indices pinned at n - 1, one interior: a linear blend of two nodes.
    const double expected =
        0.5 * (GuardFieldValue(3.0, 3.0, 1.0) + GuardFieldValue(3.0, 3.0, 2.0));
    EXPECT_NEAR(interpolate_density_at(gp, values, 3.0, 3.0, 1.5, -99.0),
                expected, 1e-4);
}
