// swig/maptitude.i
// SWIG interface file for Maptitude Python bindings
%module _maptitude

%{
// Include all necessary headers
#include "maptitude/maptitude.h"
#include "maptitude/Error.h"
#include "maptitude/Residue.h"
#include "maptitude/UnitCell.h"
#include "maptitude/SymOp.h"
#include "maptitude/ScatteringFactors.h"
#include "maptitude/Grid.h"
#include "maptitude/DensityScoreResult.h"
#include "maptitude/QScoreOptions.h"
#include "maptitude/RSCCOptions.h"
#include "maptitude/DensityCalculator.h"
#include "maptitude/Metric.h"
#include "maptitude/GridOps.h"
#include "maptitude/SpatialIndex.h"

#include <oechem.h>
#include <oegrid.h>

using namespace Maptitude;
%}

// ============================================================================
// Forward declarations for cross-module SWIG type resolution
// ============================================================================
namespace OEChem {
    class OEMolBase;
    class OEAtomBase;
}

namespace OESystem {
    class OEScalarGrid;
}

namespace OESystem {
    template <class T> class OEUnaryPredicate;
}

// ============================================================================
// OEMolBase / OEAtomBase / OEScalarGrid typemaps for cross-module SWIG
// type resolution
// ============================================================================
// OpenEye's Python bindings register molecule objects as "OEGraphMolWrapper *"
// and "OEMolWrapper *" in the SWIG runtime type table. Our module registers
// "OEChem::OEMolBase *". SWIG's automatic cross-module linking doesn't resolve
// these different type names, so we use runtime type queries to bridge them.
// OpenEye uses SWIG runtime v4; our module uses v5. Since the runtimes are
// separate, SWIG_TypeQuery and SWIG_Python_GetSwigThis cannot access OpenEye
// types. We use Python isinstance for type safety and directly extract the
// void* pointer from the SwigPyObject struct layout (stable across versions).

%{
// Minimal SwigPyObject layout compatible across SWIG runtime versions.
// The actual struct may have more fields, but ptr is always first after
// PyObject_HEAD.
struct _SwigPyObjectCompat {
    PyObject_HEAD
    void *ptr;
};

// Cached type references for OEChem types
static PyObject* _maptitude_oe_molbase_type = NULL;
static PyObject* _maptitude_oe_atombase_type = NULL;
static PyObject* _maptitude_oe_scalargrid_type = NULL;

static bool _maptitude_is_oemolbase(PyObject* obj) {
    if (!_maptitude_oe_molbase_type) {
        PyObject* mod = PyImport_ImportModule("openeye.oechem");
        if (mod) {
            _maptitude_oe_molbase_type = PyObject_GetAttrString(mod, "OEMolBase");
            Py_DECREF(mod);
        }
        if (!_maptitude_oe_molbase_type) return false;
    }
    return PyObject_IsInstance(obj, _maptitude_oe_molbase_type) == 1;
}

static bool _maptitude_is_oeatombase(PyObject* obj) {
    if (!_maptitude_oe_atombase_type) {
        PyObject* mod = PyImport_ImportModule("openeye.oechem");
        if (mod) {
            _maptitude_oe_atombase_type = PyObject_GetAttrString(mod, "OEAtomBase");
            Py_DECREF(mod);
        }
        if (!_maptitude_oe_atombase_type) return false;
    }
    return PyObject_IsInstance(obj, _maptitude_oe_atombase_type) == 1;
}

static bool _maptitude_is_oescalargrid(PyObject* obj) {
    if (!_maptitude_oe_scalargrid_type) {
        PyObject* mod = PyImport_ImportModule("openeye.oegrid");
        if (mod) {
            _maptitude_oe_scalargrid_type = PyObject_GetAttrString(mod, "OEScalarGrid");
            Py_DECREF(mod);
        }
        if (!_maptitude_oe_scalargrid_type) return false;
    }
    return PyObject_IsInstance(obj, _maptitude_oe_scalargrid_type) == 1;
}

static void* _maptitude_extract_swig_ptr(PyObject* obj) {
    // Get the .this attribute which is a SwigPyObject from any SWIG version
    PyObject* thisAttr = PyObject_GetAttrString(obj, "this");
    if (!thisAttr) {
        PyErr_Clear();
        return NULL;
    }
    // Extract the void* ptr from the SwigPyObject-compatible layout
    void* ptr = ((_SwigPyObjectCompat*)thisAttr)->ptr;
    Py_DECREF(thisAttr);
    return ptr;
}

// Wrap a heap-allocated OEScalarGrid* as an OpenEye-native Python object.
// Creates an oegrid.OEScalarGrid wrapper, replaces its internal pointer with
// the provided grid, and deletes the default-constructed grid.  When Python
// garbage-collects the returned object, OpenEye's SWIG destructor calls
// delete on the pointer, so there is no memory leak and no deep copy.
static PyObject* _maptitude_wrap_as_oe_grid(OESystem::OEScalarGrid* grid) {
    if (!grid) {
        Py_RETURN_NONE;
    }
    PyObject* oegrid_mod = PyImport_ImportModule("openeye.oegrid");
    if (!oegrid_mod) {
        delete grid;
        return NULL;
    }
    PyObject* grid_cls = PyObject_GetAttrString(oegrid_mod, "OEScalarGrid");
    Py_DECREF(oegrid_mod);
    if (!grid_cls) {
        delete grid;
        return NULL;
    }
    // Create a default-constructed OEScalarGrid Python wrapper
    PyObject* oe_grid = PyObject_CallNoArgs(grid_cls);
    Py_DECREF(grid_cls);
    if (!oe_grid) {
        delete grid;
        return NULL;
    }
    // Access the SwigPyObject (.this) inside the OE wrapper
    PyObject* thisAttr = PyObject_GetAttrString(oe_grid, "this");
    if (!thisAttr) {
        PyErr_Clear();
        Py_DECREF(oe_grid);
        delete grid;
        return NULL;
    }
    // Swap: delete the default grid, install ours
    _SwigPyObjectCompat* swig_this = (_SwigPyObjectCompat*)thisAttr;
    delete reinterpret_cast<OESystem::OEScalarGrid*>(swig_this->ptr);
    swig_this->ptr = grid;
    Py_DECREF(thisAttr);
    return oe_grid;
}
%}

// ============================================================================
// Typemap: non-const OEMolBase& (for functions that modify the molecule)
// ============================================================================
%typemap(in) OEChem::OEMolBase& (void *argp = 0, int res = 0) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
    if (!SWIG_IsOK(res)) {
        if (_maptitude_is_oemolbase($input)) {
            argp = _maptitude_extract_swig_ptr($input);
            if (argp) res = SWIG_OK;
        }
    }
    if (!SWIG_IsOK(res)) {
        SWIG_exception_fail(SWIG_ArgError(res), "Expected OEMolBase-derived object. Ensure openeye.oechem is imported.");
    }
    if (!argp) {
        SWIG_exception_fail(SWIG_NullReferenceError, "Null OEMolBase reference.");
    }
    $1 = reinterpret_cast< $1_ltype >(argp);
}

%typemap(typecheck, precedence=10) OEChem::OEMolBase& {
    void *vptr = 0;
    int res = SWIG_ConvertPtr($input, &vptr, $descriptor, SWIG_POINTER_NO_NULL);
    $1 = SWIG_IsOK(res) ? 1 : _maptitude_is_oemolbase($input) ? 1 : 0;
}

// ============================================================================
// Typemap: const OEMolBase& (for read-only molecule access)
// ============================================================================
%typemap(in) const OEChem::OEMolBase& (void *argp = 0, int res = 0) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
    if (!SWIG_IsOK(res)) {
        if (_maptitude_is_oemolbase($input)) {
            argp = _maptitude_extract_swig_ptr($input);
            if (argp) res = SWIG_OK;
        }
    }
    if (!SWIG_IsOK(res)) {
        SWIG_exception_fail(SWIG_ArgError(res), "Expected OEMolBase-derived object. Ensure openeye.oechem is imported.");
    }
    if (!argp) {
        SWIG_exception_fail(SWIG_NullReferenceError, "Null OEMolBase reference.");
    }
    $1 = reinterpret_cast< $1_ltype >(argp);
}

%typemap(typecheck, precedence=10) const OEChem::OEMolBase& {
    void *vptr = 0;
    int res = SWIG_ConvertPtr($input, &vptr, $descriptor, SWIG_POINTER_NO_NULL);
    $1 = SWIG_IsOK(res) ? 1 : _maptitude_is_oemolbase($input) ? 1 : 0;
}

// ============================================================================
// Typemap: const OEAtomBase& (for per-atom operations)
// ============================================================================
%typemap(in) const OEChem::OEAtomBase& (void *argp = 0, int res = 0) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
    if (!SWIG_IsOK(res)) {
        if (_maptitude_is_oeatombase($input)) {
            argp = _maptitude_extract_swig_ptr($input);
            if (argp) res = SWIG_OK;
        }
    }
    if (!SWIG_IsOK(res)) {
        SWIG_exception_fail(SWIG_ArgError(res), "Expected OEAtomBase-derived object. Ensure openeye.oechem is imported.");
    }
    if (!argp) {
        SWIG_exception_fail(SWIG_NullReferenceError, "Null OEAtomBase reference.");
    }
    $1 = reinterpret_cast< $1_ltype >(argp);
}

%typemap(typecheck, precedence=10) const OEChem::OEAtomBase& {
    void *vptr = 0;
    int res = SWIG_ConvertPtr($input, &vptr, $descriptor, SWIG_POINTER_NO_NULL);
    $1 = SWIG_IsOK(res) ? 1 : _maptitude_is_oeatombase($input) ? 1 : 0;
}

// ============================================================================
// Typemap: OEScalarGrid& (non-const, for functions that modify the grid)
// ============================================================================
%typemap(in) OESystem::OEScalarGrid& (void *argp = 0, int res = 0) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
    if (!SWIG_IsOK(res)) {
        if (_maptitude_is_oescalargrid($input)) {
            argp = _maptitude_extract_swig_ptr($input);
            if (argp) res = SWIG_OK;
        }
    }
    if (!SWIG_IsOK(res)) {
        SWIG_exception_fail(SWIG_ArgError(res), "Expected OEScalarGrid-derived object. Ensure openeye.oegrid is imported.");
    }
    if (!argp) {
        SWIG_exception_fail(SWIG_NullReferenceError, "Null OEScalarGrid reference.");
    }
    $1 = reinterpret_cast< $1_ltype >(argp);
}

%typemap(typecheck, precedence=10) OESystem::OEScalarGrid& {
    void *vptr = 0;
    int res = SWIG_ConvertPtr($input, &vptr, $descriptor, SWIG_POINTER_NO_NULL);
    $1 = SWIG_IsOK(res) ? 1 : _maptitude_is_oescalargrid($input) ? 1 : 0;
}

// ============================================================================
// Typemap: const OEScalarGrid& (for read-only grid access)
// ============================================================================
%typemap(in) const OESystem::OEScalarGrid& (void *argp = 0, int res = 0) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
    if (!SWIG_IsOK(res)) {
        if (_maptitude_is_oescalargrid($input)) {
            argp = _maptitude_extract_swig_ptr($input);
            if (argp) res = SWIG_OK;
        }
    }
    if (!SWIG_IsOK(res)) {
        SWIG_exception_fail(SWIG_ArgError(res), "Expected OEScalarGrid-derived object. Ensure openeye.oegrid is imported.");
    }
    if (!argp) {
        SWIG_exception_fail(SWIG_NullReferenceError, "Null OEScalarGrid reference.");
    }
    $1 = reinterpret_cast< $1_ltype >(argp);
}

%typemap(typecheck, precedence=10) const OESystem::OEScalarGrid& {
    void *vptr = 0;
    int res = SWIG_ConvertPtr($input, &vptr, $descriptor, SWIG_POINTER_NO_NULL);
    $1 = SWIG_IsOK(res) ? 1 : _maptitude_is_oescalargrid($input) ? 1 : 0;
}

// ============================================================================
// Typemap: const OEScalarGrid* (optional pointer, e.g. calc_grid)
// ============================================================================
%typemap(in) const OESystem::OEScalarGrid* (void *argp = 0, int res = 0) {
    if ($input == Py_None) {
        $1 = NULL;
    } else {
        res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
        if (!SWIG_IsOK(res)) {
            if (_maptitude_is_oescalargrid($input)) {
                argp = _maptitude_extract_swig_ptr($input);
                if (argp) res = SWIG_OK;
            }
        }
        if (!SWIG_IsOK(res)) {
            SWIG_exception_fail(SWIG_ArgError(res), "Expected OEScalarGrid or None.");
        }
        $1 = reinterpret_cast< $1_ltype >(argp);
    }
}

%typemap(typecheck, precedence=10) const OESystem::OEScalarGrid* {
    if ($input == Py_None) {
        $1 = 1;
    } else {
        void *vptr = 0;
        int res = SWIG_ConvertPtr($input, &vptr, $descriptor, 0);
        $1 = SWIG_IsOK(res) ? 1 : _maptitude_is_oescalargrid($input) ? 1 : 0;
    }
}

// ============================================================================
// Typemap: OEUnaryPredicate<OEAtomBase>* (optional atom predicate mask)
// ============================================================================
%typemap(in) const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* (void *argp = 0) {
    if ($input == Py_None) {
        $1 = NULL;
    } else {
        // Try to extract the predicate pointer - OpenEye predicates are SWIG-wrapped
        int res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
        if (!SWIG_IsOK(res)) {
            // Fall back to direct pointer extraction for openeye predicates
            argp = _maptitude_extract_swig_ptr($input);
            if (!argp) {
                SWIG_exception_fail(SWIG_ArgError(SWIG_TypeError),
                    "Expected OEUnaryAtomPred or None for mask parameter.");
            }
        }
        $1 = reinterpret_cast< $1_ltype >(argp);
    }
}

%typemap(typecheck, precedence=10) const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* {
    $1 = ($input == Py_None) ? 1 : 1;  // Accept any object; runtime check in %typemap(in)
}

// ============================================================================
// Return type typemap: OEScalarGrid* → OpenEye-native Python OEScalarGrid
//
// Functions that return OEScalarGrid* allocate on the heap.  The custom
// %typemap(out) wraps each pointer as an openeye.oegrid.OEScalarGrid Python
// object (zero-copy pointer swap), so the OpenEye SWIG module's destructor
// handles cleanup.  This avoids the "memory leak, no destructor found"
// warning that occurs when maptitude's SWIG module owns a foreign type.
// ============================================================================
%typemap(out) OESystem::OEScalarGrid* {
    $result = _maptitude_wrap_as_oe_grid($1);
    if (!$result) SWIG_fail;
}

// ============================================================================
// Include STL typemaps
// ============================================================================
%include "std_string.i"
%include "std_vector.i"
%include "std_map.i"
%include "std_array.i"
%include "stdint.i"
%include "exception.i"

// ============================================================================
// Exception handling
// ============================================================================
%exception {
    try {
        $action
    } catch (const Maptitude::StructureError& e) {
        SWIG_exception(SWIG_ValueError, e.what());
    } catch (const Maptitude::GridError& e) {
        SWIG_exception(SWIG_ValueError, e.what());
    } catch (const Maptitude::SymOpError& e) {
        SWIG_exception(SWIG_ValueError, e.what());
    } catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (...) {
        SWIG_exception(SWIG_RuntimeError, "Unknown C++ exception");
    }
}

// ============================================================================
// Template instantiations for container types
// ============================================================================
%template(DoubleVector) std::vector<double>;
%template(UnsignedIntVector) std::vector<unsigned int>;
%template(SymOpVector) std::vector<Maptitude::SymOp>;
%template(ResidueDoubleMap) std::map<Maptitude::Residue, double>;
%template(UIntDoubleMap) std::map<unsigned int, double>;
%template(Double3Array) std::array<double, 3>;
%template(Double9Array) std::array<double, 9>;

// ============================================================================
// Version macros
// ============================================================================
#define MAPTITUDE_VERSION_MAJOR 0
#define MAPTITUDE_VERSION_MINOR 1
#define MAPTITUDE_VERSION_PATCH 0

// ============================================================================
// MapOp enum
// ============================================================================
namespace Maptitude {

enum class MapOp {
    Add,
    Subtract,
    Min,
    Max
};

// ============================================================================
// Residue struct
// ============================================================================
struct Residue {
    std::string name;
    int number;
    std::string chain;
    std::string insert_code;

    Residue();
    Residue(std::string name, int number, std::string chain,
            std::string insert_code = " ");

    static Residue FromAtom(const OEChem::OEAtomBase& atom);
    std::string ToString() const;

    bool operator==(const Residue& other) const;
    bool operator!=(const Residue& other) const;
    bool operator<(const Residue& other) const;
};

// ============================================================================
// UnitCell struct
// ============================================================================
struct UnitCell {
    double a;
    double b;
    double c;
    double alpha;
    double beta;
    double gamma;

    UnitCell();
    UnitCell(double a, double b, double c,
             double alpha, double beta, double gamma);

    double Volume() const;
    std::array<double, 9> OrthogonalizationMatrix() const;
    std::array<double, 9> DeorthogonalizationMatrix() const;
    std::array<double, 3> CartesianToFractional(double x, double y, double z) const;
    std::array<double, 3> FractionalToCartesian(double u, double v, double w) const;
    std::string ToString() const;

    bool operator==(const UnitCell& other) const;
    bool operator!=(const UnitCell& other) const;
};

// ============================================================================
// SymOp struct
// ============================================================================
struct SymOp {
    std::array<double, 9> R;
    std::array<double, 3> t;

    SymOp();
    SymOp(std::array<double, 9> rotation, std::array<double, 3> translation);

    static SymOp Parse(const std::string& triplet);
    static std::vector<SymOp> ParseAll(const std::string& text);

    std::array<double, 3> Apply(double u, double v, double w) const;
    std::string ToString() const;

    bool operator==(const SymOp& other) const;
    bool operator!=(const SymOp& other) const;
};

// ============================================================================
// RadialSampling enum + QScoreOptions class
// ============================================================================
enum class RadialSampling {
    Fixed,
    Adaptive
};

class QScoreOptions {
public:
    void SetSigma(double sigma);
    double GetSigma() const;
    void SetRadialStep(double d_rad);
    double GetRadialStep() const;
    void SetMaxRadius(double to_rad);
    double GetMaxRadius() const;
    void SetNumPoints(unsigned int num_points);
    unsigned int GetNumPoints() const;
    void SetNormalizeMap(bool normalize);
    bool GetNormalizeMap() const;
    void SetIsolatePoints(bool isolate);
    bool GetIsolatePoints() const;
    void SetRadialSampling(RadialSampling method);
    RadialSampling GetRadialSampling() const;
};

// ============================================================================
// AtomRadius enum + RSCCOptions class
// ============================================================================
enum class AtomRadius {
    Fixed,
    Scaled,
    Binned
};

class RSCCOptions {
public:
    void SetAtomRadiusMethod(AtomRadius method);
    AtomRadius GetAtomRadiusMethod() const;
    void SetFixedAtomRadius(double radius);
    double GetFixedAtomRadius() const;
    void SetAtomRadiusScaling(double scaling);
    double GetAtomRadiusScaling() const;
};

// ============================================================================
// DensityScoreResult struct
// ============================================================================
struct DensityScoreResult {
    double overall;
    std::map<Residue, double> by_residue;
    std::map<unsigned int, double> by_atom;
    std::string ToString() const;
};

// ============================================================================
// GridParams struct
// ============================================================================
struct GridParams {
    double x_origin;
    double y_origin;
    double z_origin;
    unsigned int x_dim;
    unsigned int y_dim;
    unsigned int z_dim;
    double spacing;
};

// ============================================================================
// Grid utility functions
// ============================================================================
GridParams GetGridParams(const OESystem::OEScalarGrid& grid);
double InterpolateDensity(const OESystem::OEScalarGrid& grid,
                          double x, double y, double z,
                          double default_value = 0.0);
std::vector<double> InterpolateDensityBatch(
    const OESystem::OEScalarGrid& grid,
    const std::vector<double>& points,
    size_t num_points,
    double default_value = 0.0);
std::vector<unsigned int> GetAtomGridPoints(
    const OESystem::OEScalarGrid& grid,
    double x, double y, double z, double radius);
double InterpolateDensityPeriodic(
    const OESystem::OEScalarGrid& grid,
    double x, double y, double z,
    double cell_a, double cell_b, double cell_c,
    double default_value = 0.0);
std::vector<double> InterpolateDensityPeriodicBatch(
    const OESystem::OEScalarGrid& grid,
    const std::vector<double>& points,
    size_t num_points,
    double cell_a, double cell_b, double cell_c,
    double default_value = 0.0);

// ============================================================================
// DensityCalculator class
// ============================================================================
class DensityCalculator {
public:
    DensityCalculator(const UnitCell& cell, const std::vector<SymOp>& symops);
    ~DensityCalculator();

    OESystem::OEScalarGrid* Calculate(
        OEChem::OEMolBase& mol,
        const OESystem::OEScalarGrid& obs_grid,
        double resolution,
        const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask = nullptr,
        double k_sol = 0.35,
        double b_sol = 46.0,
        bool include_h = false,
        unsigned int n_scale_shells = 1) const;
};

// ============================================================================
// Density scoring functions
// ============================================================================
DensityScoreResult RSCC(
    OEChem::OEMolBase& mol,
    const OESystem::OEScalarGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask = nullptr,
    const OESystem::OEScalarGrid* calc_grid = nullptr,
    const RSCCOptions& options = RSCCOptions());

DensityScoreResult RSR(
    OEChem::OEMolBase& mol,
    const OESystem::OEScalarGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask = nullptr,
    const OESystem::OEScalarGrid* calc_grid = nullptr);

DensityScoreResult QScore(
    OEChem::OEMolBase& mol,
    const OESystem::OEScalarGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask = nullptr,
    const QScoreOptions& options = QScoreOptions());

DensityScoreResult EDIAm(
    OEChem::OEMolBase& mol,
    const OESystem::OEScalarGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask = nullptr);

DensityScoreResult Coverage(
    OEChem::OEMolBase& mol,
    const OESystem::OEScalarGrid& grid,
    double sigma = 1.0,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask = nullptr);

// ============================================================================
// Grid operations
// ============================================================================
void ScaleMap(OESystem::OEScalarGrid& grid, double factor);
OESystem::OEScalarGrid* CombineMaps(
    const OESystem::OEScalarGrid& lhs,
    const OESystem::OEScalarGrid& rhs,
    MapOp op);
OESystem::OEScalarGrid* DiffToCalc(
    const OESystem::OEScalarGrid& obs_grid,
    const OESystem::OEScalarGrid& diff_grid);
OESystem::OEScalarGrid* WrapAndPadGrid(
    const OESystem::OEScalarGrid& grid,
    OEChem::OEMolBase& mol,
    double cell_a, double cell_b, double cell_c,
    double padding = 3.0);

}  // namespace Maptitude

// ============================================================================
// Python extensions for Residue
// ============================================================================
%extend Maptitude::Residue {
%pythoncode %{
def __repr__(self):
    return f"Residue('{self.ToString()}')"

def __str__(self):
    return self.ToString()

def __hash__(self):
    return hash((self.name, self.number, self.chain, self.insert_code))
%}
}

// ============================================================================
// Python extensions for UnitCell
// ============================================================================
%extend Maptitude::UnitCell {
%pythoncode %{
def __repr__(self):
    return self.ToString()
%}
}

// ============================================================================
// Python extensions for SymOp
// ============================================================================
%extend Maptitude::SymOp {
%pythoncode %{
def __repr__(self):
    return f"SymOp('{self.ToString()}')"

def __str__(self):
    return self.ToString()
%}
}

// ============================================================================
// Python extensions for DensityScoreResult
// ============================================================================
%extend Maptitude::DensityScoreResult {
%pythoncode %{
def __repr__(self):
    return self.ToString()
%}
}

// ============================================================================
// Python properties for QScoreOptions
// ============================================================================
%extend Maptitude::QScoreOptions {
%pythoncode %{
sigma = property(GetSigma, SetSigma)
radial_step = property(GetRadialStep, SetRadialStep)
max_radius = property(GetMaxRadius, SetMaxRadius)
num_points = property(GetNumPoints, SetNumPoints)
normalize_map = property(GetNormalizeMap, SetNormalizeMap)
isolate_points = property(GetIsolatePoints, SetIsolatePoints)
radial_sampling = property(GetRadialSampling, SetRadialSampling)
%}
}

// ============================================================================
// Python properties for RSCCOptions
// ============================================================================
%extend Maptitude::RSCCOptions {
%pythoncode %{
atom_radius_method = property(GetAtomRadiusMethod, SetAtomRadiusMethod)
fixed_atom_radius = property(GetFixedAtomRadius, SetFixedAtomRadius)
atom_radius_scaling = property(GetAtomRadiusScaling, SetAtomRadiusScaling)
%}
}

// ============================================================================
// Module-level Python convenience functions
// ============================================================================
%pythoncode %{
def fc_density(obj, obs_grid, resolution, cell, mask=None,
               k_sol=0.35, b_sol=46.0, symops=None,
               include_h=False, n_scale_shells=1):
    """Compute model electron density via Fourier synthesis.

    :param obj: Input molecule (OEMolBase or OEDesignUnit).
    :param obs_grid: Observed electron density grid (OEScalarGrid).
    :param resolution: Resolution limit in Angstroms.
    :param cell: UnitCell parameters.
    :param mask: Optional atom predicate to restrict contributing atoms.
    :param k_sol: Bulk solvent scale factor (default: 0.35 e/A^3).
    :param b_sol: Bulk solvent B-factor (default: 46.0 A^2).
    :param symops: List of SymOp or string of symmetry operators.
    :param include_h: Include hydrogen atoms (default: False).
    :param n_scale_shells: Number of per-shell scaling bins (default: 1).
    :returns: OEScalarGrid with computed model density.
    """
    if symops is None:
        symops = [SymOp()]  # Identity only
    elif isinstance(symops, str):
        symops = list(SymOp.ParseAll(symops))

    if isinstance(symops, list) and len(symops) > 0 and not isinstance(symops[0], SymOp):
        # Convert from SymOpVector if needed
        symops = list(symops)

    calc = DensityCalculator(cell, SymOpVector(symops))
    return calc.Calculate(obj, obs_grid, resolution, mask,
                          k_sol, b_sol, include_h, n_scale_shells)


def rscc(obj, grid, resolution, mask=None, calc_grid=None,
         atom_radius=None, options=None):
    """Real-Space Correlation Coefficient.

    :param obj: Input molecule.
    :param grid: Observed electron density map.
    :param resolution: Resolution in Angstroms.
    :param mask: Optional atom predicate.
    :param calc_grid: Optional pre-computed calculated density.
    :param atom_radius: Atom radius method (str or AtomRadius enum value).
    :param options: RSCCOptions configuration object.
    :returns: DensityScoreResult with RSCC values.
    """
    if options is None:
        options = RSCCOptions()
    if atom_radius is not None:
        if isinstance(atom_radius, str):
            _radius_map = {
                "fixed": AtomRadius_Fixed,
                "scaled": AtomRadius_Scaled,
                "binned": AtomRadius_Binned,
            }
            atom_radius = _radius_map[atom_radius.lower()]
        options.SetAtomRadiusMethod(atom_radius)
    return RSCC(obj, grid, resolution, mask, calc_grid, options)


def rsr(obj, grid, resolution, mask=None, calc_grid=None):
    """Real-Space R-Factor.

    :param obj: Input molecule.
    :param grid: Observed electron density map.
    :param resolution: Resolution in Angstroms.
    :param mask: Optional atom predicate.
    :param calc_grid: Optional pre-computed calculated density.
    :returns: DensityScoreResult with RSR values.
    """
    return RSR(obj, grid, resolution, mask, calc_grid)


def qscore(obj, grid, resolution, mask=None, options=None):
    """Q-Score (Pintilie et al., 2020).

    :param obj: Input molecule.
    :param grid: Observed electron density map.
    :param resolution: Resolution in Angstroms.
    :param mask: Optional atom predicate.
    :param options: QScoreOptions configuration.
    :returns: DensityScoreResult with Q-score values.
    """
    if options is None:
        options = QScoreOptions()
    return QScore(obj, grid, resolution, mask, options)


def ediam(obj, grid, resolution, mask=None):
    """Electron Density Index Averaged, Modified (EDIAm).

    :param obj: Input molecule.
    :param grid: Observed electron density map.
    :param resolution: Resolution in Angstroms.
    :param mask: Optional atom predicate.
    :returns: DensityScoreResult with EDIAm values in [0, 1].
    """
    return EDIAm(obj, grid, resolution, mask)


def coverage(obj, grid, sigma=1.0, mask=None):
    """Coverage: fraction of atoms observed in density.

    :param obj: Input molecule.
    :param grid: Observed electron density map.
    :param sigma: Number of standard deviations above mean.
    :param mask: Optional atom predicate.
    :returns: DensityScoreResult with coverage fractions.
    """
    return Coverage(obj, grid, sigma, mask)


def parse_symop(s):
    """Parse a single symmetry operator from triplet notation.

    :param s: Symmetry operator string (e.g., "x,y,z" or "-x,y+1/2,-z").
    :returns: SymOp object.
    """
    return SymOp.Parse(s)


def parse_symops(text):
    """Parse multiple symmetry operators from newline/semicolon-separated text.

    :param text: Block of symmetry operators.
    :returns: List of SymOp objects.
    """
    return list(SymOp.ParseAll(text))


def scale_map(grid, factor):
    """Scale a grid by multiplying all values by a scalar.

    :param grid: Grid to scale (modified in place).
    :param factor: Scale factor.
    """
    ScaleMap(grid, factor)


def combine_maps(lhs, rhs, op):
    """Combine two grids element-wise.

    :param lhs: Left-hand side grid.
    :param rhs: Right-hand side grid.
    :param op: MapOp enum value (Add, Subtract, Min, Max).
    :returns: New OEScalarGrid with combined values.
    """
    return CombineMaps(lhs, rhs, op)


def diff_to_calc(obs_grid, diff_grid):
    """Derive calculated density from observed and difference maps.

    :param obs_grid: Observed density map (2mFo-DFc).
    :param diff_grid: Difference density map (mFo-DFc).
    :returns: New OEScalarGrid with calculated density.
    """
    return DiffToCalc(obs_grid, diff_grid)


def wrap_and_pad_grid(grid, mol, cell_a, cell_b, cell_c, padding=3.0):
    """Translate molecule into unit cell and pad grid if needed.

    :param grid: CCP4 unit-cell grid.
    :param mol: Molecule to wrap (modified in-place).
    :param cell_a: Unit cell dimension a (Angstroms).
    :param cell_b: Unit cell dimension b (Angstroms).
    :param cell_c: Unit cell dimension c (Angstroms).
    :param padding: Extra margin around atoms (Angstroms).
    :returns: Grid covering all atom coordinates (original or padded).
    """
    result = WrapAndPadGrid(grid, mol, cell_a, cell_b, cell_c, padding)
    return result if result is not None else grid


__version__ = "0.1.0"
%}
