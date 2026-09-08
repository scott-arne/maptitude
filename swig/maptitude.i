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
#include "maptitude/RsccOptions.h"
#include "maptitude/RsrOptions.h"
#include "maptitude/CoverageOptions.h"
#include "maptitude/DensityCalculator.h"
#include "maptitude/Metric.h"
#include "maptitude/GridOps.h"
#include "maptitude/MapIO.h"
#include "maptitude/SpatialIndex.h"

#include <oechem.h>
#include <oegrid.h>
#include <cstdio>
#include <new>

using namespace Maptitude;

/* Python exception objects, created on the first %init. Module-level statics:
   the interpreter owns the references for the life of the process.

   Sharing one set of classes across every interpreter makes this module unsafe
   under subinterpreters. That is a known, deliberate constraint, and these
   statics are only one reason for it -- the library also links OpenEye and
   OpenMP, neither of which is subinterpreter-safe either. */
static PyObject* g_maptitude_error = NULL;
static PyObject* g_structure_error = NULL;
static PyObject* g_grid_error = NULL;
static PyObject* g_symop_error = NULL;
static PyObject* g_cell_error = NULL;

/* Create one exception class, or reuse the cached one, then bind it on the
   module being executed. Returns 0 on success, or -1 with a Python error
   already set.

   The NULL check is load-bearing rather than defensive. SWIG emits %init into
   a Py_mod_exec slot and gives the generated PyModuleDef an m_size of 0, so
   the slot can run more than once per process -- dropping the maptitude keys
   from sys.modules and importing again is enough. Minting a fresh class on the
   second run would leave every earlier importer holding a superseded object,
   and the class %exception raises would no longer be the one their `except`
   clause names. Reusing the cached object is what keeps those identical, and
   what keeps the four subclasses sharing one MaptitudeError base.

   PyModule_AddObjectRef does not steal, so the module takes its own reference
   and the static keeps the one PyErr_NewException returned -- %exception
   dereferences that static long after %init has returned. */
static int MaptitudeAddException(PyObject* module, const char* qualified_name,
                                 const char* attribute_name, PyObject* base,
                                 PyObject** slot) {
    if (*slot == NULL) {
        /* The static keeps this reference for the life of the process; it is
           never released, because %exception dereferences it from arbitrary
           wrapper functions with no teardown hook to coordinate with. */
        *slot = PyErr_NewException(qualified_name, base, NULL);
        if (*slot == NULL) {
            return -1;
        }
    }
    return PyModule_AddObjectRef(module, attribute_name, *slot) < 0 ? -1 : 0;
}
%}

%init %{
    /* SWIG emits this block into SWIG_mod_exec, a Py_mod_exec slot returning
       int -- 0 for success, -1 for failure. Returning NULL here would be 0,
       i.e. a successful import with NULL exception statics. */
    if (MaptitudeAddException(m, "maptitude.MaptitudeError", "MaptitudeError",
                              NULL, &g_maptitude_error) < 0) {
        return -1;
    }
    if (MaptitudeAddException(m, "maptitude.StructureError", "StructureError",
                              g_maptitude_error, &g_structure_error) < 0) {
        return -1;
    }
    if (MaptitudeAddException(m, "maptitude.GridError", "GridError",
                              g_maptitude_error, &g_grid_error) < 0) {
        return -1;
    }
    if (MaptitudeAddException(m, "maptitude.SymOpError", "SymOpError",
                              g_maptitude_error, &g_symop_error) < 0) {
        return -1;
    }
    if (MaptitudeAddException(m, "maptitude.CellError", "CellError",
                              g_maptitude_error, &g_cell_error) < 0) {
        return -1;
    }
%}

// ============================================================================
// Forward declarations for cross-module SWIG type resolution
// ============================================================================
// These enable typemaps for OpenEye types whose definitions live in the
// OpenEye SWIG runtime (v4). Only types you actually use in your wrapped API
// need full #include — forward declarations suffice for the typemaps.

namespace OEChem {
    class OEMolBase;
    class OEMCMolBase;
    class OEMol;
    class OEGraphMol;
    class OEAtomBase;
    class OEBondBase;
    class OEConfBase;
    class OEMatchBase;
    class OEMolDatabase;
    class oemolistream;
    class oemolostream;
    class OEQMol;
    class OEResidue;
    class OEUniMolecularRxn;
}

namespace OEBio {
    class OEDesignUnit;
    class OEHierView;
    class OEHierResidue;
    class OEHierFragment;
    class OEHierChain;
    class OEInteractionHint;
    class OEInteractionHintContainer;
}

namespace OEDocking {
    class OEReceptor;
}

namespace OEPlatform {
    class oeifstream;
    class oeofstream;
    class oeisstream;
    class oeosstream;
}

namespace OESystem {
    class OESkewGrid;
    class OERecord;
    class OEMolRecord;
}

namespace OESystem {
    template <class T> class OEUnaryPredicate;
}

// ============================================================================
// Cross-runtime SWIG compatibility layer
// ============================================================================
// OpenEye's Python bindings use SWIG runtime v4; our module uses v5.
// Since the runtimes are separate, SWIG_TypeQuery cannot access OpenEye types.
// We validate the Python wrapper's real type (not isinstance, which honors
// __class__ properties) and directly extract the void* pointer from the
// SwigPyObject struct layout (stable across SWIG versions). Type validation
// is limited to the wrapper because the separate runtimes put OpenEye's SWIG
// type table out of reach.
//
// This approach enables passing OpenEye objects between Python and C++ without
// serialization. The macros below generate the boilerplate for each type.

%{
// Minimal SwigPyObject layout compatible across SWIG runtime versions.
// The actual struct may have more fields, but ptr is always first after
// PyObject_HEAD.
struct _SwigPyObjectCompat {
    PyObject_HEAD
    void *ptr;
};

/* Defense in depth, NOT a fix. Struct punning across two SWIG runtimes cannot
   be made safe, only unreachable from bad input: this function still casts an
   arbitrary object's `this` attribute to _SwigPyObjectCompat and reads through
   it, and no guard here can validate that layout. The type checkers in the
   typemaps are the primary defense: they validate the Python wrapper's real
   type before extraction. They cannot validate the pointer itself, because
   `this` is a writable attribute — an accepted predicate whose `this` was
   reassigned (e.g., pred.this = mol.this) will still be punned. Closing that
   would require validating the pointee against OpenEye's SWIG type table,
   which is unreachable across the v4/v5 runtime split. Do not relax a
   typemap's type check on the belief that this function is hardened -- it is not. */
static void* _maptitude_extract_swig_ptr(PyObject* obj) {
    if (obj == NULL || obj == Py_None) {
        return NULL;
    }
    PyObject* thisAttr = PyObject_GetAttrString(obj, "this");
    if (!thisAttr) {
        PyErr_Clear();
        return NULL;
    }
    void* ptr = ((_SwigPyObjectCompat*)thisAttr)->ptr;
    Py_DECREF(thisAttr);
    return ptr;
}

// ---- Type checker generator macro ----
// Generates a cached type checker for an OpenEye Python type. Uses the object's
// real type (Py_TYPE), not PyObject_IsInstance, because the extracted pointer is
// reinterpret_cast to a C++ type and a spoofed __class__ property would otherwise
// let any object claim to be one. This check covers the wrapper type only; see
// the _maptitude_extract_swig_ptr comment for the pointer validation limit.
// TAG:    identifier suffix (e.g., oemolbase)
// MODULE: Python module string (e.g., "openeye.oechem")
// CLASS:  Python class name string (e.g., "OEMolBase")
#define DEFINE_OE_TYPE_CHECKER(TAG, MODULE, CLASS) \
    static PyObject* _maptitude_oe_##TAG##_type = NULL; \
    static bool _maptitude_is_##TAG(PyObject* obj) { \
        if (!_maptitude_oe_##TAG##_type) { \
            PyObject* mod = PyImport_ImportModule(MODULE); \
            if (mod) { \
                PyObject* cls = PyObject_GetAttrString(mod, CLASS); \
                Py_DECREF(mod); \
                /* Only a real type object supports a spoof-resistant check. */ \
                if (cls && !PyType_Check(cls)) { \
                    Py_DECREF(cls); \
                    cls = NULL; \
                } \
                _maptitude_oe_##TAG##_type = cls; \
            } \
            if (!_maptitude_oe_##TAG##_type) { \
                PyErr_Clear(); \
                return false; \
            } \
        } \
        if (obj == NULL) return false; \
        return PyType_IsSubtype(Py_TYPE(obj), \
                                (PyTypeObject*)_maptitude_oe_##TAG##_type) != 0; \
    }

// ---- Molecule types (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oemolbase,    "openeye.oechem", "OEMolBase")
DEFINE_OE_TYPE_CHECKER(oemcmolbase,  "openeye.oechem", "OEMCMolBase")
DEFINE_OE_TYPE_CHECKER(oemol,        "openeye.oechem", "OEMol")
DEFINE_OE_TYPE_CHECKER(oegraphmol,   "openeye.oechem", "OEGraphMol")
DEFINE_OE_TYPE_CHECKER(oeqmol,       "openeye.oechem", "OEQMol")

// ---- Atom / bond / conformer / residue (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oeatombase,   "openeye.oechem", "OEAtomBase")
DEFINE_OE_TYPE_CHECKER(oebondbase,   "openeye.oechem", "OEBondBase")
DEFINE_OE_TYPE_CHECKER(oeconfbase,   "openeye.oechem", "OEConfBase")
DEFINE_OE_TYPE_CHECKER(oeresidue,    "openeye.oechem", "OEResidue")
DEFINE_OE_TYPE_CHECKER(oematchbase,  "openeye.oechem", "OEMatchBase")

// ---- Molecule I/O (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oemolistream, "openeye.oechem", "oemolistream")
DEFINE_OE_TYPE_CHECKER(oemolostream, "openeye.oechem", "oemolostream")
DEFINE_OE_TYPE_CHECKER(oemoldatabase,"openeye.oechem", "OEMolDatabase")

// ---- Reactions (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oeunimolecularrxn, "openeye.oechem", "OEUniMolecularRxn")

// ---- Platform streams (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oeifstream,   "openeye.oechem", "oeifstream")
DEFINE_OE_TYPE_CHECKER(oeofstream,   "openeye.oechem", "oeofstream")
DEFINE_OE_TYPE_CHECKER(oeisstream,   "openeye.oechem", "oeisstream")
DEFINE_OE_TYPE_CHECKER(oeosstream,   "openeye.oechem", "oeosstream")

// ---- Records (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oerecord,     "openeye.oechem", "OERecord")
DEFINE_OE_TYPE_CHECKER(oemolrecord,  "openeye.oechem", "OEMolRecord")

// ---- Bio / hierarchy (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oedesignunit, "openeye.oechem", "OEDesignUnit")
DEFINE_OE_TYPE_CHECKER(oehierview,   "openeye.oechem", "OEHierView")
DEFINE_OE_TYPE_CHECKER(oehierresidue,"openeye.oechem", "OEHierResidue")
DEFINE_OE_TYPE_CHECKER(oehierfragment,"openeye.oechem","OEHierFragment")
DEFINE_OE_TYPE_CHECKER(oehierchain,  "openeye.oechem", "OEHierChain")
DEFINE_OE_TYPE_CHECKER(oeinteractionhint,          "openeye.oechem", "OEInteractionHint")
DEFINE_OE_TYPE_CHECKER(oeinteractionhintcontainer, "openeye.oechem", "OEInteractionHintContainer")

// ---- Predicates (openeye.oechem) ----
DEFINE_OE_TYPE_CHECKER(oeunaryatompred, "openeye.oechem", "OEUnaryAtomPred")

// ---- Grid (openeye.oegrid) ----
DEFINE_OE_TYPE_CHECKER(oeskewgrid, "openeye.oegrid", "OESkewGrid")

// ---- Docking (openeye.oedocking) ----
DEFINE_OE_TYPE_CHECKER(oereceptor,   "openeye.oedocking", "OEReceptor")

#undef DEFINE_OE_TYPE_CHECKER

/* ---- OESkewGrid return-type helper (copy-assign into Python object) ---- */
static PyObject* _maptitude_wrap_as_oe_skew_grid(OESystem::OESkewGrid* grid) {
    if (!grid) {
        Py_RETURN_NONE;
    }
    PyObject* oegrid_mod = PyImport_ImportModule("openeye.oegrid");
    if (!oegrid_mod) {
        delete grid;
        return NULL;
    }
    PyObject* grid_cls = PyObject_GetAttrString(oegrid_mod, "OESkewGrid");
    Py_DECREF(oegrid_mod);
    if (!grid_cls) {
        delete grid;
        return NULL;
    }
    PyObject* oe_grid = PyObject_CallNoArgs(grid_cls);
    Py_DECREF(grid_cls);
    if (!oe_grid) {
        delete grid;
        return NULL;
    }
    /* Validate the constructed object's type before extracting its pointer.
       Restores the documented symmetry with the in-direction typemaps, which
       all check the type. This catches a wrong-typed object constructed from a
       rebound OESkewGrid name, but not a poisoned type cache (if the very
       first use happens after the rebind) or a correct-typed object whose
       'this' attribute has been reassigned. */
    if (!_maptitude_is_oeskewgrid(oe_grid)) {
        Py_DECREF(oe_grid);
        delete grid;
        PyErr_SetString(PyExc_TypeError,
                        "constructed object is not an OESkewGrid");
        return NULL;
    }
    PyObject* thisAttr = PyObject_GetAttrString(oe_grid, "this");
    if (!thisAttr) {
        PyErr_Clear();
        Py_DECREF(oe_grid);
        delete grid;
        PyErr_SetString(PyExc_RuntimeError,
                        "the constructed OESkewGrid has no 'this' attribute");
        return NULL;
    }

    /* Copy the value into the Python-side grid rather than swapping pointers.
       The grid under `thisAttr` was allocated inside OpenEye's shared library;
       deleting it here would run this module's operator delete on another
       runtime's allocation, and pointing it at our grid would then have
       OpenEye's destructor free our memory. Each allocation is released by the
       allocator that made it.

       Assignment, not an element loop: `oe_grid` is default-constructed above,
       so it starts 1x1x1; an element loop would copy exactly one value.
       operator= resizes the destination and copies dimensions, the unit cell,
       the space group, midpoints, title and data together.

       The out typemap is emitted outside SWIG's generated try/catch, so an
       exception here would reach CPython unhandled and abort the interpreter.
       operator= can throw std::bad_alloc during reallocation. */
    _SwigPyObjectCompat* swig_this = (_SwigPyObjectCompat*)thisAttr;
    OESystem::OESkewGrid* dest =
        reinterpret_cast<OESystem::OESkewGrid*>(swig_this->ptr);
    if (dest == NULL) {
        Py_DECREF(thisAttr);
        Py_DECREF(oe_grid);
        delete grid;
        PyErr_SetString(PyExc_RuntimeError,
                        "failed to access the wrapped OESkewGrid");
        return NULL;
    }
    try {
        *dest = *grid;
    } catch (const std::bad_alloc&) {
        Py_DECREF(thisAttr);
        Py_DECREF(oe_grid);
        delete grid;
        return PyErr_NoMemory();
    } catch (const std::exception& e) {
        Py_DECREF(thisAttr);
        Py_DECREF(oe_grid);
        delete grid;
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    } catch (...) {
        Py_DECREF(thisAttr);
        Py_DECREF(oe_grid);
        delete grid;
        PyErr_SetString(PyExc_RuntimeError,
                        "unknown C++ exception while copying the grid");
        return NULL;
    }

    /* Validate the copy. operator= routes through OpenEye's geometry
       setters, which reject values they cannot represent by returning false
       rather than throwing, and operator= ignores that return. The destination
       then silently keeps its default 1x1x1 geometry.

       same_grid_geometry cannot answer "different" for that 1x1x1 destination
       -- it throws GridError, because a one-node axis has no derivable node
       interval. So the comparison runs under a catch: on the failure this
       block exists to detect, the throw is the detection, and any other
       exception it raises is equally a reason not to hand the grid back. */
    bool copy_ok = false;
    bool cell_failure = false;
    std::string detail;
    try {
        copy_ok = Maptitude::same_grid_geometry(*dest, *grid) &&
                  dest->HasUnitCell() == grid->HasUnitCell() &&
                  dest->HasSpaceGroup() == grid->HasSpaceGroup() &&
                  dest->GetSize() == grid->GetSize();
        if (copy_ok && dest->HasSpaceGroup()) {
            copy_ok = dest->GetSpaceGroup() == grid->GetSpaceGroup();
        }
    } catch (const Maptitude::CellError& e) {
        /* The source grid's own sampling is unrepresentable -- read_map can be
           handed such a file. Nothing about the copy failed, so this is
           reported as the CellError the rest of the package documents for
           non-axis-aligned sampling rather than as a copy failure. */
        copy_ok = false;
        cell_failure = true;
        detail = e.what();
    } catch (const std::exception& e) {
        copy_ok = false;
        detail = e.what();
    } catch (...) {
        copy_ok = false;
        detail = "unknown C++ exception";
    }
    if (!copy_ok) {
        /* This runs after the %exception-protected call has returned, so the
           class is chosen here or the failure escapes as a bare RuntimeError
           that `except maptitude.MaptitudeError` cannot catch. The statics are
           set in %init and cannot be NULL at call time; the fallback keeps a
           NULL out of PyErr_SetString rather than trading one bug for a
           crash. */
        char errmsg[512];
        if (cell_failure) {
            /* The dimensions match on this path -- reporting them would point
               at a copy that did not fail. */
            std::snprintf(errmsg, sizeof(errmsg), "%s", detail.c_str());
        } else {
            std::snprintf(errmsg, sizeof(errmsg),
                          "Grid geometry copy failed: source is %ux%ux%u size=%u, "
                          "destination is %ux%ux%u size=%u%s%s",
                          grid->GetXDim(), grid->GetYDim(), grid->GetZDim(), grid->GetSize(),
                          dest->GetXDim(), dest->GetYDim(), dest->GetZDim(), dest->GetSize(),
                          detail.empty() ? "" : "; ", detail.c_str());
        }
        PyObject* error_class = cell_failure ? g_cell_error : g_grid_error;
        if (!error_class) {
            error_class = PyExc_RuntimeError;
        }
        Py_DECREF(thisAttr);
        Py_DECREF(oe_grid);
        delete grid;
        PyErr_SetString(error_class, errmsg);
        return NULL;
    }

    /* We own `grid`; the Python object owns `dest` and always has. */
    delete grid;
    Py_DECREF(thisAttr);
    return oe_grid;
}
%}

// ============================================================================
// Typemap generator macros
// ============================================================================

// Generate const-ref and non-const-ref typemaps for a cross-runtime OpenEye type.
%define OE_CROSS_RUNTIME_REF_TYPEMAPS(CPP_TYPE, CHECKER, ERR_MSG)

%typemap(in) const CPP_TYPE& (void *argp = 0, int res = 0) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
    if (!SWIG_IsOK(res)) {
        if (CHECKER($input)) {
            argp = _maptitude_extract_swig_ptr($input);
            if (argp) res = SWIG_OK;
        }
    }
    if (!SWIG_IsOK(res)) {
        SWIG_exception_fail(SWIG_ArgError(res), ERR_MSG);
    }
    if (!argp) {
        SWIG_exception_fail(SWIG_NullReferenceError, "Null reference.");
    }
    $1 = reinterpret_cast< $1_ltype >(argp);
}

%typemap(typecheck, precedence=10) const CPP_TYPE& {
    void *vptr = 0;
    int res = SWIG_ConvertPtr($input, &vptr, $descriptor, SWIG_POINTER_NO_NULL);
    $1 = SWIG_IsOK(res) ? 1 : CHECKER($input) ? 1 : 0;
}

%typemap(in) CPP_TYPE& (void *argp = 0, int res = 0) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
    if (!SWIG_IsOK(res)) {
        if (CHECKER($input)) {
            argp = _maptitude_extract_swig_ptr($input);
            if (argp) res = SWIG_OK;
        }
    }
    if (!SWIG_IsOK(res)) {
        SWIG_exception_fail(SWIG_ArgError(res), ERR_MSG);
    }
    if (!argp) {
        SWIG_exception_fail(SWIG_NullReferenceError, "Null reference.");
    }
    $1 = reinterpret_cast< $1_ltype >(argp);
}

%typemap(typecheck, precedence=10) CPP_TYPE& {
    void *vptr = 0;
    int res = SWIG_ConvertPtr($input, &vptr, $descriptor, SWIG_POINTER_NO_NULL);
    $1 = SWIG_IsOK(res) ? 1 : CHECKER($input) ? 1 : 0;
}

%enddef

// Generate nullable-pointer typemaps (accepts None) for a cross-runtime type.
%define OE_CROSS_RUNTIME_NULLABLE_PTR_TYPEMAPS(CPP_TYPE, CHECKER, ERR_MSG)

%typemap(in) const CPP_TYPE* (void *argp = 0, int res = 0) {
    if ($input == Py_None) {
        $1 = NULL;
    } else {
        res = SWIG_ConvertPtr($input, &argp, $descriptor, 0);
        if (!SWIG_IsOK(res)) {
            if (CHECKER($input)) {
                argp = _maptitude_extract_swig_ptr($input);
                if (argp) res = SWIG_OK;
            }
        }
        if (!SWIG_IsOK(res)) {
            SWIG_exception_fail(SWIG_ArgError(res), ERR_MSG);
        }
        $1 = reinterpret_cast< $1_ltype >(argp);
    }
}

%typemap(typecheck, precedence=10) const CPP_TYPE* {
    if ($input == Py_None) {
        $1 = 1;
    } else {
        void *vptr = 0;
        int res = SWIG_ConvertPtr($input, &vptr, $descriptor, 0);
        $1 = SWIG_IsOK(res) ? 1 : CHECKER($input) ? 1 : 0;
    }
}

%enddef

// ============================================================================
// Typemap declarations for all OpenEye types
// ============================================================================
// Each type gets const-ref and non-const-ref typemaps. Types that commonly
// appear as optional parameters also get nullable-pointer typemaps.
// These are inert until a wrapped function signature uses the type.

// ---- Molecule hierarchy (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMolBase,    _maptitude_is_oemolbase,    "Expected OEMolBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMCMolBase,  _maptitude_is_oemcmolbase,  "Expected OEMCMolBase-derived object (OEMCMolBase or OEMol).")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMol,        _maptitude_is_oemol,        "Expected OEMol object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEGraphMol,   _maptitude_is_oegraphmol,   "Expected OEGraphMol object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEQMol,       _maptitude_is_oeqmol,       "Expected OEQMol object.")

// ---- Atom / bond / conformer / residue / match (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEAtomBase,   _maptitude_is_oeatombase,   "Expected OEAtomBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEBondBase,   _maptitude_is_oebondbase,   "Expected OEBondBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEConfBase,   _maptitude_is_oeconfbase,   "Expected OEConfBase-derived object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEResidue,    _maptitude_is_oeresidue,    "Expected OEResidue object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMatchBase,  _maptitude_is_oematchbase,  "Expected OEMatchBase-derived object.")

// ---- Molecule I/O (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::oemolistream,  _maptitude_is_oemolistream, "Expected oemolistream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::oemolostream,  _maptitude_is_oemolostream, "Expected oemolostream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEMolDatabase, _maptitude_is_oemoldatabase,"Expected OEMolDatabase object.")

// ---- Reactions (OEChem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEChem::OEUniMolecularRxn, _maptitude_is_oeunimolecularrxn, "Expected OEUniMolecularRxn object.")

// ---- Platform streams (OEPlatform) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeifstream, _maptitude_is_oeifstream, "Expected oeifstream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeofstream, _maptitude_is_oeofstream, "Expected oeofstream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeisstream, _maptitude_is_oeisstream, "Expected oeisstream object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEPlatform::oeosstream, _maptitude_is_oeosstream, "Expected oeosstream object.")

// ---- Records (OESystem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OESystem::OERecord,    _maptitude_is_oerecord,    "Expected OERecord object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OESystem::OEMolRecord, _maptitude_is_oemolrecord, "Expected OEMolRecord object.")

// ---- Bio / hierarchy (OEBio) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEDesignUnit,   _maptitude_is_oedesignunit, "Expected OEDesignUnit object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierView,     _maptitude_is_oehierview,   "Expected OEHierView object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierResidue,  _maptitude_is_oehierresidue,"Expected OEHierResidue object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierFragment,  _maptitude_is_oehierfragment,"Expected OEHierFragment object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEHierChain,    _maptitude_is_oehierchain,  "Expected OEHierChain object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEInteractionHint,          _maptitude_is_oeinteractionhint,          "Expected OEInteractionHint object.")
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEBio::OEInteractionHintContainer, _maptitude_is_oeinteractionhintcontainer, "Expected OEInteractionHintContainer object.")

// ---- Grid (OESystem) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OESystem::OESkewGrid, _maptitude_is_oeskewgrid, "Expected OESkewGrid-derived object.")
OE_CROSS_RUNTIME_NULLABLE_PTR_TYPEMAPS(OESystem::OESkewGrid, _maptitude_is_oeskewgrid, "Expected OESkewGrid or None.")

// OESkewGrid return-type typemap (wraps C++ grid as native openeye.oegrid object)
%typemap(out) OESystem::OESkewGrid* {
    $result = _maptitude_wrap_as_oe_skew_grid($1);
    if (!$result) SWIG_fail;
}

// ---- Docking (OEDocking) ----
OE_CROSS_RUNTIME_REF_TYPEMAPS(OEDocking::OEReceptor, _maptitude_is_oereceptor, "Expected OEReceptor object.")

// ============================================================================
// Typemap: OEUnaryPredicate<OEAtomBase>* (maptitude-specific, optional atom predicate mask)
// ============================================================================
OE_CROSS_RUNTIME_NULLABLE_PTR_TYPEMAPS(OESystem::OEUnaryPredicate<OEChem::OEAtomBase>, _maptitude_is_oeunaryatompred, "Expected OEUnaryAtomPred or None for mask parameter.")

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
        PyErr_SetString(g_structure_error, e.what());
        SWIG_fail;
    } catch (const Maptitude::GridError& e) {
        PyErr_SetString(g_grid_error, e.what());
        SWIG_fail;
    } catch (const Maptitude::SymOpError& e) {
        PyErr_SetString(g_symop_error, e.what());
        SWIG_fail;
    } catch (const Maptitude::CellError& e) {
        PyErr_SetString(g_cell_error, e.what());
        SWIG_fail;
    } catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    } catch (...) {
        SWIG_exception(SWIG_RuntimeError, "Unknown C++ exception");
    }
}

%pythoncode %{
# The exception classes are created in the extension module's init function.
# SWIG's proxy module does not mirror arbitrary extension attributes, so alias
# them here; python/maptitude/__init__.py re-exports from this module.
MaptitudeError = _maptitude.MaptitudeError
StructureError = _maptitude.StructureError
GridError = _maptitude.GridError
SymOpError = _maptitude.SymOpError
CellError = _maptitude.CellError
%}

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
#define MAPTITUDE_VERSION_MINOR 5
#define MAPTITUDE_VERSION_PATCH 0

// ============================================================================
// MapOp enum
// ============================================================================
namespace Maptitude {

enum class MapOp {
    ADD,
    SUBTRACT,
    MIN,
    MAX
};

/// Which header record wins when both encode a nonzero origin.
enum class OriginSource {
    ORIGIN_RECORD,
    NXSTART
};

// SWIG's default scoped-enum conversion accepts any integer and casts it, so an
// out-of-range op reached the C++ switch and fell through to an all-zero grid
// with no error. These narrow the boundary to the declared enumerators.
%typemap(in) Maptitude::MapOp {
    if (PyBool_Check($input) || !PyLong_Check($input)) {
        SWIG_exception_fail(SWIG_TypeError,
                            "op must be one of maptitude.MapOp.ADD, "
                            ".SUBTRACT, .MIN or .MAX");
    }
    const long op_value = PyLong_AsLong($input);
    if (op_value == -1 && PyErr_Occurred()) SWIG_fail;
    if (op_value < 0 || op_value > 3) {
        PyErr_Format(PyExc_ValueError,
                     "MapOp value %ld is out of range; expected 0-3 "
                     "(ADD, SUBTRACT, MIN, MAX)", op_value);
        SWIG_fail;
    }
    $1 = static_cast<Maptitude::MapOp>(op_value);
}

%typemap(in) Maptitude::OriginSource {
    if (PyBool_Check($input) || !PyLong_Check($input)) {
        SWIG_exception_fail(SWIG_TypeError,
                            "tiebreak must be maptitude.OriginSource."
                            "ORIGIN_RECORD or .NXSTART");
    }
    const long source_value = PyLong_AsLong($input);
    if (source_value == -1 && PyErr_Occurred()) SWIG_fail;
    if (source_value < 0 || source_value > 1) {
        PyErr_Format(PyExc_ValueError,
                     "OriginSource value %ld is out of range; expected 0-1 "
                     "(ORIGIN_RECORD, NXSTART)", source_value);
        SWIG_fail;
    }
    $1 = static_cast<Maptitude::OriginSource>(source_value);
}

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
    FIXED,
    ADAPTIVE
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
// AtomRadius enum + RsccOptions / RsrOptions / CoverageOptions classes
// ============================================================================
enum class AtomRadius {
    FIXED,
    SCALED,
    BINNED,
    ADAPTIVE
};

class RsccOptions {
public:
    void SetAtomRadiusMethod(AtomRadius method);
    AtomRadius GetAtomRadiusMethod() const;
    void SetFixedAtomRadius(double radius);
    double GetFixedAtomRadius() const;
    void SetAtomRadiusScaling(double scaling);
    double GetAtomRadiusScaling() const;
};

class RsrOptions {
public:
    void SetAtomRadiusMethod(AtomRadius method);
    AtomRadius GetAtomRadiusMethod() const;
    void SetFixedAtomRadius(double radius);
    double GetFixedAtomRadius() const;
    void SetAtomRadiusScaling(double scaling);
    double GetAtomRadiusScaling() const;
};

class CoverageOptions {
public:
    void SetSigma(double sigma);
    double GetSigma() const;
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
    double x_spacing;
    double y_spacing;
    double z_spacing;
};

struct UnitCellParams {
    double a;
    double b;
    double c;
    double alpha;
    double beta;
    double gamma;
};

// ============================================================================
// Grid utility functions
// ============================================================================
GridParams get_grid_params(const OESystem::OESkewGrid& grid);
UnitCellParams get_unit_cell(const OESystem::OESkewGrid& grid);
bool grid_contains(const GridParams& gp, double x, double y, double z);
bool same_grid_geometry(const OESystem::OESkewGrid& lhs,
                        const OESystem::OESkewGrid& rhs,
                        double tol = 1e-6);
double interpolate_density(const OESystem::OESkewGrid& grid,
                          double x, double y, double z,
                          double default_value = 0.0);
std::vector<double> interpolate_density_batch(
    const OESystem::OESkewGrid& grid,
    const std::vector<double>& points,
    size_t num_points,
    double default_value = 0.0);
std::vector<unsigned int> get_atom_grid_points(
    const OESystem::OESkewGrid& grid,
    double x, double y, double z, double radius);
// Documented through %feature rather than a Python wrapper: both are hot paths
// called per point, and a wrapper would put an interpreter frame on every call
// to buy nothing but the docstring. The commensurability contract has to appear
// here because these two reach require_commensurate_cell directly, so a Python
// caller who never touches wrap_and_pad_grid would otherwise meet the CellError
// with no documented cause.
%feature("docstring") interpolate_density_periodic %{
Trilinearly interpolate density, wrapping the query point into the unit cell.

:param grid: Grid sampling one unit cell.
:param x: Query x (Angstroms).
:param y: Query y (Angstroms).
:param z: Query z (Angstroms).
:param cell_a: Unit cell dimension a (Angstroms).
:param cell_b: Unit cell dimension b (Angstroms).
:param cell_c: Unit cell dimension c (Angstroms).
:param default_value: Returned where the grid cannot supply a value.
:returns: Interpolated density at the wrapped point.
:raises GridError: If an axis has fewer than two nodes, a node has no spatial
    coordinate or a non-finite one, or a derived interval is not finite and
    positive.
:raises CellError: If the sampling is not axis-aligned, or if a cell edge does
    not round to ``n_i`` or ``n_i - 1`` of that axis's node intervals to within
    the allowance made for float node coordinates. Wrapping at a period the grid
    was never sampled on would resample the map onto a lattice it never had, so
    the caller is told rather than handed a plausible wrong number.
%}
double interpolate_density_periodic(
    const OESystem::OESkewGrid& grid,
    double x, double y, double z,
    double cell_a, double cell_b, double cell_c,
    double default_value = 0.0);
%feature("docstring") interpolate_density_periodic_batch %{
Interpolate many points at once, wrapping each into the unit cell.

:param grid: Grid sampling one unit cell.
:param points: Flat sequence of x, y, z triples.
:param num_points: Number of triples in ``points``.
:param cell_a: Unit cell dimension a (Angstroms).
:param cell_b: Unit cell dimension b (Angstroms).
:param cell_c: Unit cell dimension c (Angstroms).
:param default_value: Returned where the grid cannot supply a value.
:returns: One interpolated density per point, in input order.
:raises GridError: As :func:`interpolate_density_periodic`.
:raises CellError: As :func:`interpolate_density_periodic`. The cell is checked
    once for the whole batch, before any point is sampled.
%}
std::vector<double> interpolate_density_periodic_batch(
    const OESystem::OESkewGrid& grid,
    const std::vector<double>& points,
    size_t num_points,
    double cell_a, double cell_b, double cell_c,
    double default_value = 0.0);

// ============================================================================
// Scattering factor types and functions
// ============================================================================
struct CromerMannCoeffs {
    std::array<double, 4> a;
    std::array<double, 4> b;
    double c;
    double Evaluate(double s_squared) const;
};

struct ScatteringFactorEntry {
    uint8_t atomic_number;
    int8_t formal_charge;
    CromerMannCoeffs coeffs;
};

const CromerMannCoeffs* get_scattering_factors(
    unsigned int atomic_number, int formal_charge = 0);

// ============================================================================
// DensityCalculator class
// ============================================================================
class DensityCalculator {
public:
    DensityCalculator(const UnitCell& cell, const std::vector<SymOp>& symops);
    ~DensityCalculator();

    OESystem::OESkewGrid* Calculate(
        OEChem::OEMolBase& mol,
        const OESystem::OESkewGrid& obs_grid,
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
DensityScoreResult rscc(
    OEChem::OEMolBase& mol,
    const OESystem::OESkewGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask = nullptr,
    const OESystem::OESkewGrid* calc_grid = nullptr,
    const RsccOptions& options = RsccOptions());

DensityScoreResult rsr(
    OEChem::OEMolBase& mol,
    const OESystem::OESkewGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask = nullptr,
    const OESystem::OESkewGrid* calc_grid = nullptr,
    const RsrOptions& options = RsrOptions());

DensityScoreResult qscore(
    OEChem::OEMolBase& mol,
    const OESystem::OESkewGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask = nullptr,
    const QScoreOptions& options = QScoreOptions());

DensityScoreResult ediam(
    OEChem::OEMolBase& mol,
    const OESystem::OESkewGrid& grid,
    double resolution,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask = nullptr);

DensityScoreResult coverage(
    OEChem::OEMolBase& mol,
    const OESystem::OESkewGrid& grid,
    const OESystem::OEUnaryPredicate<OEChem::OEAtomBase>* mask = nullptr,
    const CoverageOptions& options = CoverageOptions());

// ============================================================================
// Grid operations
// ============================================================================
%feature("novaluewrapper") Maptitude::MapFile;

%typemap(out) Maptitude::MapFile {
    /* $1.grid is a unique_ptr; release it into the helper, which takes
       ownership and deletes the C++ grid after copy-assigning it into a
       Python-owned OESkewGrid. */
    PyObject* grid_obj = _maptitude_wrap_as_oe_skew_grid($1.grid.release());
    if (!grid_obj) SWIG_fail;
    PyObject* symops_obj = PyUnicode_FromStringAndSize(
        $1.symops.data(), (Py_ssize_t)$1.symops.size());
    if (!symops_obj) {
        Py_DECREF(grid_obj);
        SWIG_fail;
    }
    $result = PyTuple_New(2);
    if (!$result) {
        Py_DECREF(grid_obj);
        Py_DECREF(symops_obj);
        SWIG_fail;
    }
    PyTuple_SET_ITEM($result, 0, grid_obj);
    PyTuple_SET_ITEM($result, 1, symops_obj);
}

void scale_map(OESystem::OESkewGrid& grid, double factor);
OESystem::OESkewGrid* combine_maps(
    const OESystem::OESkewGrid& lhs,
    const OESystem::OESkewGrid& rhs,
    MapOp op);
OESystem::OESkewGrid* diff_to_calc(
    const OESystem::OESkewGrid& obs_grid,
    const OESystem::OESkewGrid& diff_grid);
OESystem::OESkewGrid* wrap_and_pad_grid(
    const OESystem::OESkewGrid& grid,
    OEChem::OEMolBase& mol,
    double cell_a, double cell_b, double cell_c,
    double padding = 3.0);

struct MapFile;

MapFile read_map(const std::string& path,
                 OriginSource tiebreak);

void write_map(const std::string& path,
               const OESystem::OESkewGrid& grid,
               const std::string& symops = "");

}  // namespace Maptitude

// ============================================================================
// Inline helper: get_scattering_factor_table as a vector
// ============================================================================
%inline %{
namespace Maptitude {
std::vector<Maptitude::ScatteringFactorEntry> _get_scattering_factor_table_vec() {
    size_t count = 0;
    const ScatteringFactorEntry* table = get_scattering_factor_table(count);
    return std::vector<ScatteringFactorEntry>(table, table + count);
}
}
%}

%template(ScatteringFactorEntryVector) std::vector<Maptitude::ScatteringFactorEntry>;

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
// Python properties for RsccOptions
// ============================================================================
%extend Maptitude::RsccOptions {
%pythoncode %{
atom_radius_method = property(GetAtomRadiusMethod, SetAtomRadiusMethod)
fixed_atom_radius = property(GetFixedAtomRadius, SetFixedAtomRadius)
atom_radius_scaling = property(GetAtomRadiusScaling, SetAtomRadiusScaling)
%}
}

// ============================================================================
// Python properties for RsrOptions
// ============================================================================
%extend Maptitude::RsrOptions {
%pythoncode %{
atom_radius_method = property(GetAtomRadiusMethod, SetAtomRadiusMethod)
fixed_atom_radius = property(GetFixedAtomRadius, SetFixedAtomRadius)
atom_radius_scaling = property(GetAtomRadiusScaling, SetAtomRadiusScaling)
%}
}

// ============================================================================
// Python properties for CoverageOptions
// ============================================================================
%extend Maptitude::CoverageOptions {
%pythoncode %{
sigma = property(GetSigma, SetSigma)
%}
}

// ============================================================================
// Module-level Python convenience functions
// ============================================================================
%pythoncode %{
import os

# Save references to SWIG-generated C++ wrappers before overriding
_cpp_rscc = rscc
_cpp_rsr = rsr
_cpp_qscore = qscore
_cpp_ediam = ediam
_cpp_coverage = coverage
_cpp_scale_map = scale_map
_cpp_combine_maps = combine_maps
_cpp_diff_to_calc = diff_to_calc
_cpp_wrap_and_pad_grid = wrap_and_pad_grid
_cpp_read_map = read_map
_cpp_write_map = write_map


def _lookup_atom_radius(metric, radius_map, name):
    """Resolve an atom-radius method name to its AtomRadius value.

    A bare ``radius_map[name.lower()]`` raised ``KeyError('vdw')`` on a
    misspelling: outside the typed exception hierarchy the README presents as the
    complete set to handle, and carrying no list of what would have worked.

    :param metric: Name of the calling metric, for the message.
    :param radius_map: Lower-case method name to AtomRadius value.
    :param name: The string the caller supplied.
    :returns: The AtomRadius value.
    :raises ValueError: If ``name`` matches no entry in ``radius_map``.
    """
    try:
        return radius_map[name.lower()]
    except KeyError:
        raise ValueError(
            "%s does not support atom_radius=%r; accepted values are %s"
            % (metric, name, ", ".join(repr(k) for k in sorted(radius_map)))
        ) from None


def _copy_rscc_options(options):
    """Return a copy of an RsccOptions so callers' objects are never mutated.

    :param options: Source options.
    :returns: An independent copy carrying the same settings.
    :raises TypeError: If ``options`` is not an RsccOptions.
    """
    if not isinstance(options, RsccOptions):
        raise TypeError(
            "options must be an RsccOptions, not %s" % type(options).__name__
        )
    copied = RsccOptions()
    copied.SetAtomRadiusMethod(options.GetAtomRadiusMethod())
    copied.SetFixedAtomRadius(options.GetFixedAtomRadius())
    copied.SetAtomRadiusScaling(options.GetAtomRadiusScaling())
    return copied


def _copy_rsr_options(options):
    """Return a copy of an RsrOptions so callers' objects are never mutated.

    :param options: Source options.
    :returns: An independent copy carrying the same settings.
    :raises TypeError: If ``options`` is not an RsrOptions.
    """
    if not isinstance(options, RsrOptions):
        raise TypeError(
            "options must be an RsrOptions, not %s" % type(options).__name__
        )
    copied = RsrOptions()
    copied.SetAtomRadiusMethod(options.GetAtomRadiusMethod())
    copied.SetFixedAtomRadius(options.GetFixedAtomRadius())
    copied.SetAtomRadiusScaling(options.GetAtomRadiusScaling())
    return copied


def _copy_coverage_options(options):
    """Return a copy of a CoverageOptions so callers' objects are never mutated.

    :param options: Source options.
    :returns: An independent copy carrying the same settings.
    :raises TypeError: If ``options`` is not a CoverageOptions.
    """
    if not isinstance(options, CoverageOptions):
        raise TypeError(
            "options must be a CoverageOptions, not %s" % type(options).__name__
        )
    copied = CoverageOptions()
    copied.SetSigma(options.GetSigma())
    return copied


def fc_density(obj, obs_grid, resolution, cell, mask=None,
               k_sol=0.35, b_sol=46.0, symops=None,
               include_h=False, n_scale_shells=1):
    """Compute model electron density via Fourier synthesis.

    :param obj: Input molecule (OEMolBase or OEDesignUnit).
        Atoms without a radius are assigned Bondi radii in place, as the
        scoring functions do.
    :param obs_grid: Observed electron density grid (OESkewGrid).
    :param resolution: Resolution limit in Angstroms.
    :param cell: UnitCell parameters.
    :param mask: Optional atom predicate to restrict contributing atoms.
    :param k_sol: Bulk solvent scale factor (default: 0.35 e/A^3).
    :param b_sol: Bulk solvent B-factor (default: 46.0 A^2).
    :param symops: Symmetry operators, as a semicolon-separated string, an
        iterable of SymOp, an iterable of operator strings, or None for the
        identity alone.
    :param include_h: Include hydrogen atoms (default: False).
    :param n_scale_shells: Number of per-shell scaling bins, in ``[1, 1000]``
        (default: 1). The upper bound is ``MAX_SCALE_SHELLS``: the bins partition
        the resolution range and even a 0.5 A dataset has far fewer independent
        shells, while a value near ``UINT_MAX`` sizes the shell-edge table into
        tens of gigabytes and at ``UINT_MAX`` itself wraps it to zero.
    :returns: OESkewGrid with computed model density.
    :raises TypeError: If ``symops`` is neither a string nor an iterable of
        SymOp or operator strings.
    :raises SymOpError: If an operator string cannot be parsed.
    :raises GridError: If ``resolution`` is not finite and positive, if
        ``n_scale_shells`` is outside ``[1, 1000]``, if ``obs_grid``'s spacing is
        at or above twice a cell edge, or if the resolution and the cell together
        need a Miller-index box of more than 2e8 points.
    :raises CellError: If ``cell`` is invalid or is not orthorhombic.
    """
    if symops is None:
        symops = [SymOp()]  # Identity only
    elif isinstance(symops, str):
        symops = list(SymOp.ParseAll(symops))
    elif isinstance(symops, SymOp):
        # A lone SymOp is not iterable, so it would reach the generic message below
        # and be reported as "not SymOp" -- a type this function does accept, inside
        # a sequence. Name the wrapping the caller has to do instead.
        raise TypeError(
            "symops must be a sequence of SymOp, not a single SymOp; pass [symops]"
        )
    elif isinstance(symops, (bytes, bytearray)):
        # bytes and bytearray iterate as ints, so the per-element check below would
        # report "not int" for input the caller never wrote as integers.
        raise TypeError(
            "symops must be a str, not %s; decode it first" % type(symops).__name__
        )
    else:
        try:
            elements = list(symops)
        except TypeError:
            raise TypeError(
                "symops must be a string, an iterable of SymOp, or an iterable of "
                "operator strings, not %s" % type(symops).__name__
            ) from None
        # A list of operator strings is the shape this docstring describes, so
        # accept it here rather than letting SymOpVector reject it with a raw SWIG
        # overload dump. ParseAll is the same path the bare-string branch takes.
        normalized = []
        for element in elements:
            if isinstance(element, SymOp):
                normalized.append(element)
            elif isinstance(element, str):
                normalized.extend(SymOp.ParseAll(element))
            else:
                raise TypeError(
                    "symops entries must be SymOp or str, not %s"
                    % type(element).__name__
                )
        symops = normalized

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
    :param atom_radius: Atom radius method, as an AtomRadius enum value or one of
        the strings ``"fixed"``, ``"scaled"``, ``"binned"`` (case-insensitive).
        ``"adaptive"`` is not accepted: rscc has no adaptive radius model.
    :param options: RsccOptions configuration object. Never mutated.
    :returns: DensityScoreResult with RSCC values.
    :raises ValueError: If ``atom_radius`` is a string naming no supported method,
        including ``"adaptive"``.
    :raises RuntimeError: If ``atom_radius`` is ``AtomRadius.ADAPTIVE``, or if
        ``options`` carries it. The two spellings of the same rejection raise
        different types because they are found in different places: the string is
        resolved against a table here, where "no such method" is a ValueError,
        while the enum value is a real method that this metric does not implement
        and is refused by the C++ layer, whose ``std::invalid_argument`` surfaces
        as RuntimeError like every other option-value rejection.
    """
    if options is None:
        options = RsccOptions()
    if atom_radius is not None:
        if isinstance(atom_radius, str):
            # Three entries, not four: rscc has no adaptive radius model and the
            # enum path rejects AtomRadius.ADAPTIVE, so omitting it here is the
            # agreeing behavior rather than an oversight.
            _radius_map = {
                "fixed": AtomRadius_FIXED,
                "scaled": AtomRadius_SCALED,
                "binned": AtomRadius_BINNED,
            }
            atom_radius = _lookup_atom_radius("rscc", _radius_map, atom_radius)
        # Copy: mutating the caller's options object would leak this call's
        # settings into their next call.
        options = _copy_rscc_options(options)
        options.SetAtomRadiusMethod(atom_radius)
    return _cpp_rscc(obj, grid, resolution, mask, calc_grid, options)


def rsr(obj, grid, resolution, mask=None, calc_grid=None,
        atom_radius=None, options=None):
    """Real-Space R-Factor.

    :param obj: Input molecule.
    :param grid: Observed electron density map.
    :param resolution: Resolution in Angstroms.
    :param mask: Optional atom predicate.
    :param calc_grid: Optional pre-computed calculated density.
    :param atom_radius: Atom radius method, as an AtomRadius enum value or one of
        the strings ``"fixed"``, ``"scaled"``, ``"binned"``, ``"adaptive"``
        (case-insensitive).
    :param options: RsrOptions configuration object. Never mutated.
    :returns: DensityScoreResult with RSR values.
    :raises ValueError: If ``atom_radius`` is a string naming no supported method.
    """
    if options is None:
        options = RsrOptions()
    if atom_radius is not None:
        if isinstance(atom_radius, str):
            _radius_map = {
                "fixed": AtomRadius_FIXED,
                "scaled": AtomRadius_SCALED,
                "binned": AtomRadius_BINNED,
                "adaptive": AtomRadius_ADAPTIVE,
            }
            atom_radius = _lookup_atom_radius("rsr", _radius_map, atom_radius)
        # Copy: mutating the caller's options object would leak this call's
        # settings into their next call.
        options = _copy_rsr_options(options)
        options.SetAtomRadiusMethod(atom_radius)
    return _cpp_rsr(obj, grid, resolution, mask, calc_grid, options)


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
    return _cpp_qscore(obj, grid, resolution, mask, options)


def ediam(obj, grid, resolution, mask=None):
    """Electron Density Index Averaged, Modified (EDIAm).

    :param obj: Input molecule.
    :param grid: Observed electron density map.
    :param resolution: Resolution in Angstroms.
    :param mask: Optional atom predicate.
    :returns: DensityScoreResult with EDIAm values in [0, 1].
    """
    return _cpp_ediam(obj, grid, resolution, mask)


def coverage(obj, grid, sigma=None, mask=None, options=None):
    """Coverage: fraction of atoms observed in density.

    :param obj: Input molecule.
    :param grid: Observed electron density map.
    :param sigma: Number of standard deviations above mean. When None, the
        value carried by ``options`` is used. An explicit value overrides it.
    :param mask: Optional atom predicate.
    :param options: CoverageOptions configuration object. Never mutated.
    :returns: DensityScoreResult with coverage fractions.
    """
    if options is None:
        options = CoverageOptions()
    if sigma is not None:
        options = _copy_coverage_options(options)
        options.SetSigma(sigma)
    return _cpp_coverage(obj, grid, mask, options)


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
    _cpp_scale_map(grid, factor)


def combine_maps(lhs, rhs, op):
    """Combine two grids element-wise.

    :param lhs: Left-hand side grid.
    :param rhs: Right-hand side grid.
    :param op: MapOp enum value (ADD, SUBTRACT, MIN, MAX).
    :returns: New OESkewGrid with combined values.
    :raises TypeError: When ``op`` is not a Python ``int``, or is a ``bool``.
        ``bool`` subclasses ``int``, but neither ``True`` nor ``False`` names an
        operation, so both are rejected. A NumPy integer scalar is rejected here
        too, since the check is on the Python type.
    :raises ValueError: When ``op`` is an ``int`` outside 0-3 that fits in a C
        ``long``. One too large to fit raises ``OverflowError`` instead.
    """
    return _cpp_combine_maps(lhs, rhs, op)


def diff_to_calc(obs_grid, diff_grid):
    """Derive calculated density from observed and difference maps.

    :param obs_grid: Observed density map (2mFo-DFc).
    :param diff_grid: Difference density map (mFo-DFc).
    :returns: New OESkewGrid with calculated density.
    """
    return _cpp_diff_to_calc(obs_grid, diff_grid)


class MapFile(tuple):
    """The contents of a map file: its density grid and its symmetry text.

    A 2-tuple of ``(grid, symops)`` with the two members also reachable by
    name, so ``result.grid`` and ``grid, symops = read_map(path)`` both work.
    """

    __slots__ = ()

    def __new__(cls, grid, symops):
        return tuple.__new__(cls, (grid, symops))

    @property
    def grid(self):
        """The density carrier, an ``OESkewGrid`` at the file's origin."""
        return self[0]

    @property
    def symops(self):
        """Newline-separated triplets; ``""`` when the file carries none."""
        return self[1]

    def __repr__(self):
        return "MapFile(grid={0!r}, symops={1!r})".format(self[0], self[1])


def read_map(path, tiebreak=OriginSource_ORIGIN_RECORD):
    """Read a CCP4 or MRC map, preserving its ORIGIN record and symmetry block.

    OEReadGrid returns the payload, the cell and the space group but never
    consults the MRC2000 ORIGIN record and does not expose the symmetry block.
    This reads both and returns a grid already moved onto the file's origin.

    :param path: Map file to read. A ``str`` or any :class:`os.PathLike`.
    :param tiebreak: Which record wins when ORIGIN and NxSTART both encode a
        nonzero origin, whether or not the two agree -- nothing here compares
        them. One of ``OriginSource.ORIGIN_RECORD`` (the default) or
        ``OriginSource.NXSTART``. Ignored when at most one is nonzero.
    :returns: A :class:`MapFile` of ``(grid, symops)``.
    :raises TypeError: If ``path`` is neither a ``str`` nor an
        :class:`os.PathLike`. A ``bytes`` path is also rejected, by the
        ``std::string`` typemap rather than by ``os.fspath``. Also when
        ``tiebreak`` is not a Python ``int``, or is a ``bool``. ``bool``
        subclasses ``int``, but neither ``True`` nor ``False`` names a tiebreak
        rule, so both are rejected. A NumPy integer scalar is rejected here
        too, since the check is on the Python type.
    :raises CellError: If the file's sampling is not axis-aligned. maptitude
        requires an orthorhombic cell; a skewed one is rejected on read.
    :raises GridError: If ``path``'s extension is not one OpenEye maps to the
        CCP4 format, or names a compressed file such as ``.ccp4.gz``: this
        function reads the file's own first 1024 bytes as a CCP4 header, and
        neither admits that. The admitted class is the one :func:`write_map`'s
        ``path`` documents. Also if the file cannot be read, its header cannot
        be parsed, or the copy that returns the grid to Python did not preserve
        the source geometry -- a degenerate cell edge reaches this. Also if
        ``path`` contains an embedded NUL, which ``os.fspath`` passes through
        unchanged: the calls that open the file stop at the NUL, so the name
        given and the name opened are not the same name.
    :raises SymOpError: If the symmetry block is present and does not parse,
        or if one of its 80-byte records still holds a ``;`` or a newline once
        its padding is stripped. A CCP4 record carries one operator, so a
        separator inside one separates nothing on disk; a record holding two
        triplets that way is refused rather than read as two operators, which
        would give a :func:`write_map` round trip one more record than the file
        has.
    :raises ValueError: When ``tiebreak`` is an ``int`` outside 0-1 that fits
        in a C ``long``. One too large to fit raises ``OverflowError``.

    Example::

        result = read_map("2fofc.ccp4")
        symops = parse_symops(result.symops) if result.symops else None
    """
    grid, symops = _cpp_read_map(os.fspath(path), tiebreak)
    return MapFile(grid, symops)


def write_map(path, grid, symops=""):
    """Write a grid as CCP4 or MRC, restoring the ORIGIN and symmetry records.

    The file is written and verified under a temporary name and renamed onto
    ``path`` only once it reads back as the grid it came from, so a grid that
    cannot be written faithfully raises with ``path`` untouched rather than
    leaving a wrong map on disk or destroying a good one.

    A relative ``path`` is resolved against the working directory once, on
    entry, and that resolved path is what the temporary is placed beside and
    what the rename targets; a thread moving the working directory during the
    write changes nothing after that point.

    An absent space group is defaulted to P1 on the way out, so a successful
    write of a grid carrying none lands at ``ISPG`` 1 rather than 0. Measured
    on ``tests/assets/mapq/390_emd_30342_A_z4.mrc``, whose ``ISPG`` is 0: the
    file this writes for it carries 1. The default is what keeps the payload
    intact -- an unset space group makes ``OEWriteGrid`` double the cell and
    regrid, and it does so silently, returning success over a self-consistent
    header describing a different grid. Defaulting happens on a copy, so
    ``grid`` itself is untouched; it is that copy the re-read map is compared
    against, not ``grid``.

    A successful write replaces the destination's inode rather than rewriting
    it in place, so three properties of an existing destination do not survive.
    Its mode resets to whatever an ordinary create gives, which under POSIX is
    0666 narrowed by the process umask, widening a permission the caller had
    narrowed. Hard links to it keep the old contents under their own names. A
    destination that was a symlink becomes a regular file, with its former
    target left untouched.

    The verification covers the bytes this function wrote, not the bytes that
    arrive at ``path``. It writes and checks a temporary beside the destination
    and then renames that onto ``path``, and every step addresses the temporary
    by path. So a process able to write the destination's directory can replace
    the temporary between the last check and the rename, and this function will
    publish its bytes and return successfully. Such a process can already
    create, replace and remove ``path`` itself; what this adds is that a
    successful return stops implying the published bytes are the ones that were
    verified. A destination directory only the caller can write closes the
    replacement, while the process umask denies others write on the temporary
    itself: under POSIX it is created at 0666 narrowed by the umask, so a umask
    that leaves it group- or world-writable lets a process able to enter the
    directory write its bytes without being able to replace it.

    Node 0 is written twice, into the MRC2000 ORIGIN record exactly and into
    NxSTART as an integer node count. A map this function wrote reproduces node
    0 to float32 under ``OriginSource.ORIGIN_RECORD`` and only to half a node
    interval per axis under ``OriginSource.NXSTART``, and the write is refused
    unless the two records agree on every axis to within that bound.

    That bound is per axis, and there is no second one over the three together:
    each axis is compared against its own node interval. So the straight-line
    distance between the two placements can reach the root-sum-square of the
    three half-intervals -- ``sqrt(3)`` times half a node interval on a grid
    sampled equally on all three axes -- and a caller budgeting one distance
    rather than three per-axis bounds needs that larger figure.

    The padded grid :func:`wrap_and_pad_grid` allocates when padding is needed
    is refused: its declared cell forces the written node count one higher per
    axis than the grid carries. That is a property of the newly allocated grid,
    not of every value that function returns: when no padding is needed it
    hands back the grid it was given.

    :param path: Destination. The format follows the extension; the admitted
        class is every spelling whose first three characters are ``ccp``,
        ``map`` or ``mrc``, compared case-insensitively. ``.ccp4``, ``.mrc``
        and ``.map`` are the spellings this writer is tested on, not the
        accepted set. A compressed destination such as ``.ccp4.gz`` is refused,
        as is a filename that is nothing but an extension. A ``str`` or any
        :class:`os.PathLike`.
    :param grid: Grid to write, an ``OESkewGrid``.
    :param symops: Symmetry text in the form :func:`read_map` returns, one
        triplet per line. Empty writes no symmetry block. The semicolon
        separator is normalized to a newline, so the operator set survives
        a round trip but that spelling does not.
    :raises TypeError: If ``path`` is neither a ``str`` nor an :class:`os.PathLike`.
        A ``bytes`` path is also rejected, by the ``std::string`` typemap rather than
        by ``os.fspath``.
    :raises CellError: If the grid's sampling is not axis-aligned. maptitude
        requires an orthorhombic cell; a skewed one is rejected on write.
    :raises GridError: If the extension is not one OpenEye maps to the CCP4
        format or names a compressed file, if the temporary this function
        writes beside ``path`` cannot be created, which an unwritable or absent
        destination directory produces, if eight attempts at a temporary name
        beside ``path`` all collide with an existing file, if the write fails,
        if the re-read map differs from the space-group-defaulted copy
        described above -- not from ``grid`` itself -- in dimensions, cell,
        per-axis spacing, node 0 or any voxel, if either the grid or the
        re-read map has an axis with fewer than two nodes, so its spacing is
        undefined, if on any axis the ORIGIN and NxSTART records place node 0
        more than half that axis's node interval apart, or if the verified
        temporary file cannot be renamed onto ``path``. Also if ``path``
        contains an embedded NUL, which ``os.fspath`` passes through
        unchanged: the extension gate reads the whole string while the calls
        that write and publish the file stop at the NUL, so the name checked
        and the name written are not the same name. That list is not closed:
        this function raises ``GridError`` from further internal checks, among
        them the header re-reads and the NSYMBT and symop-block byte
        comparisons. Catch the class rather than switching on the list.
    :raises SymOpError: If ``symops`` is non-empty and does not parse, or if
        any of its records is longer than the format's 80-character field.
    """
    _cpp_write_map(os.fspath(path), grid, symops)


def wrap_and_pad_grid(grid, mol, cell_a, cell_b, cell_c, padding=3.0):
    """Translate molecule into unit cell and pad grid if needed.

    :param grid: CCP4 unit-cell grid.
    :param mol: Molecule to wrap (modified in-place).
    :param cell_a: Unit cell dimension a (Angstroms). Must be finite and positive.
    :param cell_b: Unit cell dimension b (Angstroms). Must be finite and positive.
    :param cell_c: Unit cell dimension c (Angstroms). Must be finite and positive.
    :param padding: Extra margin around atoms (Angstroms).
    :returns: A padded grid, or the original grid when no padding is needed.
        Never ``None``: the C++ function returns a null grid to mean "no padding
        was needed", and this wrapper substitutes the original.
    :raises StructureError: If the molecule contains no heavy atoms.
    :raises CellError: If any cell edge is zero, negative, or non-finite. The
        wrap uses ``fmod`` against each edge, which is NaN for a zero divisor and
        previously filled the padded grid with NaN voxels. Also if an edge does
        not round to ``n_i`` or ``n_i - 1`` of that axis's node intervals, to
        within the allowance made for float node coordinates: the padded grid is
        filled by sampling this one periodically, which needs the grid to tile
        the cell it is given.
    """
    result = _cpp_wrap_and_pad_grid(grid, mol, cell_a, cell_b, cell_c, padding)
    return result if result is not None else grid


def get_scattering_factor_table():
    """Get the full scattering factor table.

    :returns: Tuple of (table_entries, count).
    """
    entries = _get_scattering_factor_table_vec()
    return entries, len(entries)


__version__ = "0.5.0"
%}
