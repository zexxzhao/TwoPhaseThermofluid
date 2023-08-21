#pragma once

#include <stdint.h>
typedef uint32_t FEMIndexType;
typedef double FEMScalarType;

typedef FEMIndexType IndexType;
typedef FEMScalarType ScalarType;

#include <stdbool.h>
typedef bool FEMBoolType;


/* BUGGY for double-evaluation */
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))

/* Mangling for Fortran */

#define FORTRAN_NAME(x) x##_


/* FLAG on which field to assemble*/
typedef enum {
    ASSEMBLE_FIELD_NONE = 0,
    ASSEMBLE_FIELD_NS = 1,
    ASSEMBLE_FIELD_LS = 2,
    ASSEMBLE_FIELD_VOF = ASSEMBLE_FIELD_LS,
    ASSEMBLE_FIELD_TEM = 4
} AssembleFieldType;

typedef enum {
    ASSEMBLE_TENSOR_NONE = 0,
    ASSEMBLE_TENSOR_SCALAR = 1,
    ASSEMBLE_TENSOR_VEC = 2,
    ASSEMBLE_TENSOR_MAT = 4
} AssembleTensorType;

typedef enum {
    FEM_SUCCESS = 0,
    FEM_RUNTIME_ERROR = 1,
    FEM_OUT_OF_BOUNDS = 2,
    FEM_INVALID_INPUT = 4,
    FEM_INVALID_STATE = 8,
    FEM_INVALID_OPERATION = 16
} FEMErrorCode;


#define FEM_CHECK_ERROR(run) do {\
    FEMErrorCode err = (run);\
    if (err != FEM_SUCCESS) {\
        fprintf(stderr, "Error %d at %s:%d\n", err, __FILE__, __LINE__);\
        exit(err);\
    }\
} while (0)

#ifdef PETSC_VERSION_MAJOR
#define HAVE_PETSC
#endif

#ifdef MPI_VERSION
#define HAVE_MPI
#endif