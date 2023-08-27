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
    FEM_RT_ERROR = 1,
    FEM_MEM_ERROR = 2,
    FEM_FILE_ERROR = 3,
    FEM_INVALID_INPUT = 4
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

typedef enum {
    FEM_LOG_LEVEL_NONE = 100,
    FEM_LOG_LEVEL_ERROR = 200,
    FEM_LOG_LEVEL_WARNING = 300,
    FEM_LOG_LEVEL_INFO = 400,
    FEM_LOG_LEVEL_DEBUG = 500
} FEMLogLevelType;

#define FEM_LOG_LEVEL_DEFAULT FEM_LOG_LEVEL_INFO


#define FEM_LOG(LOG_LEVEL, ...) do {\
    if (LOG_LEVEL <= FEM_LOG_LEVEL_DEFAULT) {\
        fprintf(stderr, __VA_ARGS__);\
    }\
} while (0)

