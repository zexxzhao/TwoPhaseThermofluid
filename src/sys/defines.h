#pragma once

typedef int IndexType;
typedef double ScalarType;



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
