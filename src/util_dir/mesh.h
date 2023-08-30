#pragma once

#include <sys_dir/defines.h>
#include <stdint.h>
#include <stddef.h>

#include <stdint.h>

typedef struct {
    int32_t FACE_ID;
    int32_t NFACE;

    int32_t* FACE_IEN; // This was a 2D array in Fortran
    int32_t* F2E;
    int32_t* FACE_OR;
    int32_t* NSHLB;
    int32_t* NGAUSSB;

    int32_t NNODE;
    int32_t* BNODES;

    int32_t* L2SNODE;
    int32_t* L2SELEM;

} bnd_class;

typedef struct {
    int32_t NSD, NNODE, NELEM, NBOUND;
    int32_t NPATCH, NSHLBmax, maxNSHL;
    double* xg; // 2D array, assuming real(8) is double in Fortran
    
    int32_t* IEN;  // This was a 2D array in Fortran
    int32_t* NodeID;
    int32_t* ELM_ID;

    bnd_class* bound;
    int32_t* ELMNSHL;
    int32_t* ELMNGAUSS;
} MeshData;


// Methods
void read_mesh(MeshData* mesh);
void free_mesh(MeshData* mesh);
