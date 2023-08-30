#include "util_dir/mesh.h"
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>

static void input(int id, MeshData* mesh) {
    FILE* mfid;
    char fname[64];
    char iname[9];
    int i, j, k;
    int NSD, NNODE, NSHLmax, NELEM, NBOUND, NPATCH;

    sprintf(fname, "part.%d.dat", id);
    mfid = fopen(fname, "r");

    fscanf(mfid, "%d %d %d", &mesh->NSD, &mesh->maxNSHL, &mesh->NSHLBmax);
    fscanf(mfid, "%d %d %d %d", &mesh->NNODE, &mesh->NELEM, &mesh->NBOUND, &NPATCH);
    
    NSD = mesh->NSD;
    NSHLmax = mesh->maxNSHL;
    NNODE = mesh->NNODE;
    NELEM = mesh->NELEM;
    NBOUND = mesh->NBOUND;

    // read nodes
    mesh->xg = malloc(NNODE * NSD * sizeof(double));
    mesh->NodeID = malloc(NNODE * sizeof(int32_t));
    for (i = 0; i < NNODE; i++) {
        for (j = 0; j < NSD; j++) {
            fscanf(mfid, "%lf", &mesh->xg[i * NSD + j]);
        }
        fscanf(mfid, "%d", &mesh->NodeID[i]);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // read elements
    mesh->ELMNSHL = malloc(NELEM * sizeof(int32_t));
    mesh->IEN = malloc(NELEM * NSHLmax * sizeof(int32_t));
    mesh->ELM_ID = malloc(NELEM * sizeof(int32_t));
    for (i = 0; i < NELEM; i++) {
        fscanf(mfid, "%d", &mesh->ELMNSHL[i]);
        for (j = 0; j < mesh->ELMNSHL[i]; j++) {
            fscanf(mfid, "%d", &mesh->IEN[i * NSHLmax + j]);
            mesh->IEN[i * NSHLmax + j] -= 1;
        }
        fscanf(mfid, "%d", &mesh->ELM_ID[i]);
    }


    // read faces
    mesh->bound = malloc(NBOUND * sizeof(bnd_class));
    for (i = 0; i < NBOUND; i++) {
        fscanf(mfid, "%d %d %d", &mesh->bound[i].FACE_ID, &mesh->bound[i].NFACE, &mesh->bound[i].NNODE);

        mesh->bound[i].FACE_IEN = malloc(mesh->bound[i].NFACE * mesh->NSHLBmax * sizeof(int32_t));
        mesh->bound[i].F2E = malloc(mesh->bound[i].NFACE * sizeof(int32_t));
        mesh->bound[i].FACE_OR = malloc(mesh->bound[i].NFACE * sizeof(int32_t));
        mesh->bound[i].NSHLB = malloc(mesh->bound[i].NFACE * sizeof(int32_t));
        for (j = 0; j < mesh->bound[i].NFACE; j++) {
            mesh->bound[i].NSHLB[j] = mesh->NSHLBmax;
            for (k = 0; k < mesh->NSHLBmax; k++) {
                fscanf(mfid, "%d", &mesh->bound[i].FACE_IEN[j * mesh->NSHLBmax + k]);
            }
        }
        for (j = 0; j < mesh->bound[i].NFACE; j++) {
            fscanf(mfid, "%d %d", &mesh->bound[i].F2E[j], &mesh->bound[i].FACE_OR[j]);
        }
        mesh->bound[i].BNODES = malloc(mesh->bound[i].NNODE * sizeof(int32_t));
        for (j = 0; j < mesh->bound[i].NNODE; j++) {
            fscanf(mfid, "%d", &mesh->bound[i].BNODES[j]);
        }
    }
    fclose(mfid);

    mesh->ELMNGAUSS = malloc(NELEM * sizeof(int32_t));
    for (i = 0; i < NELEM; i++) {
        mesh->ELMNGAUSS[i] = mesh->ELMNSHL[i];
    }
    for (i = 0; i < NBOUND; i++) {
        mesh->bound[i].NGAUSSB = malloc(mesh->bound[i].NFACE * sizeof(int32_t));
        for (j = 0; j < mesh->bound[i].NFACE; j++) {
            mesh->bound[i].NGAUSSB[j] = mesh->bound[i].NSHLB[j];
        }
    }

}


static void shuffle_node_id() {

}

void read_mesh(MeshData* mesh) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    input(rank+1, mesh);
}

void free_mesh(MeshData* mesh) {
    int i;
    free(mesh->xg);
    free(mesh->NodeID);
    free(mesh->ELMNSHL);
    free(mesh->IEN);
    free(mesh->ELM_ID);
    for (i = 0; i < mesh->NBOUND; i++) {
        free(mesh->bound[i].FACE_IEN);
        free(mesh->bound[i].F2E);
        free(mesh->bound[i].FACE_OR);
        free(mesh->bound[i].NSHLB);
        free(mesh->bound[i].BNODES);
        free(mesh->bound[i].NGAUSSB);
    }
    free(mesh->bound);
    free(mesh->ELMNGAUSS);
}