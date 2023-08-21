#include <stdlib.h>
#include <string.h>
#include "sparse_pattern.h"

static int compare_index(const void *a, const void *b) {
    IndexType x = *(IndexType *)a;
    IndexType y = *(IndexType *)b;
    if(x < y) {
        return -1;
    } else if(x > y) {
        return 1;
    } else {
        return 0;
    }
}

static IndexType sort_and_front_unique(IndexType *indices, IndexType length) {
    qsort(indices, length, sizeof(IndexType), compare_index);
    IndexType end = 0;
    for(IndexType i = 1; i < length; ++i) {
        if(indices[i] != indices[end]) {
            indices[++end] = indices[i];
        }
    }
    return end + 1;
}

void init_sparsity_pattern(SparsityPattern *pSparsityPattern,
                           IndexType nnode, 
                           IndexType nelem,
                           IndexType NSHLmax,
                           IndexType *ELMNSHL,
                           IndexType *ien) {

    pSparsityPattern->nnode = nnode;
    // generate sparsity pattern as a CSR graph from a mesh
    // the mesh is given by ien
    // ien is a flattened array of size nelem * NSHLmax
    // ELMNSHL is an array of size nelem
    // ELMNSHL[i] is the number of nodes in element i
    // ien[i + j * nelem] is the global node index of the j-th node in element i
    // the sparsity pattern is stored in CSR format
    // pSparsityPattern->offset is an array of size nnode + 1
    // pSparsityPattern->offset[i] is the offset of the i-th node in pSparsityPattern->indices
    // pSparsityPattern->indices is an array of size nnz
    // pSparsityPattern->indices[pSparsityPattern->offset[i] + j] is the j-th neighbor of the i-th node
    // pSparsityPattern->nnz is the number of nonzeros in the sparsity pattern
    // pSparsityPattern->offset[nnode] is the total number of nonzeros in the sparsity pattern

    // first, we need to know the number of appearances of each node
    IndexType *node_appearances = (IndexType *)malloc(sizeof(IndexType) * nnode);
    memset(node_appearances, 0, sizeof(IndexType) * nnode);
    for(IndexType i = 0; i < nelem; ++i) {
        for(IndexType j = 0; j < ELMNSHL[i]; ++j) {
            node_appearances[ien[i + j * nelem]] += ELMNSHL[i] - 1;
        }
    }
    // then, we calculate the prefix sum of node_appearances
    IndexType *node_offsets = (IndexType *)malloc(sizeof(IndexType) * (nnode + 1));
    node_offsets[0] = 0;
    for(IndexType i = 0; i < nnode; ++i) {
        node_offsets[i + 1] = node_offsets[i] + node_appearances[i];
    }
    // then, we collect the indices of the neighbors of each node
    IndexType *node_indices = (IndexType *)malloc(sizeof(IndexType) * node_offsets[nnode]);
    memset(node_appearances, 0, sizeof(IndexType) * nnode);
    for(IndexType i = 0; i < nelem; ++i) {
        for(IndexType j = 0; j < ELMNSHL[i]; ++j) {
            IndexType node = ien[i + j * nelem];
            node_indices[node_offsets[node] + node_appearances[node]] = i;
            node_appearances[node]++;
        }
    }
    // then, we sort the indices of the neighbors of each node and remove duplicates
    for(IndexType i = 0; i < nnode; ++i) {
        IndexType offset = node_offsets[i];
        IndexType length = node_appearances[i];
        IndexType *indices = node_indices + offset;
        node_appearances[i] = sort_and_front_unique(indices, length); 
    }
    // then, we calculate the total number of nonzeros
    IndexType nnz = 0;
    for(IndexType i = 0; i < nnode; ++i) {
        nnz += node_appearances[i];
    }
    // then, we allocate memory for the sparsity pattern
    pSparsityPattern->offset = (IndexType *)malloc(sizeof(IndexType) * (nnode + 1));
    pSparsityPattern->indices = (IndexType *)malloc(sizeof(IndexType) * nnz);
    // then, we calculate the sparsity pattern based on the new node_indices and node_apperances
    pSparsityPattern->offset[0] = 0;
    IndexType *offset = pSparsityPattern->offset;
    IndexType *indices = pSparsityPattern->indices;
    for(IndexType i = 0; i < nnode; ++i) {
        offset[i + 1] = offset[i] + node_appearances[i];
        memcpy(indices + offset[i], node_indices + node_offsets[i], sizeof(IndexType) * node_appearances[i]);
    }


    free(node_appearances);
    free(node_offsets);
    free(node_indices);

}

void free_sparsity_pattern(SparsityPattern *pSparsityPattern) {
    free(pSparsityPattern->offset);
    free(pSparsityPattern->indices);
}