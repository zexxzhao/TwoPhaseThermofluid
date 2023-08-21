#pragma once

#include "../sys/defines.h"

struct SparsityPattern {
    IndexType nnode;
    IndexType nnz;
    IndexType *indices;
    IndexType *offsets;
};

void init_sparsity_pattern(SparsityPattern *pSparsityPattern,
                           IndexType nnode, 
                           IndexType nelem,
                           IndexType NSHLmax,
                           IndexType *ELMNSHL,
                           IndexType *ien);

void free_sparsity_pattern(SparsityPattern *pSparsityPattern);