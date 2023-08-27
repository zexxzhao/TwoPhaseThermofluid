#pragma once

#include <sys_dir/defines.h>

struct SparsityPattern {
    IndexType nnode;
    IndexType nnz;
    IndexType *indices;
    IndexType *offset;
};

typedef struct SparsityPattern SparsityPattern;

void init_sparsity_pattern(SparsityPattern *pSparsityPattern,
                           IndexType nnode, 
                           IndexType nelem,
                           IndexType NSHLmax,
                           IndexType *ELMNSHL,
                           IndexType *ien);

void free_sparsity_pattern(SparsityPattern *pSparsityPattern);