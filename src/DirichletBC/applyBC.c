#include <stdlib.h>
#include "applyBC.h"

void MatApplyDirichletBC(IndexType n, IndexType nshl, 
                         IndexType nvar, IndexType iel,
                         const IndexType __restrict__ *ien,
                         const IndexType __restrict__ *ibc,
                         ScalarType *A) {

    // Allocate memory for local ibc array
    IndexType * ibc_loc_ptr = (IndexType *) malloc(nshl * nvar * sizeof(IndexType));
    // Copy ibc to local array
    for(IndexType ivar = 0; ivar < nvar; ++ivar) {
        for(IndexType ish = 0; ish < nshl; ++ish) {
            IndexType innode = ien[ROW_MAJOR_INDEX_2D(ish, iel, n)];
            ibc_loc_ptr[ivar * nshl + ish] = ibc[ROW_MAJOR_INDEX_2D(ivar, innode, n)];
        }
    }

    // Apply BC
    for(IndexType ish = 0; ish < nshl; ++ish) {
        for(IndexType ivar = 0; ivar < nvar; ++ivar) {
            if(! ibc_loc_ptr[ivar * nshl + ish]) continue;
            for(IndexType jsh = 0; jsh < nshl; ++jsh) {
                for(IndexType jvar = 0; jvar < nvar; ++jvar) {
                    A[ROW_MAJOR_INDEX_2D(ish, jsh, nshl) * nvar * nvar + ivar * nvar + jvar] = 0.0;
                    A[ROW_MAJOR_INDEX_2D(jsh, ish, nshl) * nvar * nvar + jvar * nvar + ivar] = 0.0;
                }
            }
            A[ROW_MAJOR_INDEX_2D(ish, ish, nshl) * nvar * nvar + ivar * nvar + ivar] = 1.0;
        }
    }

    free(ibc_loc_ptr);
}