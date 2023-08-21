#include "utility.h"
#include <stdio.h>
#include <string.h>
#include <math.h>


#include <mpi.h>

#define NRES (4)

void init(int *argv, char ***argc) {
#if defined(HAVE_PETSC)
    PetscInitialize(argv, argc, NULL, NULL);
#elif defined(HAVE_MPI)
    MPI_Init(argv, argc);
#endif
}

void finalize() {
#if defined(HAVE_PETSC)
    PetscFinalize();
#elif defined(HAVE_MPI)
    MPI_Finalize();
#endif
}

void calculate_residual(ScalarType *pResidual, 
                        const ScalarType *pRHSGU,
                        const ScalarType *pRHSGP,
                        const ScalarType *pRHSGLS,
                        const ScalarType *pRHSGTem,
                        IndexType nNodes,
                        IndexType nDim) {

    memset(pResidual, 0, sizeof(ScalarType) * NRES);
    for(IndexType i = 0; i < nNodes * nDim; ++i) {
        pResidual[0] += pRHSGU[i] * pRHSGU[i];
    }
    for(IndexType i = 0; i < nNodes; ++i) {
        pResidual[1] += pRHSGP[i] * pRHSGP[i];
    }
    for(IndexType i = 0; i < nNodes; ++i) {
        pResidual[2] += pRHSGLS[i] * pRHSGLS[i];
    }
    for(IndexType i = 0; i < nNodes; ++i) {
        pResidual[3] += pRHSGTem[i] * pRHSGTem[i];
    }
    ScalarType tmp[NRES];
    memset(tmp, 0, sizeof(ScalarType) * NRES);
    int mpi_err = MPI_Allreduce(pResidual, tmp, NRES, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(mpi_err != MPI_SUCCESS) {
        printf("MPI_Allreduce failed in calculate_residual\n");
        MPI_Abort(MPI_COMM_WORLD, mpi_err);
    }
    memcpy(pResidual, tmp, sizeof(ScalarType) * NRES);
    for(IndexType i = 0; i < NRES; ++i) {
        pResidual[i] = sqrt(pResidual[i]);
    }
}

void check_convergence(IndexType *pConverged, 
                       const ScalarType *pResidual,
                       const ScalarType *pResidual0,
                       const ScalarType *pResidualTol,
                       IndexType assembleFieldFlag) {
    memset(pConverged, 0, sizeof(IndexType)* NRES);
    if(assembleFieldFlag & ASSEMBLE_FIELD_NS) {
        if(pResidual[0] < pResidual0[0] * pResidualTol[0] ||
            pResidual[0] < 1e-15) {
            pConverged[0] = 1;
        }
        if(pResidual[1] < pResidual0[1] * pResidualTol[1] ||
            pResidual[1] < 1e-15) {
            pConverged[1] = 1;
        }
    }
    if(assembleFieldFlag & ASSEMBLE_FIELD_LS) {
        if(pResidual[2] < pResidual0[2] * pResidualTol[2] ||
            pResidual[2] < 1e-15) {
            pConverged[2] = 1;
        }
    }
    if(assembleFieldFlag & ASSEMBLE_FIELD_TEM) {
        if(pResidual[3] < pResidual0[3] * pResidualTol[3] ||
            pResidual[3] < 1e-15) {
            pConverged[3] = 1;
        }
    }
}


void print_residual(const ScalarType *pResidual,
                     const ScalarType *pResidual0,
                     const ScalarType *pResidualTol,
                     IndexType assembleFieldFlag,
                     IndexType iter) {
    int mpi_rank;
    int mpi_err = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    if(mpi_err != MPI_SUCCESS) {
        printf("MPI_Comm_rank failed in print_residual\n");
        MPI_Abort(MPI_COMM_WORLD, mpi_err);
    }
    if(mpi_rank != 0) {
        return;
    }

    if(assembleFieldFlag & ASSEMBLE_FIELD_NS) {
        printf("%d) MOM: abs = %e,  rel = %e, tol = %e\n", 
               iter, pResidual[0], pResidual[0] / MAX(pResidual0[0], 1e-15), pResidualTol[0]);
        printf("%d) CON: abs = %e,  rel = %e, tol = %e\n",
               iter, pResidual[1], pResidual[1] / MAX(pResidual0[1], 1e-15), pResidualTol[1]);
    }
    if(assembleFieldFlag & ASSEMBLE_FIELD_LS) {
        printf("%d) VOF: abs = %e,  rel = %e, tol = %e\n",
               iter, pResidual[2], pResidual[2] / MAX(pResidual0[2], 1e-15), pResidualTol[2]);
    }
    if(assembleFieldFlag & ASSEMBLE_FIELD_TEM) {
        printf("%d) TEM: abs = %e,  rel = %e, tol = %e\n",
               iter, pResidual[3], pResidual[3] / MAX(pResidual0[3], 1e-15), pResidualTol[3]);
    }
}




void FORTRAN_NAME(calculate_residual)(ScalarType *pResidual, 
                                      const ScalarType *pRHSGU,
                                      const ScalarType *pRHSGP,
                                      const ScalarType *pRHSGLS,
                                      const ScalarType *pRHSGTem,
                                      IndexType *nNodes,
                                      IndexType *nDim) {
    calculate_residual(pResidual, pRHSGU, pRHSGP, pRHSGLS, pRHSGTem, *nNodes, *nDim);
}

void FORTRAN_NAME(check_convergence)(IndexType *pConverged, 
                                     const ScalarType *pResidual,
                                     const ScalarType *pResidual0,
                                     const ScalarType *pResidualTol,
                                     IndexType *assembleFieldFlag) {
    check_convergence(pConverged, pResidual, pResidual0, pResidualTol, *assembleFieldFlag);
}

void FORTRAN_NAME(print_residual)(const ScalarType *pResidual,
                                    const ScalarType *pResidual0,
                                    const ScalarType *pResidualTol,
                                    IndexType *assembleFieldFlag,
                                    IndexType *iter) {
    print_residual(pResidual, pResidual0, pResidualTol, *assembleFieldFlag, *iter);
}
