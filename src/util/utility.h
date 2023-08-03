#pragma once

#include "../sys/defines.h"

void calculate_residual(ScalarType *pResidual, 
                        const ScalarType *pRHSGU,
                        const ScalarType *pRHSGP,
                        const ScalarType *pRHSGLS,
                        const ScalarType *pRHSGTem,
                        IndexType nNodes,
                        IndexType nDim);


void check_convergence(IndexType *pConverged,
                       const ScalarType *pResidual,
                       const ScalarType *pResidual0,
                       const ScalarType *pResidualTol,
                       IndexType assembleFieldFlag);


void print_residual(const ScalarType *pResidual,
                    const ScalarType *pResidual0,
                    const ScalarType *pResidualTol,
                    IndexType assembleFieldFlag,
                    IndexType iter);

/* Fortran function prototypes*/

void FORTRAN_NAME(calculate_residual)(ScalarType *pResidual, 
                                      const ScalarType *pRHSGU,
                                      const ScalarType *pRHSGP,
                                      const ScalarType *pRHSGLS,
                                      const ScalarType *pRHSGTem,
                                      IndexType *nNodes,
                                      IndexType *nDim);

void FORTRAN_NAME(check_convergence)(IndexType *pConverged,
                                        const ScalarType *pResidual,
                                        const ScalarType *pResidual0,
                                        const ScalarType *pResidualTol,
                                        IndexType *assembleFieldFlag);

void FORTRAN_NAME(print_residual)(const ScalarType *pResidual,
                                  const ScalarType *pResidual0,
                                  const ScalarType *pResidualTol,
                                  IndexType *assembleFieldFlag,
                                  IndexType *iter);