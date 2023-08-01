#progma once

typedef IndexType int;
typedef ScalarType double;

#define ROW_MAJOR_INDEX_2D(i, j, n) (i * n + j)
#define COL_MAJOR_INDEX_2D(i, j, n) (j * n + i)

#define ROW_MAJOR_INDEX_4D(i, j, k, l, n, m, o) (i * n * m * o + j * m * o + k * o + l)
#define COL_MAJOR_INDEX_4D(i, j, k, l, n, m, o) (l * n * m * o + k * n * m + j * n + i)

/**
 * @brief Apply Dirichlet boundary conditions to a matrix
 * 
 * @param n Number of nodes
 * @param nshl Number of shape functions
 * @param nvar Number of variables
 * @param iel Index of element
 * @param ien ien[*][n] is the connectivity array
 * @param ibc ibc[*][n] is the boundary condition array
 * @param A A[nshl][nshl][nvar][nvar] is the matrix to which the boundary conditions are applied
 */
void MatApplyDirichletBC(IndexType n, IndexType nshl,
                         IndexType nvar, IndexType iel,
                         const IndexType __restrict__ *ien,
                         const IndexType __restrict__ *ibc,
                         ScalarType *A);