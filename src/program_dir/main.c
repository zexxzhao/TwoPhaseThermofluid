
#include <util_dir/config.h>
#include <util_dir/utility.h>
#include <util_dir/mesh.h>

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
/*
 * TODO:
 * transpose the mesh data xg and IEN
 * transpose the residual to ensure DOFs belonging to the same node are contiguous
 */

double tet_vol(const double *xl) {
    double cache[9];
    for(IndexType i = 0; i < 3; ++i) {
        for(IndexType j = 0; j < 3; ++j) {
            cache[i * 3 + j] = xl[(i + 1) * 3 + j] - xl[j];
        }
    }
    double tmp = cache[0] * (cache[4] * cache[8] - cache[5] * cache[7]) - cache[1] * (cache[3] * cache[8] - cache[5] * cache[6]) + cache[2] * (cache[3] * cache[7] - cache[4] * cache[6]);
    return tmp / 6.0;
}


int main(int argc, char **argv) {
    init(&argc, &argv);
    // ConfigType *config = Config_new("config.dat");
    // Config_print(config, stdout);
    // Config_free(config);
    // getparam_();
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MeshData mesh;
    read_mesh(&mesh);
    
    double sum = 0.0;
    // integral x[0] on the whole domain
    // loop over elements
    double *xl = malloc(mesh.NSD * sizeof(double) * mesh.maxNSHL);
    printf("mesh.NELEM = %d\n", mesh.NELEM);
    for(IndexType i = 0; i < mesh.NELEM; ++i) {
        // copy coordinates of nodes of element i to xl
        for(IndexType j = 0; j < mesh.ELMNSHL[i]; ++j) {
            for(IndexType k = 0; k < mesh.NSD; ++k) {
                xl[j * mesh.NSD + k] = mesh.xg[mesh.IEN[i * mesh.maxNSHL + j] * mesh.NSD + k];
            }
        }
        // Use one quadrature point at the center of the element
        // evaluate the integral of x[0] on the element
        double xi = (xl[0+2] + xl[3+2] + xl[6+2] + xl[9+2]) / 4.0;
        sum += tet_vol(xl) * xi;
    }

    double sum_all = 0.0;
    MPI_Allreduce(&sum, &sum_all, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    if(rank == 0) {
        printf("sum = %12.11lf\n", sum_all);
    }
    if(rank == -1) {
        printf("mesh.NNODE = %d\n", mesh.NNODE);
        // print IEN
        for(IndexType i = 0; i < mesh.NELEM; ++i) {
            printf("IEN[%d] = ", i);
            for(IndexType j = 0; j < mesh.ELMNSHL[i]; ++j) {
                printf("%d ", mesh.IEN[i * mesh.maxNSHL + j]);
            }
            printf("\n");
        }
    }
    free_mesh(&mesh);
    finalize();
    return 0;
}