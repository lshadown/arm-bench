#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <math.h>

namespace Zuker {

#define min(a, b) (((a)<(b))?(a):(b))
#define MIN(a, b) (((a)<(b))?(a):(b))
#define max(a, b) (((a)>(b))?(a):(b))
#define floord(n, d) floor(((double)(n))/((double)(d)))
#define ceild(n, d) ceil(((double)(n))/((double)(d)))

#define CHECK_VALID 0

    int **W;
    int **V;
    int **V1;
    int **tmp_V;
    int **W1;
    int **tmp_W;
    int **EFL;
    int **EHF;
    int **VMI;
    int **VBI;
    int **WZ;

    int DIM = N+2;

#include "zuker_codes.h"
#include "mem.h"


    void zuker() {

        int i, j, k;
//

        int num_proc = num_threads;

        DIM = N + 2;


        W = mem();
        V = mem();
        V1 = mem();
        W1 = mem();
        EFL = mem();
        EHF = mem();
        VMI = mem();
        VBI = mem();
        WZ = mem();

        for (i = 0; i < N; i++)
            for (j = 0; j < N; j++) {
                W[i][j] = i * j;
                V[i][j] = i + 1;
                EHF[i][j] = i + 1;
                EFL[i][j] = i + 1;
                V1[i][j] = V[i][j];
                W1[i][j] = W[i][j];
            }
        omp_set_num_threads(num_proc);

        double start = omp_get_wtime();

        if (Comp == "Dapt") {
            printf("dapt\n");
            zuker_dapt();
        }


        if (Comp == "Traco") {
            printf("traco tilecorr\n");
            zuker_traco();
        }

        if (Comp == "Pluto") {
            printf("pluto\n");
            zuker_pluto();
        }

        if (Comp == "Original") {
            tmp_V = V;
            V = V1;
            tmp_W = W;
            W = W1;

            zuker_seq();
            V = tmp_V;
            W = tmp_W;

        }

    }

}