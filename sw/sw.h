#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <omp.h>
#include <math.h>

namespace SW {


#define min(a, b) (((a)<(b))?(a):(b))
#define MIN(a, b) (((a)<(b))?(a):(b))
#define max(a, b) (((a)>(b))?(a):(b))
#define MAX(a, b) (((a)>(b))?(a):(b))
#define floord(n, d) floor(((double)(n))/((double)(d)))
#define ceild(n, d) ceil(((double)(n))/((double)(d)))
#define CHECK_VALID 0


    int **H;
    int **H1;
    int **tmp_H;
    int **m1;
    int **m2;

    int *W;
    unsigned char *a;
    unsigned char *b;


    int DIM = N + 2;

//Similarity score of the elements that constituted the two sequences
    int s(unsigned char x, unsigned char z) {
        return (x == z) ? 1 : -1;
    }


#include "mem.h"
#include "sw_codes.h"


    void sw(){

        int i, j, k;

        int num_proc = num_threads;

        DIM = 2 * N + 2;

        // H is the scoring matrix
        H = mem();
        H1 = mem();
        m1 = mem();
        m2 = mem();

        W = (int *) malloc(DIM * sizeof(int));
        a = (unsigned char *) malloc(DIM * sizeof(unsigned char));
        b = (unsigned char *) malloc(DIM * sizeof(unsigned char));


        for (i = 0; i <= N; i++) {
            H[i][0] = 0;
            H[0][i] = 0;
            H1[i][0] = 0;
            H1[0][i] = 0;

        }


        // W is the gap alignment
        W[0] = 2;
        for (i = 1; i <= N; i++)
            W[i] = i * W[0];

        rand_seq(a, N);
        rand_seq(b, N);


        omp_set_num_threads(num_proc);

        double start = omp_get_wtime();

        if (Comp=="Pluto")
            sw_pluto();


        if (Comp=="Traco")
            sw_traco();


        if (Comp=="Dapt")
            sw_dapt();

        if (Comp=="Original") {
            tmp_H = H;
            H = H1;
            sw_seq();
            H = tmp_H;

        }

    }

}