// See the Cormen book for details of the following algorithm

// Accelerating Minimum Cost Polygon Triangulation Code with the TRACO Compiler, Palkowski Bielecki FedCsis 2018
// https://annals-csis.org/Volume_17/drp/pdf/8.pdf

#include<stdio.h>
#include<limits.h>
#include <math.h>
#include <omp.h>

namespace Triang {

    int DIM;

#define min(a, b) (((a)<(b))?(a):(b))
#define MIN(a, b) (((a)<(b))?(a):(b))
#define max(a, b) (((a)>(b))?(a):(b))
#define MAX(a, b) (((a)>(b))?(a):(b))
#define floord(n, d) floor(((double)(n))/((double)(d)))
#define ceild(n, d) ceil(((double)(n))/((double)(d)))


#include "mem.h"


    int **points;

// A utility function to find distance between two points in a plane
    double dist(int *p1, int *p2) {
        return sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) +
                    (p1[1] - p2[1]) * (p1[1] - p2[1]));
    }

// A utility function to find cost of a triangle. The cost is considered
// as perimeter (sum of lengths of all edges) of the triangle
    double cost(int i, int j, int k) {
        int *p1 = points[i];
        int *p2 = points[j];
        int *p3 = points[k];
        return dist(p1, p2) + dist(p2, p3) + dist(p3, p1);
    }


// Matrix Ai has dimension p[i-1] x p[i] for i = 1..n
    void triang() {

        int num_proc = num_threads, i;
        DIM = N + 2;

        points = (int **) malloc(DIM * sizeof(int *));

        for (i = 0; i < DIM; i++)
            points[i] = (int *) malloc(2 * sizeof(int));


        double **table = memd();


        int j, k, gap, q;


        if (Comp == "Original") {

            printf("original\n");
            for (gap = 0; gap < N; gap++) {
                if (gap % 100 == 0) printf("%i\n", gap);
                for (j = gap; j < N; j++)    // i = j - gap
                {
                    if (gap < 2)
                        table[j - gap][j] = 0.0;
                    else {
                        table[j - gap][j] = INT_MAX;
                        for (k = j - gap + 1; k < j; k++) {
                            table[j - gap][j] = MIN(table[j - gap][j],
                                                    table[j - gap][k] + table[k][j] +
                                                    cost(j - gap, j, k));

                        }
                    }
                }
            }


        }


        if (Comp == "Pluto")  // sprawdzic bez maxfuse
        {
            int t1, t2, t3, t4, t5;
            int lb, ub, lbp, ubp, lb2, ub2;
            int lbv, ubv;

            printf("pluto\n");

            if (N >= 1) {
                for (t1 = 0; t1 <= floord(N - 1, 8); t1++) {
                    lbp = ceild(t1, 2);
                    ubp = min(floord(N - 1, 16), t1);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5, t2)
                    for (t2 = lbp; t2 <= ubp; t2++) {
                        if (t1 == t2) {
                            for (t3 = 0; t3 <= min(1, N - 1); t3++) {
                                for (t4 = max(16 * t1, t3); t4 <= min(N - 1, 16 * t1 + 15); t4++) {
                                    table[t4 - t3][t4] = 0.0;;
                                }
                            }
                        }
                        for (t3 = max(2, 16 * t1 - 16 * t2);
                             t3 <= min(N - 1, 16 * t1 - 16 * t2 + 15); t3++) {
                            for (t4 = max(16 * t2, t3); t4 <= min(N - 1, 16 * t2 + 15); t4++) {
                                table[t4 - t3][t4] = INT_MAX;
                                for (t5 = -t3 + t4 + 1; t5 <= t4 - 1; t5++) {
                                    table[t4 - t3][t4] = MIN(table[t4 - t3][t4],
                                                             table[t4 - t3][t5] + table[t5][t4] +
                                                             cost(t4 - t3, t4, t5));
                                }
                            }
                        }
                    }
                }
            }


        }


        if (Comp == "Traco") {
            int c1, c3, c4, c5, c9, c11;
            int t1, t2, t3, t4, t5, t6;
            int lb, ub, lbp, ubp, lb2, ub2;
            int lbv, ubv;



            if (N >= 1) {
                lbp = 0;
                ubp = floord(N - 1, 16);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5, t6)
                for (t2 = lbp; t2 <= ubp; t2++) {
                    for (t3 = t2; t3 <= floord(N - 1, 16); t3++) {
                        if (t2 == 0) {
                            for (t4 = 0; t4 <= min(1, N - 1); t4++) {
                                lbv = max(16 * t3, t4);
                                ubv = min(N - 1, 16 * t3 + 15);
#pragma ivdep
#pragma vector always
                                for (t5 = lbv; t5 <= ubv; t5++) {
                                    table[t5 - t4][t5] = 0.0;;
                                }
                            }
                        }
                        for (t4 = max(2, 16 * t2); t4 <= min(N - 1, 16 * t2 + 15); t4++) {
                            lbv = max(16 * t3, t4);
                            ubv = min(N - 1, 16 * t3 + 15);
#pragma ivdep
#pragma vector always
                            for (t5 = lbv; t5 <= ubv; t5++) {
                                table[t5 - t4][t5] = INT_MAX;
                            }
                        }
                    }
                }

                // 1 x 32 x 16
                for (c1 = 4; c1 < 2 * N - 1; c1 += 1)
#pragma omp parallel for schedule(dynamic, 1) private(c3, c5, c9, c11) shared(c1)
                        for (c3 = -((c1 - 1) % 2) + 1;
                             c3 <= min(c1 - 4, (2 * N - c1 - 2) / 31); c3 += 2)
                            for (c5 = 0; c5 <= (c1 - c3 - 4) / 32; c5 += 1)
                                for (c9 = (c1 + 31 * c3) / 2;
                                     c9 <= min(N - 1, ((c1 + 31 * c3) / 2) + 15); c9 += 1)
                                    for (c11 = ((-c1 + c3) / 2) + 16 * c5 + c9 + 1; c11 <=
                                                                                    min(c9 - 1,
                                                                                        ((-c1 +
                                                                                          c3) / 2) +
                                                                                        16 * c5 +
                                                                                        c9 +
                                                                                        16); c11 += 1)
                                        table[c9 - ((c1 - c3) / 2)][c9] = MIN(
                                                table[c9 - ((c1 - c3) / 2)][c9],
                                                table[c9 - ((c1 - c3) / 2)][c11] + table[c11][c9] +
                                                cost(c9 - ((c1 - c3) / 2), c9, c11));


            }
        }


        if (Comp == "Dapt") {
            printf("dapt\n");
            if (N >= 3) {
                for (int w0 = 0; w0 <= (N - 1) / 8; w0 += 1) {
                    {
                        for (int i0 = 0; i0 <= 1; i0 += 1) {
                            for (int i1 = max(16 * w0, i0);
                                 i1 <= min(N - 1, 16 * w0 + 15); i1 += 1) {
                                table[-i0 + i1][i1] = 0.0;
                            }
                        }
#pragma omp parallel for
                        for (int h0 = max(0, w0 - (N + 15) / 16 + 1); h0 <= w0 / 2; h0 += 1) {
                            for (int i0 = max(2, 16 * h0);
                                 i0 <= min(N - 1, 16 * h0 + 15); i0 += 1) {
                                for (int i1 = max(16 * w0 - 16 * h0, i0);
                                     i1 <= min(N - 1, 16 * w0 - 16 * h0 + 15); i1 += 1) {
                                    {
                                        table[-i0 + i1][i1] = (INT_MAX);
                                        for (int i4 = -i0 + i1 + 1; i4 < i1; i4 += 1) {
                                            table[-i0 + i1][i1] = MIN(table[-i0 + i1][i1],
                                                                      (table[-i0 + i1][i4] +
                                                                       table[i4][i1]) +
                                                                      cost((-i0 + i1), (i1), (i4)));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                for (int i0 = 0; i0 < N; i0 += 1) {
                    for (int i1 = i0; i1 < N; i1 += 1) {
                        table[-i0 + i1][i1] = 0.0;
                    }
                }
            }


        }  // kind



        //     return table[0][N - 1];
    }


}