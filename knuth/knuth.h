#define min(a,b) (((a)<(b))?(a):(b))
#define MIN(a,b) (((a)<(b))?(a):(b))
#define max(a,b) (((a)>(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define ceild(n,d) ceil(((double)(n))/((double)(d)))


namespace Knuth{

int **c;
int **ck;
int **w;

int DIM;

#include "mem.h"

void knuth() {


    int num_proc = num_threads;
    int i, j, k, ll, p, q;
    int c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12;

    int t1, t2, t3, t4, t5, t6;
    int lb, ub, lbp, ubp, lb2, ub2;
    int lbv, ubv;

    srand(time(NULL));


    int kind = 1;


    DIM = N + 10;
    int n = N;


    omp_set_num_threads(num_proc);


    c = mem();
    ck = mem();
    w = mem();

    for (i = 0; i < DIM; i++)
        for (j = 0; j < DIM; j++) {
            ck[i][j] = i + j;
            c[i][j] = ck[i][j];
            w[i][j] = i - j;
        }


    double start = omp_get_wtime();
    if (Comp == "Original") {
        //printf("serial check\n");
#pragma scop
        for (i = n - 1; i >= 1; i--) {
            if (!(i % 100)) printf("%i\n", i);
            for (j = i + 1; j <= n; j++)
                for (k = i + 1; k < j; k++)
                    ck[i][j] = MIN(ck[i][j], w[i][j] + ck[i][k] + ck[k][j]);
        }
#pragma endscop
    }

    if (Comp == "Pluto") {// pluto
        int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
        int lb, ub, lbp, ubp, lb2, ub2;
        int lbv, ubv;
        /* Start of CLooG code */
        if (n >= 3) {
            for (t2 = -1; t2 <= floord(n - 16, 16); t2++) {
                lbp = t2 + 1;
                ubp = min(floord(n, 16), floord(16 * t2 + n + 13, 16));
#pragma omp parallel for private(lbv, ubv, t5, t6, t7, t8, t9, t10, t4) shared(t2)
                for (t4 = lbp; t4 <= ubp; t4++) {
                    for (t5 = max(max(-n + 2, 16 * t2 - 16 * t4), -16 * t4 - 13);
                         t5 <= 16 * t2 - 16 * t4 + 15; t5++) {
                        for (t7 = max(16 * t4, -t5 + 2); t7 <= min(n, 16 * t4 + 15); t7++) {
                            for (t9 = -t5 + 1; t9 <= t7 - 1; t9++) {
                                c[-t5][t7] = MIN(c[-t5][t7], w[-t5][t7] + c[-t5][t9] + c[t9][t7]);;
                            }
                        }
                    }
                }
            }
        }
        /* End of CLooG code */
    }


    if (Comp == "Traco") // traco
    {

        {
            for (c1 = 0; c1 < n + floord(-3 * n - 3, 8); c1 += 1)
#pragma omp parallel for shared(c1) private(c3, c5, c9, c11) schedule(dynamic, 1)
                    for (c3 = max(0, c1 - (n + 6) / 8 + 1);
                         c3 <= min(n / 2 - 1, c1 - (c1 + 6) / 5 + 1); c3 += 1)
                        for (c5 = 0; c5 <= c3 / 128; c5 += 1)
                            for (c7 = max(max(-n + 2 * c3 + 1, -n + 8 * c1 - 8 * c3 + 1),
                                          -n + c3 + 128 * c5 + 2);
                                 c7 <= min(-1, -n + 8 * c1 - 8 * c3 + 8); c7 += 1) {
                                if (n + 8 * c3 + c7 >= 8 * c1 + 2) {
                                    for (c11 = 256 * c5 - c7 + 1;
                                         c11 <= min(2 * c3 - c7, 256 * c5 - c7 + 256); c11 += 1)
                                        c[(-c7)][(2 * c3 - c7 + 1)] = MIN(
                                                c[(-c7)][(2 * c3 - c7 + 1)],
                                                w[(-c7)][(2 * c3 - c7 + 1)] + c[(-c7)][c11] +
                                                c[c11][(2 * c3 - c7 + 1)]);
                                    if (128 * c5 + 128 >= c3 && n + c7 >= 2 * c3 + 2) {
                                        if (c3 >= 128 * c5 + 1)
                                            for (c11 = -c7 + 1; c11 <= 256 * c5 - c7; c11 += 1)
                                                c[(-c7)][(2 * c3 - c7 + 2)] = MIN(
                                                        c[(-c7)][(2 * c3 - c7 + 2)],
                                                        w[(-c7)][(2 * c3 - c7 + 2)] +
                                                        c[(-c7)][c11] + c[c11][(2 * c3 - c7 + 2)]);
                                        for (c11 = 256 * c5 - c7 + 1; c11 <= min(2 * c3 - c7 + 1,
                                                                                 256 * c5 - c7 +
                                                                                 256); c11 += 1)
                                            c[(-c7)][(2 * c3 - c7 + 2)] = MIN(
                                                    c[(-c7)][(2 * c3 - c7 + 2)],
                                                    w[(-c7)][(2 * c3 - c7 + 2)] + c[(-c7)][c11] +
                                                    c[c11][(2 * c3 - c7 + 2)]);
                                    }
                                } else {
                                    for (c9 = max(n - 8 * c1 + 10 * c3,
                                                  n - 8 * c1 + 8 * c3 + 256 * c5 + 1);
                                         c9 <= min(n, n - 8 * c1 + 10 * c3 + 1); c9 += 1)
                                        for (c11 = n - 8 * c1 + 8 * c3 + 256 * c5; c11 <=
                                                                                   min(n - 8 * c1 +
                                                                                       8 * c3 +
                                                                                       256 * c5 +
                                                                                       255, c9 -
                                                                                            1); c11 += 1)
                                            c[(n - 8 * c1 + 8 * c3 - 1)][c9] = MIN(
                                                    c[(n - 8 * c1 + 8 * c3 - 1)][c9],
                                                    w[(n - 8 * c1 + 8 * c3 - 1)][c9] +
                                                    c[(n - 8 * c1 + 8 * c3 - 1)][c11] + c[c11][c9]);
                                }
                            }
            if ((n - 2) % 8 == 0)
                for (c5 = 0; c5 <= floord(n - 10, 256); c5 += 1)
                    for (c11 = 256 * c5 + 2; c11 <= min(n - 1, 256 * c5 + 257); c11 += 1)
                        c[1][n] = MIN(c[1][n], w[1][n] + c[1][c11] + c[c11][n]);
        }


    }

    // dapt
    if (Comp == "Dapt") {
        int h0, i0, i1, i2, w0;
        for (int w0 = -1; w0 < floord(n, 16); w0 += 1) {
#pragma omp parallel for schedule(dynamic, 1) shared(w0) private(h0, i0, i1, i2)
            for (int h0 = max(w0 - (n + 16) / 16 + 1, -((n + 13) / 16)); h0 < 0; h0 += 1) {
                for (int i0 = max(max(-n + 2, -16 * w0 + 16 * h0 - 13), 16 * h0);
                     i0 <= 16 * h0 + 15; i0 += 1) {
                    for (int i1 = max(16 * w0 - 16 * h0, -i0 + 2);
                         i1 <= min(n, 16 * w0 - 16 * h0 + 15); i1 += 1) {
                        for (int i2 = -i0 + 1; i2 < i1; i2 += 1) {
                            ck[-i0][i1] = MIN(ck[-i0][i1], (w[-i0][i1] + ck[-i0][i2]) + ck[i2][i1]);
                        }
                    }
                }
            }
        }
    }

}
}