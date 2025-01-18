void pb() {


    int i, j, k, ll, p, q;
    int c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c15;

    int t1, t2, t3, t4, t5, t6, t7;
    int lb, ub, lbp, ubp, lb2, ub2;
    int lbv, ubv;



    if (kind == 1) {
#pragma scop
        for (i = 0; i < N; i++) {
            if (i % 50 == 0) printf("%i\n", i);
            for (j = i + 1; j < N; j++) {
                Pbp[i][j] = (Q[0][i] * Q[j][N - 1] * Qbp[i][j]) /
                            Q[0][N - 1];   //  Pbp[i][j] = (Q[1][i]*Q[j+1][N]*Qbp[i][j])/Q[0][N-1];
                for (p = 0; p < i; p++) {
                    for (q = j + 1; q < N; q++) {
                        Pbp[i][j] += (Pbp[p][q] * ERT * Q[p + 1][i] * Qbp[i][j] * Q[j + 1][q - 1]) /
                                     (Qbp[p][q] == 0 ? 1 : Qbp[p][q]);

                    }
                }
            }
        }
#pragma endscop
    }
    if (kind == 2) // pluto
    {
        printf("pluto\n");
        lbp = 0;
        ubp = floord(N - 2, 16);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5, t6, t7)
        for (t2 = lbp; t2 <= ubp; t2++) {
            for (t3 = t2; t3 <= floord(N - 1, 16); t3++) {
                for (t4 = 16 * t2; t4 <= min(min(N - 2, 16 * t2 + 15), 16 * t3 + 14); t4++) {
                    lbv = max(16 * t3, t4 + 1);
                    ubv = min(N - 1, 16 * t3 + 15);
#pragma ivdep
#pragma vector always
                    for (t5 = lbv; t5 <= ubv; t5++) {
                        Pbp[t4][t5] = (Q[0][t4] * Q[t5][N - 1] * Qbp[t4][t5]) / Q[0][N - 1];;
                    }
                }
            }
        }
/*
  for (t2=0;t2<=floord(N-4,8);t2++) {
    lbp=max(0,ceild(16*t2-N+3,16));
    ubp=floord(t2,2);
#pragma omp parallel for private(lbv,ubv,t4,t5,t6,t7)
    for (t3=lbp;t3<=ubp;t3++) {
      for (t4=max(16*t2-16*t3,16*t3+1);t4<=min(N-3,16*t2-16*t3+15);t4++) {
        for (t5=16*t3;t5<=min(16*t3+15,t4-1);t5++) {
          for (t6=t4+1;t6<=N-2;t6++) {
            for (t7=t6+1;t7<=N-1;t7++) {
              Pbp[t4][t6] += (Pbp[t5][t7] * ERT * Q[t5+1][t4] * Qbp[t4][t6] * Q[t6+1][t7-1]) / (Qbp[t5][t7]  ==0 ? 1 : Qbp[t5][t7]);
            }
          }
        }
      }
    }
  }
*/
        if (N >= 4) {
            for (t1 = 1; t1 <= floord(17 * N - 52, 16); t1++) {
                lbp = max(0, t1 - N + 3);
                ubp = floord(t1 - 1, 17);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5, t6)
                for (t2 = lbp; t2 <= ubp; t2++) {
                    for (t4 = 16 * t2; t4 <= min(16 * t2 + 15, t1 - t2 - 1); t4++) {
                        for (t5 = t1 - t2 + 1; t5 <= N - 2; t5++) {
                            for (t6 = t5 + 1; t6 <= N - 1; t6++) {
                                Pbp[(t1 - t2)][t5] += (Pbp[t4][t6] * ERT * Q[t4 + 1][(t1 - t2)] *
                                                       Qbp[(t1 - t2)][t5] * Q[t5 + 1][t6 - 1]) /
                                                      (Qbp[t4][t6] == 0 ? 1 : Qbp[t4][t6]);
                            }
                        }
                    }
                }
            }
        }


    }
    if (kind == 3) // traco
    {
        printf("traco\n");
        lbp = 0;
        ubp = floord(N - 2, 16);
#pragma omp parallel for private(lbv, ubv, t3, t4, t5, t6, t7)
        for (t2 = lbp; t2 <= ubp; t2++) {
            for (t3 = t2; t3 <= floord(N - 1, 16); t3++) {
                for (t4 = 16 * t2; t4 <= min(min(N - 2, 16 * t2 + 15), 16 * t3 + 14); t4++) {
                    lbv = max(16 * t3, t4 + 1);
                    ubv = min(N - 1, 16 * t3 + 15);
#pragma ivdep
#pragma vector always
                    for (t5 = lbv; t5 <= ubv; t5++) {
                        Pbp[t4][t5] = (Q[0][t4] * Q[t5][N - 1] * Qbp[t4][t5]) / Q[0][N - 1];;
                    }
                }
            }
        }

        for (c1 = 1; c1 < N - 2; c1 += 1)
#pragma omp parallel for schedule(dynamic, 1)
                for (c3 = 0; c3 <= (N - c1 - 3) / 16; c3 += 1)
                    for (c5 = 0; c5 <= (c1 - 1) / 16; c5 += 1)
                        for (c7 = 0; c7 <= -c3 + (N - c1 - 3) / 16; c7 += 1)
                            for (c11 = c1 + 16 * c3 + 1;
                                 c11 <= min(c1 + 16 * c3 + 16, N - 16 * c7 - 2); c11 += 1) {
                                if (N >= 16 * c7 + c11 + 18) {
                                    for (c15 = 16 * c7 + c11 + 1;
                                         c15 <= 16 * c7 + c11 + 16; c15 += 1)
                                        Pbp[c1][c11] +=
                                                (Pbp[16 * c5][c15] * ERT * Q[16 * c5 + 1][c1] *
                                                 Qbp[c1][c11] * Q[c11 + 1][c15 - 1]) /
                                                (Qbp[16 * c5][c15] == 0 ? 1 : Qbp[16 * c5][c15]);
                                } else {
                                    for (c13 = 16 * c5;
                                         c13 <= min(c1 - 1, 16 * c5 + 15); c13 += 1) {
                                        if (c13 >= 16 * c5 + 1)
                                            for (c15 = c11 + 1; c15 <= 16 * c7 + c11; c15 += 1)
                                                Pbp[c1][c11] +=
                                                        (Pbp[c13][c15] * ERT * Q[c13 + 1][c1] *
                                                         Qbp[c1][c11] * Q[c11 + 1][c15 - 1]) /
                                                        (Qbp[c13][c15] == 0 ? 1 : Qbp[c13][c15]);
                                        for (c15 = 16 * c7 + c11 + 1; c15 < N; c15 += 1)
                                            Pbp[c1][c11] += (Pbp[c13][c15] * ERT * Q[c13 + 1][c1] *
                                                             Qbp[c1][c11] * Q[c11 + 1][c15 - 1]) /
                                                            (Qbp[c13][c15] == 0 ? 1
                                                                                : Qbp[c13][c15]);
                                    }
                                }
                            }


    }
    if (kind == 5) // dapt
    {

        for (int w0 = floord(-N + 1, 16); w0 < 0; w0 += 1) {
#pragma omp parallel for
            for (int h0 = 0; h0 <= w0 + (N - 2) / 16 + 1; h0 += 1) {
                {
                    for (int h1 = max(w0 - 2 * h0, -((N + 13) / 16)); h1 < w0 - h0; h1 += 1) {
                        for (int i0 = max(16 * h0, 16 * w0 - 16 * h0 - 16 * h1 + 1);
                             i0 <= 16 * h0 + 15; i0 += 1) {
                            for (int i1 = -16 * h1 - 15; i1 <= min(N - 2, -16 * h1); i1 += 1) {
                                for (int i3 = 16 * w0 - 16 * h0 - 16 * h1;
                                     i3 <= min(16 * w0 - 16 * h0 - 16 * h1 + 15, i0 - 1); i3 += 1) {
                                    for (int i4 = i1 + 1; i4 < N; i4 += 1) {
                                        Pbp[i0][i1] += (((((Pbp[i3][i4] * (ERT)) * Q[i3 + 1][i0]) *
                                                          Qbp[i0][i1]) * Q[i1 + 1][i4 - 1]) /
                                                        ((Qbp[i3][i4] == 0) ? 1 : Qbp[i3][i4]));
                                    }
                                }
                            }
                        }
                    }
                    for (int i0 = 16 * h0; i0 <= min(N - 2, 16 * h0 + 15); i0 += 1) {
                        for (int i1 = max(-16 * w0 + 16 * h0 - 15, i0 + 1);
                             i1 <= min(N - 1, -16 * w0 + 16 * h0); i1 += 1) {
                            {
                                Pbp[i0][i1] = (((Q[0][i0] * Q[i1][N - 1]) * Qbp[i0][i1]) /
                                               Q[0][N - 1]);
                                for (int i3 = 0; i3 <= min(15, i0 - 1); i3 += 1) {
                                    for (int i4 = i1 + 1; i4 < N; i4 += 1) {
                                        Pbp[i0][i1] += (((((Pbp[i3][i4] * (ERT)) * Q[i3 + 1][i0]) *
                                                          Qbp[i0][i1]) * Q[i1 + 1][i4 - 1]) /
                                                        ((Qbp[i3][i4] == 0) ? 1 : Qbp[i3][i4]));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        for (int w0 = 0; w0 < floord(N - 4, 16); w0 += 1) {
#pragma omp parallel for
            for (int h0 = w0 + 1; h0 <= (N - 3) / 16; h0 += 1) {
                for (int h1 = max(w0 - 2 * h0, -((N + 13) / 16)); h1 < -h0; h1 += 1) {
                    for (int i0 = max(16 * h0, 16 * w0 - 16 * h0 - 16 * h1 + 1);
                         i0 <= min(N - 3, 16 * h0 + 15); i0 += 1) {
                        for (int i1 = max(-16 * h1 - 15, i0 + 1);
                             i1 <= min(N - 2, -16 * h1); i1 += 1) {
                            for (int i3 = 16 * w0 - 16 * h0 - 16 * h1;
                                 i3 <= min(16 * w0 - 16 * h0 - 16 * h1 + 15, i0 - 1); i3 += 1) {
                                for (int i4 = i1 + 1; i4 < N; i4 += 1) {
                                    Pbp[i0][i1] += (((((Pbp[i3][i4] * (ERT)) * Q[i3 + 1][i0]) *
                                                      Qbp[i0][i1]) * Q[i1 + 1][i4 - 1]) /
                                                    ((Qbp[i3][i4] == 0) ? 1 : Qbp[i3][i4]));
                                }
                            }
                        }
                    }
                }
            }
        }


    }
}

