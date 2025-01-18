

void mcc() {


    int i, j, k, ll, p, q;
    int c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12;

    int t1, t2, t3, t4, t5, t6;
    int lb, ub, lbp, ubp, lb2, ub2;
    int lbv, ubv;

    ERT = exp((float) -Ebp / (float) RT);


    srand(time(NULL));




    RNA = (char *) malloc(DIM * sizeof(char *));  //read from FASTA file
    //rand_seq(RNA, N);
    //RNA = "GGUCCAC";


    //printf("Sequence: ");
    //for(i=0; i<N; i++)
    //   printf("%c", RNA[i]);
    //printf("\n\n");




    Q = memd();
    Q1 = memd();
    Qbp = memd();
    Pbp = memd();
    Pu = memd();
    M = memd();

    rna_array_init(Q, 1, 1);
    rna_array_init(Q1, 1, 1);
    rna_array_init(Qbp, 0, 0);
    rna_array_init(Pbp, 0, 0);
    rna_array_init(Pu, 0, 0);
    rna_array_init(M, 0, 0);


    double start = omp_get_wtime();
    //  compute the partition functions Q and Qbp
    if (kind == 1) {
#pragma scop
        if (N >= 1 && l >= 0 && l <= 5)
            for (i = N - 1; i >= 0; i--) {
                if (i % 100 == 0)
                    printf("%i\n", i);
                for (j = i + 1; j <= N; j++) {
//        printf("%i\n",j );
                    Q[i][j] = Q[i][j - 1];
                    for (k = 0; k < j - i - l; k++) {
                        Qbp[k + i][j] = Q[k + i + 1][j - 1] * ERT * paired(k + i, j - 1);
                        Q[i][j] += Q[i][k + i] * Qbp[k + i][j];
                    }

                }
            }
#pragma endscop
    }


    if (kind == 2) // pluto
    {

        if ((N >= 2) && (l >= 0) && (l <= 5)) {
            for (t1 = 1; t1 <= N - 1; t1++) {
                for (t2 = 0; t2 <= floord(t1 - 1, 16); t2++) {
                    for (t3 = 0; t3 <= t2; t3++) {
                        if ((t1 >= l + 1) && (t2 == 0) && (t3 == 0)) {
                            Qbp[0 + 0][t1] = Q[0 + 0 + 1][t1 - 1] * ERT * paired(0 + 0, t1 - 1);;
                            Q[0][t1] = Q[0][t1 - 1];;
                            Q[0][t1] += Q[0][0 + 0] * Qbp[0 + 0][t1];;
                        }
                        if (t3 == 0) {
                            for (t4 = max(1, 16 * t2); t4 <= min(16 * t2 + 15, t1 - l - 1); t4++) {
                                Qbp[0 + t4][t1] =
                                        Q[0 + t4 + 1][t1 - 1] * ERT * paired(0 + t4, t1 - 1);;
                                Q[t4][t1] = Q[t4][t1 - 1];;
                                Q[t4][t1] += Q[t4][0 + t4] * Qbp[0 + t4][t1];;
                                for (t5 = 1; t5 <= min(15, t4); t5++) {
                                    Qbp[t5 + (t4 - t5)][t1] = Q[t5 + (t4 - t5) + 1][t1 - 1] * ERT *
                                                              paired(t5 + (t4 - t5), t1 - 1);;
                                    Q[(t4 - t5)][t1] +=
                                            Q[(t4 - t5)][t5 + (t4 - t5)] * Qbp[t5 + (t4 - t5)][t1];;
                                }
                            }
                        }
                        if (t3 == 0) {
                            for (t4 = max(16 * t2, t1 - l); t4 <= min(t1 - 1, 16 * t2 + 15); t4++) {
                                Q[t4][t1] = Q[t4][t1 - 1];;
                            }
                        }
                        if (t3 >= 1) {
                            for (t4 = 16 * t2; t4 <= min(16 * t2 + 15, t1 - l - 1); t4++) {
                                for (t5 = 16 * t3; t5 <= min(t4, 16 * t3 + 15); t5++) {
                                    Qbp[t5 + (t4 - t5)][t1] = Q[t5 + (t4 - t5) + 1][t1 - 1] * ERT *
                                                              paired(t5 + (t4 - t5), t1 - 1);;
                                    Q[(t4 - t5)][t1] +=
                                            Q[(t4 - t5)][t5 + (t4 - t5)] * Qbp[t5 + (t4 - t5)][t1];;
                                }
                            }
                        }
                    }
                }
            }
        }

    }
    if (kind == 3) // traco
    {


        if (N >= 10 && l >= 0 && l <= 5)
            for (c1 = 1; c1 < N + (N - 2) / 16; c1 += 1)
#pragma omp parallel for schedule(dynamic, 1)
                    for (c3 = max(0, -N + c1 + 1); c3 <= (c1 - 1) / 17; c3 += 1)
                        for (c4 = 0; c4 <= 1; c4 += 1) {
                            if (c4 == 1) {
                                for (c5 = 0; c5 <= c3; c5 += 1)
                                    for (c9 = N - c1 + 17 * c3;
                                         c9 <= min(N - 1, N - c1 + 17 * c3 + 15); c9 += 1) {
                                        if (c5 == c3 && c1 + c9 >= N + 17 * c3 + 1)
                                            Q[(N - c1 + c3 - 1)][c9] = Q[(N - c1 + c3 - 1)][c9 - 1];
                                        if (c5 == c3 && c1 + c9 >= N + 17 * c3 + 1)
                                            for (c11 = 0; c11 < 16 * c3; c11 += 1)
                                                Q[(N - c1 + c3 - 1)][c9] +=
                                                        Q[(N - c1 + c3 - 1)][c11 +
                                                                             (N - c1 + c3 - 1)] *
                                                        Qbp[c11 + (N - c1 + c3 - 1)][c9];
                                        for (c11 = 16 * c5; c11 <= min(16 * c5 + 15, -N + c1 - c3 +
                                                                                     c9); c11 += 1) {
                                            Qbp[c11 + (N - c1 + c3 - 1)][c9] =
                                                    Q[c11 + (N - c1 + c3 - 1) + 1][c9 - 1] * ERT *
                                                    paired(c11 + (N - c1 + c3 - 1), c9 - 1);
                                            if (c5 == c3) {
                                                Q[(N - c1 + c3 - 1)][c9] +=
                                                        Q[(N - c1 + c3 - 1)][c11 +
                                                                             (N - c1 + c3 - 1)] *
                                                        Qbp[c11 + (N - c1 + c3 - 1)][c9];
                                            } else if (c1 + c9 == N + 17 * c3) {
                                                Q[(N - c1 + c3 - 1)][(N - c1 + 17 * c3)] +=
                                                        Q[(N - c1 + c3 - 1)][c11 +
                                                                             (N - c1 + c3 - 1)] *
                                                        Qbp[c11 + (N - c1 + c3 - 1)][(N - c1 +
                                                                                      17 * c3)];
                                            }
                                        }
                                    }
                            } else {
                                Q[(N - c1 + c3 - 1)][(N - c1 + 17 * c3)] = Q[(N - c1 + c3 - 1)][
                                        (N - c1 + 17 * c3) - 1];
                            }
                        }


    }


    if (kind == 5) // dapt
    {

        if (l >= 0 && l <= 5) {
            if (l + 1 >= N) {
                for (int w0 = floord(-N + 2, 16); w0 <= 0; w0 += 1) {
                    for (int i0 = max(-N + 2, 16 * w0); i0 <= min(0, 16 * w0 + 15); i0 += 1) {
                        for (int i1 = -i0 + 1; i1 < N; i1 += 1) {
                            Q[-i0][i1] = Q[-i0][i1 - 1];
                        }
                    }
                }
            } else {
                for (int w0 = -1; w0 <= (N - 1) / 16; w0 += 1) {
#pragma omp parallel for
                    for (int h0 = max(w0 - (N + 15) / 16 + 1, -((N + 13) / 16));
                         h0 <= min(0, w0); h0 += 1) {
                        for (int i0 = max(max(-N + 2, -16 * w0 + 16 * h0 - 14), 16 * h0);
                             i0 <= min(0, 16 * h0 + 15); i0 += 1) {
                            for (int i1 = max(16 * w0 - 16 * h0, -i0 + 1);
                                 i1 <= min(N - 1, 16 * w0 - 16 * h0 + 15); i1 += 1) {
                                {
                                    Q[-i0][i1] = Q[-i0][i1 - 1];
                                    for (int i3 = 0; i3 < -l + i0 + i1; i3 += 1) {
                                        {
                                            Qbp[-i0 + i3][i1] = ((Q[-i0 + i3 + 1][i1 - 1] * (ERT)) *
                                                                 paired((-i0 + i3), (i1 - 1)));
                                            Q[-i0][i1] += (Q[-i0][-i0 + i3] * Qbp[-i0 + i3][i1]);
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }


    }




}
