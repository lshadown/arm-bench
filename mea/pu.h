void pu(){
    int i, j, k, ll, p, q;
    int c0, c1, c2, c3, c4, c5, c6, c7, c8, c9, c10, c11, c12, c13, c15;

    int t1, t2, t3, t4, t5, t6, t7;
    int lb, ub, lbp, ubp, lb2, ub2;
    int lbv, ubv;




    if(kind==1){
#pragma scop
        for(i=N-1; i>=0; i--){
            for(j=i+1; j<N; j++){
                Pu[i][j] = (Q[0][i]*Q[j][N-1]*1)/Q[0][N-1];
                for(p=0; p<i; p++){
                    for(q=j+1; q<N; q++){
                        Pu[i][j] += (Pbp[p][q] * ERT * Q[p+1][i] * 1 * Q[j+1][q-1]) /  (Qbp[p][q] ==0 ? 1 : Qbp[p][q]) ;
                    }
                }
            }
        }
#pragma endscop

    }
    if(kind==2) // pluto
    {
        printf("pluto\n");

/* We do not support C11 <threads.h>.  */
        int t1, t2, t3, t4, t5, t6, t7, t8;
        int lb, ub, lbp, ubp, lb2, ub2;

/* Start of CLooG code */
        if (N >= 2) {
            lbp=0;
            ubp=floord(N-2,16);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
            for (t2=lbp;t2<=ubp;t2++) {
                for (t3=t2;t3<=floord(N-1,16);t3++) {
                    for (t4=16*t2;t4<=min(min(N-2,16*t2+15),16*t3+14);t4++) {
                        lbv=max(16*t3,t4+1);
                        ubv=min(N-1,16*t3+15);
#pragma ivdep
#pragma vector always
                        for (t5=lbv;t5<=ubv;t5++) {
                            Pu[t4][t5] = (Q[0][t4]*Q[t5][N-1]*1)/Q[0][N-1];;
                        }
                    }
                }
            }
            if (N >= 4) {
                lbp=0;
                ubp=floord(N-3,16);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7,t8)
                for (t2=lbp;t2<=ubp;t2++) {
                    for (t3=t2;t3<=floord(N-2,16);t3++) {
                        for (t4=0;t4<=min(floord(N-4,16),t2);t4++) {
                            for (t5=max(16*t2,16*t4+1);t5<=min(min(N-3,16*t2+15),16*t3+14);t5++) {
                                for (t6=max(16*t3,t5+1);t6<=min(N-2,16*t3+15);t6++) {
                                    for (t7=16*t4;t7<=min(16*t4+15,t5-1);t7++) {
                                        for (t8=t6+1;t8<=N-1;t8++) {
                                            Pu[t5][t6] += (Pbp[t7][t8] * ERT * Q[t7+1][t5] * 1 * Q[t6+1][t8-1]) / (Qbp[t7][t8] ==0 ? 1 :Qbp[t7][t8]) ;;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
/* End of CLooG code */

    }
    if(kind==3) // traco
    {
        printf("traco\n");


    }

    if(kind==5)
    {
        if (N >= 2) {
            for (int w0 = max(-1, -((N + 13) / 16)); w0 <= (N - 1) / 16; w0 += 1) {
                {
#pragma omp parallel for
                    for (int h0 = -((N + 12) / 16); h0 <= w0 - (N + 15) / 16; h0 += 1) {
                        for (int h1 = max(w0 + 1, -h0 - 1); h1 <= (N - 2) / 16; h1 += 1) {
                            for (int i0 = max(max(-N + 3, 16 * h0), -16 * h1 - 14); i0 <= 16 * h0 + 15; i0 += 1) {
                                for (int i1 = max(16 * h1, -i0 + 1); i1 <= min(N - 2, 16 * h1 + 15); i1 += 1) {
                                    for (int i3 = 16 * w0 - 16 * h0 - 16 * h1; i3 <= min(16 * w0 - 16 * h0 - 16 * h1 + 15, -i0 - 1); i3 += 1) {
                                        for (int i4 = i1 + 1; i4 < N; i4 += 1) {
                                            Pu[-i0][i1] += (((((Pbp[i3][i4] * (ERT)) * Q[i3 + 1][-i0]) * 1) * Q[i1 + 1][i4 - 1]) / ((Qbp[i3][i4] == 0) ? 1 : Qbp[i3][i4]));
                                        }
                                    }
                                }
                            }
                        }
                    }
#pragma omp parallel for
                    for (int h0 = max(w0 - (N + 15) / 16 + 1, -((N + 13) / 16)); h0 <= min(0, w0); h0 += 1) {
                        {
                            for (int h1 = max(w0 + 1, -h0 - 1); h1 < w0 - h0; h1 += 1) {
                                for (int i0 = max(16 * h0, -16 * h1 - 14); i0 <= 16 * h0 + 15; i0 += 1) {
                                    for (int i1 = max(16 * h1, -i0 + 1); i1 <= 16 * h1 + 15; i1 += 1) {
                                        for (int i3 = 16 * w0 - 16 * h0 - 16 * h1; i3 <= min(16 * w0 - 16 * h0 - 16 * h1 + 15, -i0 - 1); i3 += 1) {
                                            for (int i4 = i1 + 1; i4 < N; i4 += 1) {
                                                Pu[-i0][i1] += (((((Pbp[i3][i4] * (ERT)) * Q[i3 + 1][-i0]) * 1) * Q[i1 + 1][i4 - 1]) / ((Qbp[i3][i4] == 0) ? 1 : Qbp[i3][i4]));
                                            }
                                        }
                                    }
                                }
                            }
                            for (int i0 = max(max(-N + 2, -16 * w0 + 16 * h0 - 14), 16 * h0); i0 <= min(0, 16 * h0 + 15); i0 += 1) {
                                for (int i1 = max(16 * w0 - 16 * h0, -i0 + 1); i1 <= min(N - 1, 16 * w0 - 16 * h0 + 15); i1 += 1) {
                                    {
                                        Pu[-i0][i1] = (((Q[0][-i0] * Q[i1][N - 1]) * 1) / Q[0][N - 1]);
                                        for (int i3 = 0; i3 <= min(15, -i0 - 1); i3 += 1) {
                                            for (int i4 = i1 + 1; i4 < N; i4 += 1) {
                                                Pu[-i0][i1] += (((((Pbp[i3][i4] * (ERT)) * Q[i3 + 1][-i0]) * 1) * Q[i1 + 1][i4 - 1]) / ((Qbp[i3][i4] == 0) ? 1 : Qbp[i3][i4]));
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


}