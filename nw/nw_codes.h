void nw_dapt(){



    for (int w0 = 0; w0 <= floord(N, 8); w0 += 1) {
#pragma omp parallel for
        for (int h0 = max(0, w0 - (N + 16) / 16 + 1); h0 <= min(w0, N / 16); h0 += 1) {
            for (int i0 = max(1, 16 * h0); i0 <= min(N, 16 * h0 + 15); i0 += 1) {
                for (int i1 = max(1, 16 * w0 - 16 * h0); i1 <= min(N, 16 * w0 - 16 * h0 + 15); i1 += 1) {
                    {
                        m1[i0][i1] = (INT_MIN);
                        for (int i3 = 1; i3 <= i0; i3 += 1) {
                            m1[i0][i1] = MAX(m1[i0][i1], F[i0 - i3][i1] - W[i3]);
                        }
                        m2[i0][i1] = (INT_MIN);
                        for (int i3 = 1; i3 <= i1; i3 += 1) {
                            m2[i0][i1] = MAX(m2[i0][i1], F[i0][i1 - i3] - W[i3]);
                        }
                        F[i0][i1] = MAX(0, MAX(F[i0 - 1][i1 - 1] + s(a[i0], b[i0]), MAX(m1[i0][i1], m2[i0][i1])));
                    }
                }
            }
        }
    }

}



void nw_pluto(){

    int t1, t2, t3, t4, t5, t6, t7;
    int lb, ub, lbp, ubp, lb2, ub2;
    int lbv, ubv;
    if (N >= 1) {
        lbp=0;
        ubp=floord(N,71);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7)
        for (t2=lbp;t2<=ubp;t2++) {
            for (t3=0;t3<=floord(N,177);t3++) {
                for (t4=max(1,71*t2);t4<=(min(N,71*t2+70))-7;t4+=8) {
                    lbv=max(1,177*t3);
                    ubv=min(N,177*t3+176);
#pragma ivdep
#pragma vector always
                    for (t5=lbv;t5<=ubv;t5++) {
                        m2[t4][t5] = INT_MIN;;
                        m1[t4][t5] = INT_MIN;;
                        m2[(t4+1)][t5] = INT_MIN;;
                        m1[(t4+1)][t5] = INT_MIN;;
                        m2[(t4+2)][t5] = INT_MIN;;
                        m1[(t4+2)][t5] = INT_MIN;;
                        m2[(t4+3)][t5] = INT_MIN;;
                        m1[(t4+3)][t5] = INT_MIN;;
                        m2[(t4+4)][t5] = INT_MIN;;
                        m1[(t4+4)][t5] = INT_MIN;;
                        m2[(t4+5)][t5] = INT_MIN;;
                        m1[(t4+5)][t5] = INT_MIN;;
                        m2[(t4+6)][t5] = INT_MIN;;
                        m1[(t4+6)][t5] = INT_MIN;;
                        m2[(t4+7)][t5] = INT_MIN;;
                        m1[(t4+7)][t5] = INT_MIN;;
                    }
                }
                for (;t4<=min(N,71*t2+70);t4++) {
                    lbv=max(1,177*t3);
                    ubv=min(N,177*t3+176);
#pragma ivdep
#pragma vector always
                    for (t5=lbv;t5<=ubv;t5++) {
                        m2[t4][t5] = INT_MIN;;
                        m1[t4][t5] = INT_MIN;;
                    }
                }
            }
        }
        for (t2=0;t2<=floord(248*N,12567);t2++) {
            lbp=max(0,ceild(71*t2-N,71));
            ubp=min(floord(N,177),t2);
#pragma omp parallel for private(lbv,ubv,t4,t5,t6,t7)
            for (t3=lbp;t3<=ubp;t3++) {
                if ((N >= 2) && (t2 == 0) && (t3 == 0)) {
                    m2[1][1] = MAX(m2[1][1] ,F[1][1 -1] - W[1]);;
                    m1[1][1] = MAX(m1[1][1] ,F[1 -1][1] - W[1]);;
                    F[1][1] = MAX(0, MAX( F[1 -1][1 -1] + s(a[1], b[1]), MAX(m1[1][1], m2[1][1])));;
                    for (t5=2;t5<=min(176,N);t5++) {
                        for (t6=2;t6<=t5;t6++) {
                            m2[1][t5] = MAX(m2[1][t5] ,F[1][t5-(t6-1)] - W[(t6-1)]);;
                        }
                        m2[1][t5] = MAX(m2[1][t5] ,F[1][t5-t5] - W[t5]);;
                        m1[1][t5] = MAX(m1[1][t5] ,F[1 -1][t5] - W[1]);;
                        F[1][t5] = MAX(0, MAX( F[1 -1][t5-1] + s(a[1], b[1]), MAX(m1[1][t5], m2[1][t5])));;
                    }
                }
                if ((N == 1) && (t2 == 0) && (t3 == 0)) {
                    m2[1][1] = MAX(m2[1][1] ,F[1][1 -1] - W[1]);;
                    m1[1][1] = MAX(m1[1][1] ,F[1 -1][1] - W[1]);;
                    F[1][1] = MAX(0, MAX( F[1 -1][1 -1] + s(a[1], b[1]), MAX(m1[1][1], m2[1][1])));;
                }
                for (t4=max(1,71*t2-71*t3);t4<=min(177*t3-1,71*t2-71*t3+70);t4++) {
                    for (t5=177*t3;t5<=min(N,177*t3+176);t5++) {
                        for (t6=t4+1;t6<=t5;t6++) {
                            m2[t4][t5] = MAX(m2[t4][t5] ,F[t4][t5-(-t4+t6)] - W[(-t4+t6)]);;
                        }
                        for (t6=t5+1;t6<=t4+t5;t6++) {
                            m2[t4][t5] = MAX(m2[t4][t5] ,F[t4][t5-(-t4+t6)] - W[(-t4+t6)]);;
                            m1[t4][t5] = MAX(m1[t4][t5] ,F[t4-(-t5+t6)][t5] - W[(-t5+t6)]);;
                        }
                        F[t4][t5] = MAX(0, MAX( F[t4-1][t5-1] + s(a[t4], b[t4]), MAX(m1[t4][t5], m2[t4][t5])));;
                    }
                }
                if ((t2 <= floord(248*t3,71)) && (t2 >= ceild(248*t3-70,71)) && (t3 >= 1) && (t3 <= floord(N-1,177))) {
                    for (t6=177*t3+1;t6<=354*t3;t6++) {
                        m2[177*t3][177*t3] = MAX(m2[177*t3][177*t3] ,F[177*t3][177*t3-(-177*t3+t6)] - W[(-177*t3+t6)]);;
                        m1[177*t3][177*t3] = MAX(m1[177*t3][177*t3] ,F[177*t3-(-177*t3+t6)][177*t3] - W[(-177*t3+t6)]);;
                    }
                    F[177*t3][177*t3] = MAX(0, MAX( F[177*t3-1][177*t3-1] + s(a[177*t3], b[177*t3]), MAX(m1[177*t3][177*t3], m2[177*t3][177*t3])));;
                    for (t5=177*t3+1;t5<=min(N,177*t3+176);t5++) {
                        for (t6=177*t3+1;t6<=t5;t6++) {
                            m2[177*t3][t5] = MAX(m2[177*t3][t5] ,F[177*t3][t5-(-177*t3+t6)] - W[(-177*t3+t6)]);;
                        }
                        for (t6=t5+1;t6<=177*t3+t5;t6++) {
                            m2[177*t3][t5] = MAX(m2[177*t3][t5] ,F[177*t3][t5-(-177*t3+t6)] - W[(-177*t3+t6)]);;
                            m1[177*t3][t5] = MAX(m1[177*t3][t5] ,F[177*t3-(-t5+t6)][t5] - W[(-t5+t6)]);;
                        }
                        F[177*t3][t5] = MAX(0, MAX( F[177*t3-1][t5-1] + s(a[177*t3], b[177*t3]), MAX(m1[177*t3][t5], m2[177*t3][t5])));;
                    }
                }
                if ((t2 >= ceild(248*N-12390,12567)) && (177*t3 == N)) {
                    for (t6=N+1;t6<=2*N;t6++) {
                        if (176*N%177 == 0) {
                            m2[N][N] = MAX(m2[N][N] ,F[N][N-(t6-N)] - W[(t6-N)]);;
                        }
                        if (176*N%177 == 0) {
                            m1[N][N] = MAX(m1[N][N] ,F[N-(t6-N)][N] - W[(t6-N)]);;
                        }
                    }
                    if (N%177 == 0) {
                        F[N][N] = MAX(0, MAX( F[N-1][N-1] + s(a[N], b[N]), MAX(m1[N][N], m2[N][N])));;
                    }
                }
                for (t4=max(max(2,71*t2-71*t3),177*t3+1);t4<=min(min(N-1,177*t3+175),71*t2-71*t3+70);t4++) {
                    for (t5=max(1,177*t3);t5<=t4-1;t5++) {
                        for (t6=t5+1;t6<=t4;t6++) {
                            m1[t4][t5] = MAX(m1[t4][t5] ,F[t4-(-t5+t6)][t5] - W[(-t5+t6)]);;
                        }
                        for (t6=t4+1;t6<=t4+t5;t6++) {
                            m2[t4][t5] = MAX(m2[t4][t5] ,F[t4][t5-(-t4+t6)] - W[(-t4+t6)]);;
                            m1[t4][t5] = MAX(m1[t4][t5] ,F[t4-(-t5+t6)][t5] - W[(-t5+t6)]);;
                        }
                        F[t4][t5] = MAX(0, MAX( F[t4-1][t5-1] + s(a[t4], b[t4]), MAX(m1[t4][t5], m2[t4][t5])));;
                    }
                    for (t6=t4+1;t6<=2*t4;t6++) {
                        m2[t4][t4] = MAX(m2[t4][t4] ,F[t4][t4-(-t4+t6)] - W[(-t4+t6)]);;
                        m1[t4][t4] = MAX(m1[t4][t4] ,F[t4-(-t4+t6)][t4] - W[(-t4+t6)]);;
                    }
                    F[t4][t4] = MAX(0, MAX( F[t4-1][t4-1] + s(a[t4], b[t4]), MAX(m1[t4][t4], m2[t4][t4])));;
                    for (t5=t4+1;t5<=min(N,177*t3+176);t5++) {
                        for (t6=t4+1;t6<=t5;t6++) {
                            m2[t4][t5] = MAX(m2[t4][t5] ,F[t4][t5-(-t4+t6)] - W[(-t4+t6)]);;
                        }
                        for (t6=t5+1;t6<=t4+t5;t6++) {
                            m2[t4][t5] = MAX(m2[t4][t5] ,F[t4][t5-(-t4+t6)] - W[(-t4+t6)]);;
                            m1[t4][t5] = MAX(m1[t4][t5] ,F[t4-(-t5+t6)][t5] - W[(-t5+t6)]);;
                        }
                        F[t4][t5] = MAX(0, MAX( F[t4-1][t5-1] + s(a[t4], b[t4]), MAX(m1[t4][t5], m2[t4][t5])));;
                    }
                }
                if ((t2 <= floord(248*t3+176,71)) && (t2 >= ceild(248*t3+106,71)) && (t3 <= floord(N-176,177))) {
                    for (t5=max(1,177*t3);t5<=177*t3+175;t5++) {
                        for (t6=t5+1;t6<=177*t3+176;t6++) {
                            m1[(177*t3+176)][t5] = MAX(m1[(177*t3+176)][t5] ,F[(177*t3+176)-(-t5+t6)][t5] - W[(-t5+t6)]);;
                        }
                        for (t6=177*t3+177;t6<=177*t3+t5+176;t6++) {
                            m2[(177*t3+176)][t5] = MAX(m2[(177*t3+176)][t5] ,F[(177*t3+176)][t5-(-177*t3+t6-176)] - W[(-177*t3+t6-176)]);;
                            m1[(177*t3+176)][t5] = MAX(m1[(177*t3+176)][t5] ,F[(177*t3+176)-(-t5+t6)][t5] - W[(-t5+t6)]);;
                        }
                        F[(177*t3+176)][t5] = MAX(0, MAX( F[(177*t3+176)-1][t5-1] + s(a[(177*t3+176)], b[(177*t3+176)]), MAX(m1[(177*t3+176)][t5], m2[(177*t3+176)][t5])));;
                    }
                    for (t6=177*t3+177;t6<=354*t3+352;t6++) {
                        m2[(177*t3+176)][(177*t3+176)] = MAX(m2[(177*t3+176)][(177*t3+176)] ,F[(177*t3+176)][(177*t3+176)-(-177*t3+t6-176)] - W[(-177*t3+t6-176)]);;
                        m1[(177*t3+176)][(177*t3+176)] = MAX(m1[(177*t3+176)][(177*t3+176)] ,F[(177*t3+176)-(-177*t3+t6-176)][(177*t3+176)] - W[(-177*t3+t6-176)]);;
                    }
                    F[(177*t3+176)][(177*t3+176)] = MAX(0, MAX( F[(177*t3+176)-1][(177*t3+176)-1] + s(a[(177*t3+176)], b[(177*t3+176)]), MAX(m1[(177*t3+176)][(177*t3+176)], m2[(177*t3+176)][(177*t3+176)])));;
                }
                if ((N >= 2) && (t2 >= ceild(71*t3+N-70,71)) && (t3 <= floord(N-1,177)) && (t3 >= ceild(N-175,177))) {
                    for (t5=max(1,177*t3);t5<=N-1;t5++) {
                        for (t6=t5+1;t6<=N;t6++) {
                            m1[N][t5] = MAX(m1[N][t5] ,F[N-(-t5+t6)][t5] - W[(-t5+t6)]);;
                        }
                        for (t6=N+1;t6<=t5+N;t6++) {
                            m2[N][t5] = MAX(m2[N][t5] ,F[N][t5-(t6-N)] - W[(t6-N)]);;
                            m1[N][t5] = MAX(m1[N][t5] ,F[N-(-t5+t6)][t5] - W[(-t5+t6)]);;
                        }
                        F[N][t5] = MAX(0, MAX( F[N-1][t5-1] + s(a[N], b[N]), MAX(m1[N][t5], m2[N][t5])));;
                    }
                    for (t6=N+1;t6<=2*N;t6++) {
                        m2[N][N] = MAX(m2[N][N] ,F[N][N-(t6-N)] - W[(t6-N)]);;
                        m1[N][N] = MAX(m1[N][N] ,F[N-(t6-N)][N] - W[(t6-N)]);;
                    }
                    F[N][N] = MAX(0, MAX( F[N-1][N-1] + s(a[N], b[N]), MAX(m1[N][N], m2[N][N])));;
                }
                for (t4=max(71*t2-71*t3,177*t3+177);t4<=min(N,71*t2-71*t3+70);t4++) {
                    for (t5=max(1,177*t3);t5<=177*t3+176;t5++) {
                        for (t6=t5+1;t6<=t4;t6++) {
                            m1[t4][t5] = MAX(m1[t4][t5] ,F[t4-(-t5+t6)][t5] - W[(-t5+t6)]);;
                        }
                        for (t6=t4+1;t6<=t4+t5;t6++) {
                            m2[t4][t5] = MAX(m2[t4][t5] ,F[t4][t5-(-t4+t6)] - W[(-t4+t6)]);;
                            m1[t4][t5] = MAX(m1[t4][t5] ,F[t4-(-t5+t6)][t5] - W[(-t5+t6)]);;
                        }
                        F[t4][t5] = MAX(0, MAX( F[t4-1][t5-1] + s(a[t4], b[t4]), MAX(m1[t4][t5], m2[t4][t5])));;
                    }
                }
            }
        }
    }
}




void nw_traco()
{
    int c1,c3,c4,c5,c7,c9,c11,c10;
    for( c1 = 0; c1 <= floord(N - 1, 8); c1 += 1)
#pragma omp parallel for schedule(dynamic, 1) shared(c1) private(c3,c4,c7,c9,c10,c5,c11)
            for( c3 = max(0, c1 - (N + 15) / 16 + 1); c3 <= min(c1, (N - 1) / 16); c3 += 1)
                for( c4 = 0; c4 <= 4; c4 += 1) {
                    if (c4 == 4) {
                        for( c7 = 16 * c1 - 16 * c3 + 1; c7 <= min(N, 16 * c1 - 16 * c3 + 16); c7 += 1)
                            for( c9 = 16 * c3 + 1; c9 <= min(N, 16 * c3 + 16); c9 += 1) {
                                if (16 * c3 + c7 >= 16 * c1 + 2)
                                    for( c11 = 1; c11 <= c7; c11 += 1)
                                        m1[c7][c9] = MAX(m1[c7][c9] ,F[c7-c11][c9] - W[c11]);
                                for( c10 = max(3, 16 * c3 - c9 + 5); c10 <= 4; c10 += 1) {
                                    if (c10 == 4) {
                                        F[c7][c9] = MAX(0, MAX( F[c7-1][c9-1] + s(a[c7], b[c7]), MAX(m1[c7][c9], m2[c7][c9])));
                                    } else {
                                        for( c11 = 1; c11 <= c9; c11 += 1)
                                            m2[c7][c9] = MAX(m2[c7][c9] ,F[c7][c9-c11] - W[c11]);
                                    }
                                }
                            }
                    } else if (c4 == 3) {
                        for( c5 = 0; c5 <= c3; c5 += 1)
                            for( c7 = 16 * c1 - 16 * c3 + 1; c7 <= min(N, 16 * c1 - 16 * c3 + 16); c7 += 1)
                                for( c11 = 16 * c5 + 1; c11 <= min(16 * c3 + 1, 16 * c5 + 16); c11 += 1)
                                    m2[c7][(16*c3+1)] = MAX(m2[c7][(16*c3+1)] ,F[c7][(16*c3+1)-c11] - W[c11]);
                    } else if (c4 == 2) {
                        for( c7 = 16 * c1 - 16 * c3 + 1; c7 <= min(N, 16 * c1 - 16 * c3 + 16); c7 += 1)
                            for( c9 = 16 * c3 + 1; c9 <= min(N, 16 * c3 + 16); c9 += 1)
                                m2[c7][c9] = INT_MIN;
                    } else if (c4 == 1) {
                        for( c5 = 0; c5 <= c1 - c3; c5 += 1)
                            for( c9 = 16 * c3 + 1; c9 <= min(N, 16 * c3 + 16); c9 += 1)
                                for( c11 = 16 * c5 + 1; c11 <= min(16 * c1 - 16 * c3 + 1, 16 * c5 + 16); c11 += 1)
                                    m1[(16*c1-16*c3+1)][c9] = MAX(m1[(16*c1-16*c3+1)][c9] ,F[(16*c1-16*c3+1)-c11][c9] - W[c11]);
                    } else {
                        for( c7 = 16 * c1 - 16 * c3 + 1; c7 <= min(N, 16 * c1 - 16 * c3 + 16); c7 += 1)
                            for( c9 = 16 * c3 + 1; c9 <= min(N, 16 * c3 + 16); c9 += 1)
                                m1[c7][c9] = INT_MIN;
                    }
                }




}


void nw_seq()
{
    int i,j,k;

    printf("- oryginal code - \n\n");

#pragma scop
    for (i=1; i <=N; i++)
        for (j=1; j <=N; j++){
            // Block S
            m1[i][j] = INT_MIN;
            for (k=1; k <=i; k++)
                m1[i][j] = MAX(m1[i][j] ,F[i-k][j] - W[k]);
            m2[i][j] = INT_MIN;
            for (k=1; k <=j; k++)
                m2[i][j] = MAX(m2[i][j] ,F[i][j-k] - W[k]);
            F[i][j] = MAX(0, MAX( F[i-1][j-1] + s(a[i], b[i]), MAX(m1[i][j], m2[i][j])));
        }
#pragma endscop


}