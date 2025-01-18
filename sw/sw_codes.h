
void sw_seq()
{
    int i,j,k;

//printf("- oryginal code - \n\n");

#pragma scop
    for (i=1; i <=N; i++){
        if(i%100==0)
            printf("%i\n", i);
        for (j=1; j <=N; j++){
            // Block S
            m1[i][j] = INT_MIN;
            for (k=1; k <=i; k++)
                m1[i][j] = MAX(m1[i][j] ,H[i-k][j] + W[k]);
            m2[i][j] = INT_MIN;
            for (k=1; k <=j; k++)
                m2[i][j] = MAX(m2[i][j] ,H[i][j-k] + W[k]);
            H[i][j] = MAX(0, MAX( H[i-1][j-1] + s(a[i], b[i]), MAX(m1[i][j], m2[i][j])));
        }
    }
#pragma endscop


}


void sw_pluto(){

//printf("- pluto [16x16x16] - \n\n");

/* We do not support C11 <threads.h>.  */
    int t1, t2, t3, t4, t5, t6, t7;
    int lb, ub, lbp, ubp, lb2, ub2;
    int lbv, ubv;
/* Start of CLooG code */
    if (N >= 1) {
        lbp=0;
        ubp=floord(N,16);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5,t6,t7)
        for (t2=lbp;t2<=ubp;t2++) {
            for (t3=0;t3<=floord(N,16);t3++) {
                for (t4=max(1,16*t2);t4<=min(N,16*t2+15);t4++) {
                    lbv=max(1,16*t3);
                    ubv=min(N,16*t3+15);
#pragma ivdep
#pragma vector always
                    for (t5=lbv;t5<=ubv;t5++) {
                        m2[t4][t5] = INT_MIN;;
                        m1[t4][t5] = INT_MIN;;
                    }
                }
            }
        }
        for (t2=0;t2<=floord(N,8);t2++) {
            lbp=max(0,ceild(16*t2-N,16));
            ubp=min(floord(N,16),t2);
#pragma omp parallel for private(lbv,ubv,t4,t5,t6,t7)
            for (t3=lbp;t3<=ubp;t3++) {
                if (t2 >= 2*t3+1) {
                    for (t4=16*t2-16*t3;t4<=min(N,16*t2-16*t3+15);t4++) {
                        for (t5=max(1,16*t3);t5<=16*t3+15;t5++) {
                            for (t6=t5+1;t6<=t4;t6++) {
                                m1[t4][t5] = MAX(m1[t4][t5] ,H[t4-(-t5+t6)][t5] + W[(-t5+t6)]);;
                            }
                            for (t6=t4+1;t6<=t4+t5;t6++) {
                                m2[t4][t5] = MAX(m2[t4][t5] ,H[t4][t5-(-t4+t6)] + W[(-t4+t6)]);;
                                m1[t4][t5] = MAX(m1[t4][t5] ,H[t4-(-t5+t6)][t5] + W[(-t5+t6)]);;
                            }
                            H[t4][t5] = MAX(0, MAX( H[t4-1][t5-1] + s(a[t4], b[t4]), MAX(m1[t4][t5], m2[t4][t5])));;
                        }
                    }
                }
                if ((N >= 2) && (t2 == 0) && (t3 == 0)) {
                    m2[1][1] = MAX(m2[1][1] ,H[1][1 -1] + W[1]);;
                    m1[1][1] = MAX(m1[1][1] ,H[1 -1][1] + W[1]);;
                    H[1][1] = MAX(0, MAX( H[1 -1][1 -1] + s(a[1], b[1]), MAX(m1[1][1], m2[1][1])));;
                    for (t5=2;t5<=min(15,N);t5++) {
                        for (t6=2;t6<=t5;t6++) {
                            m2[1][t5] = MAX(m2[1][t5] ,H[1][t5-(t6-1)] + W[(t6-1)]);;
                        }
                        m2[1][t5] = MAX(m2[1][t5] ,H[1][t5-t5] + W[t5]);;
                        m1[1][t5] = MAX(m1[1][t5] ,H[1 -1][t5] + W[1]);;
                        H[1][t5] = MAX(0, MAX( H[1 -1][t5-1] + s(a[1], b[1]), MAX(m1[1][t5], m2[1][t5])));;
                    }
                }
                if ((N == 1) && (t2 == 0) && (t3 == 0)) {
                    m2[1][1] = MAX(m2[1][1] ,H[1][1 -1] + W[1]);;
                    m1[1][1] = MAX(m1[1][1] ,H[1 -1][1] + W[1]);;
                    H[1][1] = MAX(0, MAX( H[1 -1][1 -1] + s(a[1], b[1]), MAX(m1[1][1], m2[1][1])));;
                }
                if ((t2 >= 2) && (t2 == 2*t3) && (t2 <= floord(N-1,8))) {
                    for (t6=8*t2+1;t6<=16*t2;t6++) {
                        if (t2%2 == 0) {
                            m2[8*t2][8*t2] = MAX(m2[8*t2][8*t2] ,H[8*t2][8*t2-(-8*t2+t6)] + W[(-8*t2+t6)]);;
                        }
                        if (t2%2 == 0) {
                            m1[8*t2][8*t2] = MAX(m1[8*t2][8*t2] ,H[8*t2-(-8*t2+t6)][8*t2] + W[(-8*t2+t6)]);;
                        }
                    }
                    if (t2%2 == 0) {
                        H[8*t2][8*t2] = MAX(0, MAX( H[8*t2-1][8*t2-1] + s(a[8*t2], b[8*t2]), MAX(m1[8*t2][8*t2], m2[8*t2][8*t2])));;
                    }
                    for (t5=8*t2+1;t5<=min(N,8*t2+15);t5++) {
                        for (t6=8*t2+1;t6<=t5;t6++) {
                            if (t2%2 == 0) {
                                m2[8*t2][t5] = MAX(m2[8*t2][t5] ,H[8*t2][t5-(-8*t2+t6)] + W[(-8*t2+t6)]);;
                            }
                        }
                        for (t6=t5+1;t6<=8*t2+t5;t6++) {
                            if (t2%2 == 0) {
                                m2[8*t2][t5] = MAX(m2[8*t2][t5] ,H[8*t2][t5-(-8*t2+t6)] + W[(-8*t2+t6)]);;
                            }
                            if (t2%2 == 0) {
                                m1[8*t2][t5] = MAX(m1[8*t2][t5] ,H[8*t2-(-t5+t6)][t5] + W[(-t5+t6)]);;
                            }
                        }
                        if (t2%2 == 0) {
                            H[8*t2][t5] = MAX(0, MAX( H[8*t2-1][t5-1] + s(a[8*t2], b[8*t2]), MAX(m1[8*t2][t5], m2[8*t2][t5])));;
                        }
                    }
                }
                if ((8*t2 == N) && (16*t3 == N)) {
                    for (t6=N+1;t6<=2*N;t6++) {
                        if (N%16 == 0) {
                            m2[N][N] = MAX(m2[N][N] ,H[N][N-(t6-N)] + W[(t6-N)]);;
                        }
                        if (N%16 == 0) {
                            m1[N][N] = MAX(m1[N][N] ,H[N-(t6-N)][N] + W[(t6-N)]);;
                        }
                    }
                    if (N%16 == 0) {
                        H[N][N] = MAX(0, MAX( H[N-1][N-1] + s(a[N], b[N]), MAX(m1[N][N], m2[N][N])));;
                    }
                }
                if (t2 <= 2*t3-1) {
                    for (t4=max(1,16*t2-16*t3);t4<=16*t2-16*t3+15;t4++) {
                        for (t5=16*t3;t5<=min(N,16*t3+15);t5++) {
                            for (t6=t4+1;t6<=t5;t6++) {
                                m2[t4][t5] = MAX(m2[t4][t5] ,H[t4][t5-(-t4+t6)] + W[(-t4+t6)]);;
                            }
                            for (t6=t5+1;t6<=t4+t5;t6++) {
                                m2[t4][t5] = MAX(m2[t4][t5] ,H[t4][t5-(-t4+t6)] + W[(-t4+t6)]);;
                                m1[t4][t5] = MAX(m1[t4][t5] ,H[t4-(-t5+t6)][t5] + W[(-t5+t6)]);;
                            }
                            H[t4][t5] = MAX(0, MAX( H[t4-1][t5-1] + s(a[t4], b[t4]), MAX(m1[t4][t5], m2[t4][t5])));;
                        }
                    }
                }
                if (t2 == 2*t3) {
                    for (t4=max(2,8*t2+1);t4<=min(N-1,8*t2+14);t4++) {
                        for (t5=max(1,8*t2);t5<=t4-1;t5++) {
                            for (t6=t5+1;t6<=t4;t6++) {
                                if (t2%2 == 0) {
                                    m1[t4][t5] = MAX(m1[t4][t5] ,H[t4-(-t5+t6)][t5] + W[(-t5+t6)]);;
                                }
                            }
                            for (t6=t4+1;t6<=t4+t5;t6++) {
                                if (t2%2 == 0) {
                                    m2[t4][t5] = MAX(m2[t4][t5] ,H[t4][t5-(-t4+t6)] + W[(-t4+t6)]);;
                                }
                                if (t2%2 == 0) {
                                    m1[t4][t5] = MAX(m1[t4][t5] ,H[t4-(-t5+t6)][t5] + W[(-t5+t6)]);;
                                }
                            }
                            if (t2%2 == 0) {
                                H[t4][t5] = MAX(0, MAX( H[t4-1][t5-1] + s(a[t4], b[t4]), MAX(m1[t4][t5], m2[t4][t5])));;
                            }
                        }
                        for (t6=t4+1;t6<=2*t4;t6++) {
                            if (t2%2 == 0) {
                                m2[t4][t4] = MAX(m2[t4][t4] ,H[t4][t4-(-t4+t6)] + W[(-t4+t6)]);;
                            }
                            if (t2%2 == 0) {
                                m1[t4][t4] = MAX(m1[t4][t4] ,H[t4-(-t4+t6)][t4] + W[(-t4+t6)]);;
                            }
                        }
                        if (t2%2 == 0) {
                            H[t4][t4] = MAX(0, MAX( H[t4-1][t4-1] + s(a[t4], b[t4]), MAX(m1[t4][t4], m2[t4][t4])));;
                        }
                        for (t5=t4+1;t5<=min(N,8*t2+15);t5++) {
                            for (t6=t4+1;t6<=t5;t6++) {
                                if (t2%2 == 0) {
                                    m2[t4][t5] = MAX(m2[t4][t5] ,H[t4][t5-(-t4+t6)] + W[(-t4+t6)]);;
                                }
                            }
                            for (t6=t5+1;t6<=t4+t5;t6++) {
                                if (t2%2 == 0) {
                                    m2[t4][t5] = MAX(m2[t4][t5] ,H[t4][t5-(-t4+t6)] + W[(-t4+t6)]);;
                                }
                                if (t2%2 == 0) {
                                    m1[t4][t5] = MAX(m1[t4][t5] ,H[t4-(-t5+t6)][t5] + W[(-t5+t6)]);;
                                }
                            }
                            if (t2%2 == 0) {
                                H[t4][t5] = MAX(0, MAX( H[t4-1][t5-1] + s(a[t4], b[t4]), MAX(m1[t4][t5], m2[t4][t5])));;
                            }
                        }
                    }
                }
                if ((N >= 2) && (t2 == 2*t3) && (t2 <= floord(N-1,8)) && (t2 >= ceild(N-14,8))) {
                    for (t5=max(1,8*t2);t5<=N-1;t5++) {
                        for (t6=t5+1;t6<=N;t6++) {
                            if (t2%2 == 0) {
                                m1[N][t5] = MAX(m1[N][t5] ,H[N-(-t5+t6)][t5] + W[(-t5+t6)]);;
                            }
                        }
                        for (t6=N+1;t6<=t5+N;t6++) {
                            if (t2%2 == 0) {
                                m2[N][t5] = MAX(m2[N][t5] ,H[N][t5-(t6-N)] + W[(t6-N)]);;
                            }
                            if (t2%2 == 0) {
                                m1[N][t5] = MAX(m1[N][t5] ,H[N-(-t5+t6)][t5] + W[(-t5+t6)]);;
                            }
                        }
                        if (t2%2 == 0) {
                            H[N][t5] = MAX(0, MAX( H[N-1][t5-1] + s(a[N], b[N]), MAX(m1[N][t5], m2[N][t5])));;
                        }
                    }
                    for (t6=N+1;t6<=2*N;t6++) {
                        if (t2%2 == 0) {
                            m2[N][N] = MAX(m2[N][N] ,H[N][N-(t6-N)] + W[(t6-N)]);;
                        }
                        if (t2%2 == 0) {
                            m1[N][N] = MAX(m1[N][N] ,H[N-(t6-N)][N] + W[(t6-N)]);;
                        }
                    }
                    if (t2%2 == 0) {
                        H[N][N] = MAX(0, MAX( H[N-1][N-1] + s(a[N], b[N]), MAX(m1[N][N], m2[N][N])));;
                    }
                }
                if ((t2 == 2*t3) && (t2 <= floord(N-15,8))) {
                    for (t5=max(1,8*t2);t5<=8*t2+14;t5++) {
                        for (t6=t5+1;t6<=8*t2+15;t6++) {
                            if (t2%2 == 0) {
                                m1[(8*t2+15)][t5] = MAX(m1[(8*t2+15)][t5] ,H[(8*t2+15)-(-t5+t6)][t5] + W[(-t5+t6)]);;
                            }
                        }
                        for (t6=8*t2+16;t6<=8*t2+t5+15;t6++) {
                            if (t2%2 == 0) {
                                m2[(8*t2+15)][t5] = MAX(m2[(8*t2+15)][t5] ,H[(8*t2+15)][t5-(-8*t2+t6-15)] + W[(-8*t2+t6-15)]);;
                            }
                            if (t2%2 == 0) {
                                m1[(8*t2+15)][t5] = MAX(m1[(8*t2+15)][t5] ,H[(8*t2+15)-(-t5+t6)][t5] + W[(-t5+t6)]);;
                            }
                        }
                        if (t2%2 == 0) {
                            H[(8*t2+15)][t5] = MAX(0, MAX( H[(8*t2+15)-1][t5-1] + s(a[(8*t2+15)], b[(8*t2+15)]), MAX(m1[(8*t2+15)][t5], m2[(8*t2+15)][t5])));;
                        }
                    }
                    for (t6=8*t2+16;t6<=16*t2+30;t6++) {
                        if (t2%2 == 0) {
                            m2[(8*t2+15)][(8*t2+15)] = MAX(m2[(8*t2+15)][(8*t2+15)] ,H[(8*t2+15)][(8*t2+15)-(-8*t2+t6-15)] + W[(-8*t2+t6-15)]);;
                        }
                        if (t2%2 == 0) {
                            m1[(8*t2+15)][(8*t2+15)] = MAX(m1[(8*t2+15)][(8*t2+15)] ,H[(8*t2+15)-(-8*t2+t6-15)][(8*t2+15)] + W[(-8*t2+t6-15)]);;
                        }
                    }
                    if (t2%2 == 0) {
                        H[(8*t2+15)][(8*t2+15)] = MAX(0, MAX( H[(8*t2+15)-1][(8*t2+15)-1] + s(a[(8*t2+15)], b[(8*t2+15)]), MAX(m1[(8*t2+15)][(8*t2+15)], m2[(8*t2+15)][(8*t2+15)])));;
                    }
                }
            }
        }
    }
/* End of CLooG code */
}

void sw_traco(){

//printf("- traco [16x16x16 - \n\n");

    int c1,c2,c3,c4,c5,c6,c7,c8,c9,c11,c10,c12,c13,c14,c15;

#pragma omp parallel for
    for( c1 = 0; c1 <= (N - 1)/16; c1 += 1)
        for( c3 = 0; c3 <= (N - 1) / 16; c3 += 1)
            for( c5 = 16 * c1 + 1; c5 <= min(N, 16 * c1 + 16); c5 += 1)
                for( c7 = 16 * c3 + 1; c7 <= min(N, 16 * c3 + 16); c7 += 1){
                    m1[c5][c7] = INT_MIN;
                    m2[c5][c7] = INT_MIN;
                }



    for( c1 = 0; c1 <= floord(N - 1, 8); c1 += 1)
#pragma omp parallel for schedule(dynamic, 1) shared(c1) private(c3,c4,c7,c9,c10,c5,c11)
            for( c3 = max(0, c1 - (N + 15) / 16 + 1); c3 <= min(c1, (N - 1) / 16); c3 += 1)
                for( c4 = 0; c4 <= 2; c4 += 1) {
                    if (c4 == 2) {
                        for( c7 = 16 * c1 - 16 * c3 + 1; c7 <= min(N, 16 * c1 - 16 * c3 + 16); c7 += 1)
                            for( c9 = 16 * c3 + 1; c9 <= min(N, 16 * c3 + 16); c9 += 1) {
                                for( c10 = max(0, 16 * c1 - 16 * c3 - c7 + 2); c10 <= min(1, -16 * c3 + c9 - 1); c10 += 1) {
                                    if (c10 == 1) {
                                        for( c11 = 1; c11 <= c9; c11 += 1)
                                            m2[c7][c9] = MAX(m2[c7][c9] + H[c7][c9-c11], W[c11]);
                                    } else
                                        for( c11 = 1; c11 <= c7; c11 += 1)
                                            m1[c7][c9] = MAX(m1[c7][c9] + H[c7-c11][c9], W[c11]);
                                }
                                H[c7][c9] = MAX(0, MAX(H[c7-1][c9-1] + s(a[c7],b[c9]), MAX(m1[c7][c9], m2[c7][c9])));
                            }
                    } else if (c4 == 1) {
                        for( c5 = 0; c5 <= c3; c5 += 1)
                            for( c7 = 16 * c1 - 16 * c3 + 1; c7 <= min(N, 16 * c1 - 16 * c3 + 16); c7 += 1)
                                for( c11 = 16 * c5 + 1; c11 <= min(16 * c3 + 1, 16 * c5 + 16); c11 += 1)
                                    m2[c7][(16*c3+1)] = MAX(m2[c7][(16*c3+1)], H[c7][(16*c3+1)-c11] + W[c11]);
                    } else
                        for( c5 = 0; c5 <= c1 - c3; c5 += 1)
                            for( c9 = 16 * c3 + 1; c9 <= min(N, 16 * c3 + 16); c9 += 1)
                                for( c11 = 16 * c5 + 1; c11 <= min(16 * c1 - 16 * c3 + 1, 16 * c5 + 16); c11 += 1)
                                    m1[(16*c1-16*c3+1)][c9] = MAX(m1[(16*c1-16*c3+1)][c9], H[(16*c1-16*c3+1)-c11][c9] + W[c11]);
                }



/*
for (i=1; i <=N; i++)
    for (j=1; j <=N; j++){
    // Block S
        m1[i][j] = INT_MIN;
        for (k=1; k <=i; k++)
            m1[i][j] = MAX(m1[i][j] ,H[i-k][j] + W[k]);
        m2[i][j] = INT_MIN;
        for (k=1; k <=j; k++)
            m2[i][j] = MAX(m2[i][j] ,H[i][j-k] + W[k]);
        H[i][j] = MAX(0, MAX( H[i-1][j-1] + s(a[i], b[i]), MAX(m1[i][j], m2[i][j])));
    }
*/


}

void sw_dapt(){

//printf("- dapt - \n\n");


    for (int w0 = 0; w0 <= floord(N, 8); w0 += 1) {
#pragma omp parallel for
        for (int h0 = max(0, w0 - (N + 16) / 16 + 1); h0 <= min(w0, N / 16); h0 += 1) {
            for (int i0 = max(1, 16 * h0); i0 <= min(N, 16 * h0 + 15); i0 += 1) {
                for (int i1 = max(1, 16 * w0 - 16 * h0); i1 <= min(N, 16 * w0 - 16 * h0 + 15); i1 += 1) {
                    {
                        m1[i0][i1] = (INT_MIN);
                        for (int i3 = 1; i3 <= i0; i3 += 1) {
                            m1[i0][i1] = MAX(m1[i0][i1], H[i0 - i3][i1] + W[i3]);
                        }
                        m2[i0][i1] = (INT_MIN);
                        for (int i3 = 1; i3 <= i1; i3 += 1) {
                            m2[i0][i1] = MAX(m2[i0][i1], H[i0][i1 - i3] + W[i3]);
                        }
                        H[i0][i1] = MAX(0, MAX(H[i0 - 1][i1 - 1] + s(a[i0], b[i0]), MAX(m1[i0][i1], m2[i0][i1])));
                    }
                }
            }
        }
    }

}