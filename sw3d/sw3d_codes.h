void sw_seq()
{
    int i,j,k,l;

    printf("- oryginal code - \n\n");

#pragma scop
    for (i=1; i <=N; i++)
        for (j=1; j <=N; j++){
            for (l=1; l<=N; l++){
                // Block S
                m1[i][j][l] = INT_MIN;
                for (k=1; k <=i; k++)
                    m1[i][j][l] = MAX(m1[i][j][l] ,H[i-k][j][l] - 2*W[k]);
                m2[i][j][l] = INT_MIN;
                for (k=1; k <=j; k++)
                    m2[i][j][l] = MAX(m2[i][j][l], H[i][j-k][l] - 2*W[k]);
                m3[i][j][l] = INT_MIN;
                for (k=1; k <=l; k++)
                    m3[i][j][l] = MAX(m3[i][j][l], H[i][j][l-k] - 2*W[k]);
                m4[i][j][l] = INT_MIN;
                for (k=1; k <=min(i,j); k++)
                    m4[i][j][l] = MAX(m4[i][j][l], H[i-k][j-k][l] - W[k] + s(a[i], b[j]));
                m5[i][j][l] = INT_MIN;
                for (k=1; k <=min(j,l); k++)
                    m5[i][j][l] = MAX(m5[i][j][l], H[i][j-k][l-k] - W[k] + s(b[j], c[l]));
                m6[i][j][l] = INT_MIN;
                for (k=1; k <=min(i,l); k++)
                    m6[i][j][l] = MAX(m6[i][j][l], H[i-k][j][l-k] - W[k] + s(a[i], c[l]));
                H[i][j][l] = MAX(0, MAX( H[i-1][j-1][l-1] + s(a[i], b[j]) + s(a[i], c[l]) + s(b[j], c[l]), MAX(m1[i][j][l], MAX(m2[i][j][l], MAX(m3[i][j][l], MAX(m4[i][j][l], MAX(m5[i][j][l], m6[i][j][l])))))));


            }
        }
#pragma endscop

}


void sw3d_dapt(){

    for (int w0 = 0; w0 <= floord(3 * N, 16); w0 += 1) {
#pragma omp parallel for
        for (int h0 = max(0, w0 - (N + 8) / 8 + 1); h0 <= min(w0, N / 16); h0 += 1) {
            for (int h1 = max(0, w0 - h0 - (N + 16) / 16 + 1); h1 <= min(w0 - h0, N / 16); h1 += 1) {
                for (int i0 = max(1, 16 * h0); i0 <= min(N, 16 * h0 + 15); i0 += 1) {
                    for (int i1 = max(1, 16 * h1); i1 <= min(N, 16 * h1 + 15); i1 += 1) {
                        for (int i2 = max(1, 16 * w0 - 16 * h0 - 16 * h1); i2 <= min(N, 16 * w0 - 16 * h0 - 16 * h1 + 15); i2 += 1) {
                            {
                                m1[i0][i1][i2] = (INT_MIN);
                                for (int i4 = 1; i4 <= i0; i4 += 1) {
                                    m1[i0][i1][i2] = MAX(m1[i0][i1][i2], H[i0 - i4][i1][i2] - (2 * W[i4]));
                                }
                                m2[i0][i1][i2] = (INT_MIN);
                                for (int i4 = 1; i4 <= i1; i4 += 1) {
                                    m2[i0][i1][i2] = MAX(m2[i0][i1][i2], H[i0][i1 - i4][i2] - (2 * W[i4]));
                                }
                                m3[i0][i1][i2] = (INT_MIN);
                                for (int i4 = 1; i4 <= i2; i4 += 1) {
                                    m3[i0][i1][i2] = MAX(m3[i0][i1][i2], H[i0][i1][i2 - i4] - (2 * W[i4]));
                                }
                                m4[i0][i1][i2] = (INT_MIN);
                                for (int i4 = 1; i4 <= min(i0, i1); i4 += 1) {
                                    m4[i0][i1][i2] = MAX(m4[i0][i1][i2], (H[i0 - i4][i1 - i4][i2] - W[i4]) + s(a[i0], b[i1]));
                                }
                                if (i1 >= i0 + 1) {
                                    m5[i0][i1][i2] = (INT_MIN);
                                } else {
                                    m5[i0][i1][i2] = (INT_MIN);
                                }
                                for (int i4 = 1; i4 <= min(i1, i2); i4 += 1) {
                                    m5[i0][i1][i2] = MAX(m5[i0][i1][i2], (H[i0][i1 - i4][i2 - i4] - W[i4]) + s(b[i1], c[i2]));
                                }
                                if (i2 >= i1 + 1) {
                                    m6[i0][i1][i2] = (INT_MIN);
                                } else {
                                    m6[i0][i1][i2] = (INT_MIN);
                                }
                                for (int i4 = 1; i4 <= min(i0, i2); i4 += 1) {
                                    m6[i0][i1][i2] = MAX(m6[i0][i1][i2], (H[i0 - i4][i1][i2 - i4] - W[i4]) + s(a[i0], c[i2]));
                                }
                                if (i2 >= i0 + 1) {
                                    H[i0][i1][i2] = MAX(0, MAX(((H[i0 - 1][i1 - 1][i2 - 1] + s(a[i0], b[i1])) + s(a[i0], c[i2])) + s(b[i1], c[i2]), MAX(m1[i0][i1][i2], MAX(m2[i0][i1][i2], MAX(m3[i0][i1][i2], MAX(m4[i0][i1][i2], MAX(m5[i0][i1][i2], m6[i0][i1][i2])))))));
                                } else {
                                    H[i0][i1][i2] = MAX(0, MAX(((H[i0 - 1][i1 - 1][i2 - 1] + s(a[i0], b[i1])) + s(a[i0], c[i2])) + s(b[i1], c[i2]), MAX(m1[i0][i1][i2], MAX(m2[i0][i1][i2], MAX(m3[i0][i1][i2], MAX(m4[i0][i1][i2], MAX(m5[i0][i1][i2], m6[i0][i1][i2])))))));
                                }
                            }
                        }
                    }
                }
            }
        }
    }


}

void sw_traco3d(){

    printf("- traco [16x16x16x16] - \n\n");

    int c1,c2,c3,c4,c5,c6,c7,c8,c9,c11,c10,c12,c13,c14,c15;

#pragma omp parallel for shared(N) private(c1, c2,c3,c4,c5,c6,c7,c8,c9,c11,c10,c12,c13,c14,c15)
    for( c1 = 0; c1 <= (N - 1)/16; c1 += 1)
        for( c3 = 0; c3 <= (N - 1) / 16; c3 += 1)
            for( c5 = 0; c5 <= (N - 1) / 16; c5 += 1)
                for( c7 = 16 * c1 + 1; c7 <= min(N, 16 * c1 + 16); c7 += 1)
                    for( c9 = 16 * c3 + 1; c9 <= min(N, 16 * c3 + 16); c9 += 1)
                        for( c11 = 16 * c5 + 1; c11 <= min(N, 16 * c5 + 16); c11 += 1)

                        {
                            m1[c7][c9][c11] = INT_MIN;
                            m2[c7][c9][c11] = INT_MIN;
                            m3[c7][c9][c11] = INT_MIN;
                            m4[c7][c9][c11] = INT_MIN;
                            m5[c7][c9][c11] = INT_MIN;
                            m6[c7][c9][c11] = INT_MIN;
                        }



    for( c1 = 0; c1 <= floord(N - 1, 8); c1 += 1)
#pragma omp parallel for shared(c1, N) private(c2,c3,c4,c5,c6,c7,c8,c9,c11,c10,c12,c13,c14,c15) schedule(dynamic, 1)
            for( c3 = max(0, c1 - (N + 15) / 16 + 1); c3 <= min(c1, (N - 1) / 16); c3 += 1)
                for( c5 = 0; c5 <= (N - 1) / 16; c5 += 1)
                    for( c6 = 0; c6 <= 6; c6 += 1) {
                        if (c6 == 6) {
                            for( c9 = 16 * c1 - 16 * c3 + 1; c9 <= min(N, 16 * c1 - 16 * c3 + 16); c9 += 1)
                                for( c11 = 16 * c3 + 1; c11 <= min(N, 16 * c3 + 16); c11 += 1)
                                    for( c13 = 16 * c5 + 1; c13 <= min(N, 16 * c5 + 16); c13 += 1) {
                                        for( c14 = max(0, 16 * c1 - 16 * c3 - c9 + 2); c14 <= min(1, -16 * c3 + c11 - 1); c14 += 1) {
                                            if (c14 == 1) {
                                                for( c15 = 1; c15 <= c11; c15 += 1)
                                                    m2[c9][c11][c13] = MAX(m2[c9][c11][c13], H[c9][c11-c15][c13] - 2*W[c15]);
                                            } else
                                                for( c15 = 1; c15 <= c9; c15 += 1)
                                                    m1[c9][c11][c13] = MAX(m1[c9][c11][c13] ,H[c9-c15][c11][c13] - 2*W[c15]);
                                        }
                                        for( c14 = max(2, 16 * c5 - c13 + 4); c14 <= min(min(3, -16 * c1 + 16 * c3 + c9 + 1), -16 * c3 + c11 + 1); c14 += 1) {
                                            if (c14 == 3) {
                                                for( c15 = 1; c15 <= min(c9, c11); c15 += 1)
                                                    m4[c9][c11][c13] = MAX(m4[c9][c11][c13], H[c9-c15][c11-c15][c13] - W[c15] + s(c9,c11));
                                            } else
                                                for( c15 = 1; c15 <= c13; c15 += 1)
                                                    m3[c9][c11][c13] = MAX(m3[c9][c11][c13], H[c9][c11][c13-c15] - 2*W[c15]);
                                        }
                                        if (c13 >= 16 * c5 + 2)
                                            for( c14 = max(4, 16 * c3 - c11 + 6); c14 <= min(5, -16 * c1 + 16 * c3 + c9 + 3); c14 += 1) {
                                                if (c14 == 5) {
                                                    for( c15 = 1; c15 <= min(c9, c13); c15 += 1)
                                                        m6[c9][c11][c13] = MAX(m6[c9][c11][c13], H[c9-c15][c11][c13-c15] - W[c15] + s(c9,c13));
                                                } else
                                                    for( c15 = 1; c15 <= min(c11, c13); c15 += 1)
                                                        m5[c9][c11][c13] = MAX(m5[c9][c11][c13], H[c9][c11-c15][c13-c15] - W[c15] + s(c11,c13));
                                            }
                                        H[c9][c11][c13] = MAX(0, MAX( H[c9-1][c11-1][c13-1] + s(a[c9], b[c11]) + s(a[c9], c[c13]) + s(b[c11], c[13]), MAX(m1[c9][c11][c13], MAX(m2[c9][c11][c13], MAX(m3[c9][c11][c13], MAX(m4[c9][c11][c13], MAX(m5[c9][c11][c13], m6[c9][c11][c13])))))));
                                    }
                        } else if (c6 == 5) {
                            for( c7 = 0; c7 <= min(c1 - c3, c5); c7 += 1)
                                for( c9 = 16 * c1 - 16 * c3 + 1; c9 <= min(N, 16 * c1 - 16 * c3 + 16); c9 += 1)
                                    for( c11 = 16 * c3 + 1; c11 <= min(N, 16 * c3 + 16); c11 += 1) {
                                        if (16 * c3 + c9 >= 16 * c1 + 2) {
                                            for( c15 = 16 * c7 + 1; c15 <= min(min(16 * c5 + 1, 16 * c7 + 16), c9); c15 += 1)
                                                m6[c9][c11][(16*c5+1)] = MAX(m6[c9][c11][(16*c5+1)], H[c9-c15][c11][(16*c5+1)-c15] - W[c15] + s(c9,16*c5+1));
                                        } else
                                            for( c13 = 16 * c5 + 1; c13 <= min(N, 16 * c5 + 16); c13 += 1)
                                                for( c15 = 16 * c7 + 1; c15 <= min(min(16 * c1 - 16 * c3 + 1, 16 * c7 + 16), c13); c15 += 1)
                                                    m6[(16*c1-16*c3+1)][c11][c13] = MAX(m6[(16*c1-16*c3+1)][c11][c13], H[(16*c1-16*c3+1)-c15][c11][c13-c15] - W[c15] + s(16*c1-16*c3+1,c13));
                                    }
                        } else if (c6 == 4) {
                            for( c7 = 0; c7 <= min(c3, c5); c7 += 1)
                                for( c9 = 16 * c1 - 16 * c3 + 1; c9 <= min(N, 16 * c1 - 16 * c3 + 16); c9 += 1)
                                    for( c11 = 16 * c3 + 1; c11 <= min(N, 16 * c3 + 16); c11 += 1) {
                                        if (c11 >= 16 * c3 + 2) {
                                            for( c15 = 16 * c7 + 1; c15 <= min(min(16 * c5 + 1, 16 * c7 + 16), c11); c15 += 1)
                                                m5[c9][c11][(16*c5+1)] = MAX(m5[c9][c11][(16*c5+1)], H[c9][c11-c15][(16*c5+1)-c15] - W[c15] + s(c11,16*c5+1));
                                        } else
                                            for( c13 = 16 * c5 + 1; c13 <= min(N, 16 * c5 + 16); c13 += 1)
                                                for( c15 = 16 * c7 + 1; c15 <= min(min(16 * c3 + 1, 16 * c7 + 16), c13); c15 += 1)
                                                    m5[c9][(16*c3+1)][c13] = MAX(m5[c9][(16*c3+1)][c13], H[c9][(16*c3+1)-c15][c13-c15] - W[c15] + s(16*c3+1,c13));
                                    }
                        } else if (c6 == 3) {
                            for( c7 = 0; c7 <= min(c3, c1 - c3); c7 += 1)
                                for( c9 = 16 * c1 - 16 * c3 + 1; c9 <= min(N, 16 * c1 - 16 * c3 + 16); c9 += 1) {
                                    if (16 * c3 + c9 >= 16 * c1 + 2) {
                                        for( c13 = 16 * c5 + 1; c13 <= min(N, 16 * c5 + 16); c13 += 1)
                                            for( c15 = 16 * c7 + 1; c15 <= min(min(16 * c3 + 1, 16 * c7 + 16), c9); c15 += 1)
                                                m4[c9][(16*c3+1)][c13] = MAX(m4[c9][(16*c3+1)][c13], H[c9-c15][(16*c3+1)-c15][c13] - W[c15] + s(c9,16*c3+1));
                                    } else
                                        for( c11 = 16 * c3 + 1; c11 <= min(N, 16 * c3 + 16); c11 += 1)
                                            for( c13 = 16 * c5 + 1; c13 <= min(N, 16 * c5 + 16); c13 += 1)
                                                for( c15 = 16 * c7 + 1; c15 <= min(min(16 * c1 - 16 * c3 + 1, 16 * c7 + 16), c11); c15 += 1)
                                                    m4[(16*c1-16*c3+1)][c11][c13] = MAX(m4[(16*c1-16*c3+1)][c11][c13], H[(16*c1-16*c3+1)-c15][c11-c15][c13] - W[c15] + s(16*c1-16*c3+1, c11));
                                }
                        } else if (c6 == 2) {
                            for( c7 = 0; c7 <= c5; c7 += 1)
                                for( c9 = 16 * c1 - 16 * c3 + 1; c9 <= min(N, 16 * c1 - 16 * c3 + 16); c9 += 1)
                                    for( c11 = 16 * c3 + 1; c11 <= min(N, 16 * c3 + 16); c11 += 1)
                                        for( c15 = 16 * c7 + 1; c15 <= min(16 * c5 + 1, 16 * c7 + 16); c15 += 1)
                                            m3[c9][c11][(16*c5+1)] = MAX(m3[c9][c11][(16*c5+1)], H[c9][c11][(16*c5+1)-c15] - 2*W[c15]);
                        } else if (c6 == 1) {
                            for( c7 = 0; c7 <= c3; c7 += 1)
                                for( c9 = 16 * c1 - 16 * c3 + 1; c9 <= min(N, 16 * c1 - 16 * c3 + 16); c9 += 1)
                                    for( c13 = 16 * c5 + 1; c13 <= min(N, 16 * c5 + 16); c13 += 1)
                                        for( c15 = 16 * c7 + 1; c15 <= min(16 * c3 + 1, 16 * c7 + 16); c15 += 1)
                                            m2[c9][(16*c3+1)][c13] = MAX(m2[c9][(16*c3+1)][c13], H[c9][(16*c3+1)-c15][c13] - 2*W[c15]);
                        } else
                            for( c7 = 0; c7 <= c1 - c3; c7 += 1)
                                for( c11 = 16 * c3 + 1; c11 <= min(N, 16 * c3 + 16); c11 += 1)
                                    for( c13 = 16 * c5 + 1; c13 <= min(N, 16 * c5 + 16); c13 += 1)
                                        for( c15 = 16 * c7 + 1; c15 <= min(16 * c1 - 16 * c3 + 1, 16 * c7 + 16); c15 += 1)
                                            m1[(16*c1-16*c3+1)][c11][c13] = MAX(m1[(16*c1-16*c3+1)][c11][c13] ,H[(16*c1-16*c3+1)-c15][c11][c13] - 2*W[c15]);
                    }
}

void sw3d_pluto(){
    printf("pluto\n");
    int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12, t13, t14;
    int lb, ub, lbp, ubp, lb2, ub2;
    int lbv, ubv;
/* Start of CLooG code */
    if (N >= 1) {
        for (t2=0;t2<=floord(N,8);t2++) {
            lbp=max(0,ceild(16*t2-N,16));
            ubp=min(floord(N,16),t2);
#pragma omp parallel for private(lbv,ubv,t5,t6,t7,t8,t9,t10,t11,t12,t13,t14)
            for (t4=lbp;t4<=ubp;t4++) {
                for (t6=0;t6<=floord(N,16);t6++) {
                    for (t7=max(1,16*t2-16*t4);t7<=min(N,16*t2-16*t4+15);t7++) {
                        for (t9=max(1,16*t4);t9<=min(N,16*t4+15);t9++) {
                            for (t11=max(1,16*t6);t11<=min(N,16*t6+15);t11++) {
                                m1[t7][t9][t11] = INT_MIN;;
                                for (t13=1;t13<=t7;t13++) {
                                    m1[t7][t9][t11] = MAX(m1[t7][t9][t11] ,H[t7-t13][t9][t11] - 2*W[t13]);;
                                }
                                m2[t7][t9][t11] = INT_MIN;;
                                for (t13=1;t13<=t9;t13++) {
                                    m2[t7][t9][t11] = MAX(m2[t7][t9][t11], H[t7][t9-t13][t11] - 2*W[t13]);;
                                }
                                m3[t7][t9][t11] = INT_MIN;;
                                for (t13=1;t13<=t11;t13++) {
                                    m3[t7][t9][t11] = MAX(m3[t7][t9][t11], H[t7][t9][t11-t13] - 2*W[t13]);;
                                }
                                m4[t7][t9][t11] = INT_MIN;;
                                for (t13=1;t13<=min(t7,t9);t13++) {
                                    m4[t7][t9][t11] = MAX(m4[t7][t9][t11], H[t7-t13][t9-t13][t11] - W[t13] + s(a[t7], b[t9]));;
                                }
                                m5[t7][t9][t11] = INT_MIN;;
                                for (t13=1;t13<=min(t11,t9);t13++) {
                                    m5[t7][t9][t11] = MAX(m5[t7][t9][t11], H[t7][t9-t13][t11-t13] - W[t13] + s(b[t9], c[t11]));;
                                }
                                m6[t7][t9][t11] = INT_MIN;;
                                for (t13=1;t13<=min(t11,t7);t13++) {
                                    m6[t7][t9][t11] = MAX(m6[t7][t9][t11], H[t7-t13][t9][t11-t13] - W[t13] + s(a[t7], c[t11]));;
                                }
                                H[t7][t9][t11] = MAX(0, MAX( H[t7-1][t9-1][t11-1] + s(a[t7], b[t9]) + s(a[t7], c[t11]) + s(b[t9], c[t11]), MAX(m1[t7][t9][t11], MAX(m2[t7][t9][t11], MAX(m3[t7][t9][t11], MAX(m4[t7][t9][t11], MAX(m5[t7][t9][t11], m6[t7][t9][t11])))))));;
                            }
                        }
                    }
                }
            }
        }
    }
/* End of CLooG code */


}