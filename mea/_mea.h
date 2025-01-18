int _mea(){

    int i,j,k,ll,p,q;
    int c0, c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12,c13,c15;

    int t1, t2, t3, t4, t5, t6,t7;
    int lb, ub, lbp, ubp, lb2, ub2;
    int lbv, ubv;



    RNA =  (char*) malloc(DIM * sizeof(char*));  //read from FASTA file
    Puu = (float *) malloc(DIM * sizeof(float*));



    //printf("Sequence: ");
    //for(i=0; i<N; i++)
    //   printf("%c", RNA[i]);
    //printf("\n\n");




    rna_array_init(M, 0, 0);
    int a = 0;

    double start = omp_get_wtime();
    //  compute the partition functions Q and Qbp
    if(kind==1){
#pragma scop
        for(i=N-1; i>=0; i--){
            for(j=i+1; j<N; j++){
                for(k=0; k<j-i-l; k++){
                    //      a++;
                    M[i][j] = MAX(M[i][j], M[i][k+i-1] + M[k+i+1][j-1] + delta*Pbp[k+i][j])*paired(k+i,j-1);
                }
                M[i][j] = MAX(M[i][j], M[i][j-1] + Puu[j-1]);
            }
        }
#pragma endscop
    }
    if(kind==2) // pluto
    {
        printf("pluto\n");
/*      if (N >= 2) {
  for (t1=1;t1<=N-1;t1++) {
    lbp=0;
    ubp=t1-1;
#pragma omp parallel for private(lbv,ubv,t3,t4,t5)
    for (t2=lbp;t2<=ubp;t2++) {
      for (t3=0;t3<=floord(t1-t2-1,16);t3++) {
        for (t5=16*t3;t5<=min(16*t3+15,t1-t2-1);t5++) {
          M[t2][t1] = MAX(M[t2][t1], M[t2][t5+t2-1] + M[t5+t2+1][t1-1] + delta*Pbp[t5+t2][t1])*paired(t5+t2,t1-1);;
        }
      }
      t3 = floord(t1,16);
      M[t2][t1] = MAX(M[t2][t1], M[t2][t1-1] + Puu[t1-1]);;
    }
  }
}
*/
/* We do not support C11 <threads.h>.  */
        int t1, t2, t3, t4, t5;
        int lb, ub, lbp, ubp, lb2, ub2;
        int lbv, ubv;
        /* Start of CLooG code */
        if (N >= 2) {
            for (t1=1;t1<=N-1;t1++) {
                lbp=0;
                ubp=floord(t1-1,16);
#pragma omp parallel for private(lbv,ubv,t3,t4,t5)
                for (t2=lbp;t2<=ubp;t2++) {
                    for (t3=0;t3<=floord(t1,16);t3++) {
                        if ((t1 <= 16*t3+15) && (t2 == 0)) {
                            for (t4=0;t4<=t1-16*t3-1;t4++) {
                                for (t5=16*t3;t5<=t1-t4-1;t5++) {
                                    M[t4][t1] = MAX(M[t4][t1], M[t4][t5+t4-1] + M[t5+t4+1][t1-1] + delta*Pbp[t5+t4][t1])*paired(t5+t4,t1-1);;
                                }
                                M[t4][t1] = MAX(M[t4][t1], M[t4][t1-1] + Puu[t1-1]);;
                            }
                        }
                        if (t1 >= 16*t3+16) {
                            for (t4=16*t2;t4<=min(16*t2+15,t1-16*t3-1);t4++) {
                                for (t5=16*t3;t5<=min(16*t3+15,t1-t4-1);t5++) {
                                    M[t4][t1] = MAX(M[t4][t1], M[t4][t5+t4-1] + M[t5+t4+1][t1-1] + delta*Pbp[t5+t4][t1])*paired(t5+t4,t1-1);;
                                }
                            }
                        }
                        if (t1 <= 16*t3+15) {
                            for (t4=max(16*t2,t1-16*t3);t4<=min(t1-1,16*t2+15);t4++) {
                                M[t4][t1] = MAX(M[t4][t1], M[t4][t1-1] + Puu[t1-1]);;
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
/*
for( c1 = 2; c1 < 2 * N + floord(N - 2, 16) - 1; c1 += 1)
  #pragma omp parallel for
  for( c3 = max(-2 * N + c1 + 2, -((c1 - 1) % 2) + 1); c3 <= (c1 - 2) / 33; c3 += 2) {
    for( c5 = 0; c5 <= c3; c5 += 1)
      for( c9 = ((-c1 + 33 * c3) / 2) + N; c9 <= min(N - 1, ((-c1 + 33 * c3) / 2) + N + 15); c9 += 1)
        for( c11 = 16 * c5; c11 <= min(min(16 * c3 + 1, 16 * c5 + 15), ((c1 - c3) / 2) - N + c9); c11 += 1)
          M[(((-c1+c3)/2)+N-1)][c9] = MAX(M[(((-c1+c3)/2)+N-1)][c9], M[(((-c1+c3)/2)+N-1)][c11+(((-c1+c3)/2)+N-1)-1] + M[c11+(((-c1+c3)/2)+N-1)+1][c9-1] + delta*Pbp[c11+(((-c1+c3)/2)+N-1)][c9])*paired(c11+(((-c1+c3)/2)+N-1),c9-1);
    for( c9 = ((-c1 + 33 * c3) / 2) + N; c9 <= min(N - 1, ((-c1 + 33 * c3) / 2) + N + 15); c9 += 1)
      for( c10 = max(0, 8 * c3 - c9 + (2 * N - c1 + c3 + 2 * c9 - 1) / 4 + 2); c10 <= 1; c10 += 1) {
        if (c10 == 1) {
          M[(((-c1+c3)/2)+N-1)][c9] = MAX(M[(((-c1+c3)/2)+N-1)][c9], M[(((-c1+c3)/2)+N-1)][c9-1] + Puu[c9-1]);
        } else {
          for( c11 = 16 * c3 + 2; c11 <= ((c1 - c3) / 2) - N + c9; c11 += 1)
            M[(((-c1+c3)/2)+N-1)][c9] = MAX(M[(((-c1+c3)/2)+N-1)][c9], M[(((-c1+c3)/2)+N-1)][c11+(((-c1+c3)/2)+N-1)-1] + M[c11+(((-c1+c3)/2)+N-1)+1][c9-1] + delta*Pbp[c11+(((-c1+c3)/2)+N-1)][c9])*paired(c11+(((-c1+c3)/2)+N-1),c9-1);
        }
      }
  }
  */
        for( c0 = 1; c0 < N + floord(N - 2, 16); c0 += 1)
#pragma omp parallel for
                for( c1 = c0 - (c0 + 16) / 17 + 1; c1 <= min(N - 1, c0); c1 += 1)
                    for( c3 = 16 * c0 - 16 * c1 + 1; c3 <= min(c1, 16 * c0 - 16 * c1 + 16); c3 += 1) {
                        for( c4 = 0; c4 <= c0 - c1; c4 += 1)
                            for( c10 = 16 * c4; c10 <= min(c3 - 1, 16 * c4 + 15); c10 += 1){
                                //    a++;
                                M[(N-c1-1)][(N-c1+c3-1)] = MAX(M[(N-c1-1)][(N-c1+c3-1)], M[(N-c1-1)][c10+(N-c1-1)-1] + M[c10+(N-c1-1)+1][(N-c1+c3-1)-1] + delta*Pbp[c10+(N-c1-1)][(N-c1+c3-1)])*paired(c10+(N-c1-1),(N-c1+c3-1)-1);
                            } M[(N-c1-1)][(N-c1+c3-1)] = MAX(M[(N-c1-1)][(N-c1+c3-1)], M[(N-c1-1)][(N-c1+c3-1)-1] + Puu[(N-c1+c3-1)-1]);
                    }




    }
    if(kind==4) // traco tstile
    {
        printf("traco corr\n");
        for( c1 = 1; c1 < N + floord(N - 2, 128); c1 += 1)
#pragma omp parallel for schedule(dynamic, 1)
                for( c3 = max(0, -N + c1 + 1); c3 <= (c1 - 1) / 129; c3 += 1)
                    for( c4 = 0; c4 <= 1; c4 += 1) {
                        if (c4 == 1) {
                            for( c9 = N - c1 + 129 * c3; c9 <= min(N - 1, N - c1 + 129 * c3 + 127); c9 += 1)
                                for( c10 = max(0, -c1 + 64 * c3 - c9 + (N + c1 + c3 + c9 + 1) / 2 + 1); c10 <= 1; c10 += 1) {
                                    if (c10 == 1) {
                                        M[(N-c1+c3-1)][c9] = MAX(M[(N-c1+c3-1)][c9], M[(N-c1+c3-1)][c9-1] + Puu[c9-1]);
                                    } else {
                                        for( c11 = 128 * c3 + 2; c11 <= -N + c1 - c3 + c9; c11 += 1)
                                            M[(N-c1+c3-1)][c9] = MAX(M[(N-c1+c3-1)][c9], M[(N-c1+c3-1)][c11+(N-c1+c3-1)-1] + M[c11+(N-c1+c3-1)+1][c9-1] + delta*Pbp[c11+(N-c1+c3-1)][c9])*paired(c11+(N-c1+c3-1),c9-1);
                                    }
                                }
                        } else {
                            for( c5 = 0; c5 <= 8 * c3; c5 += 1)
                                for( c9 = N - c1 + 129 * c3; c9 <= min(N - 1, N - c1 + 129 * c3 + 127); c9 += 1)
                                    for( c11 = 16 * c5; c11 <= min(min(128 * c3 + 1, 16 * c5 + 15), -N + c1 - c3 + c9); c11 += 1)
                                        M[(N-c1+c3-1)][c9] = MAX(M[(N-c1+c3-1)][c9], M[(N-c1+c3-1)][c11+(N-c1+c3-1)-1] + M[c11+(N-c1+c3-1)+1][c9-1] + delta*Pbp[c11+(N-c1+c3-1)][c9])*paired(c11+(N-c1+c3-1),c9-1);
                        }
                    }



    }

    if(kind==5){  // dapt
        for (int w0 = max(floord(-l + 1, 16) - 1, floord(-N - l + 4, 16) - 1); w0 <= floord(2 * N - l - 2, 16); w0 += 1) {
#pragma omp parallel for
            for (int h0 = max(-((N + 13) / 16), -((2 * N - l - 16 * w0 + 45) / 32) + 1); h0 <= min(0, w0  + 1 +  (l - 2)/16); h0 += 1) {
                for (int h1 = max(max(max(w0 - h0 + floord(-N + l, 16) + 1, w0 - 2 * h0 + floord(-N + l + 1, 16)), -h0 + floord(l + 16 * w0 + 1, 32)), -h0 + floord(l + 16 * w0 + 16 * h0 + 16, 32)); h1 <= min(min(floord(N - 1, 16), w0 - h0 + floord(l - 2, 16) + 1), -h0 + floord(l + 16 * w0 + 15, 32)); h1 += 1) {
                    for (int i0 = max(max(max(max(-N + 2, 16 * h0), l + 16 * w0 - 16 * h0 - 32 * h1 - 15), -16 * h1 - 14), -N + l + 16 * w0 - 16 * h0 - 16 * h1 + 1); i0 <= min(min(0, 16 * h0 + 15), l + 16 * w0 - 16 * h0 - 32 * h1 + 15); i0 += 1) {
                        for (int i1 = max(max(16 * h1, l + 16 * w0 - 16 * h0 - 16 * h1 - i0), -i0 + 1); i1 <= min(min(N - 1, 16 * h1 + 15), l + 16 * w0 - 16 * h0 - 16 * h1 - i0 + 15); i1 += 1) {
                            M[-i0][i1] = MAX(M[-i0][i1], M[-i0][i1 - 1] + Puu[i1 - 1]);
                        }
                    }
                }
            }
        }
    }

    double stop = omp_get_wtime();
    printf("%.4f   %i\n",stop - start, a);

    //printf("Q\n");
    //rna_array_print(Q);
    //printf("Qbp\n");
    //rna_array_print(Qbp);



    return 0;

}