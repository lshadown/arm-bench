#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <omp.h>
#include <math.h>


namespace MEA{





    #define min(a,b) (((a)<(b))?(a):(b))
    #define MIN(a,b) (((a)<(b))?(a):(b))
    #define max(a,b) (((a)>(b))?(a):(b))
    #define MAX(a,b) (((a)>(b))?(a):(b))
    #define floord(n,d) floor(((double)(n))/((double)(d)))
    #define ceild(n,d) ceil(((double)(n))/((double)(d)))

    float ** Q;
    float ** Q1;
    float ** Qbp;
    float ** Pbp;
    float ** Pu;
    float *Puu;
    float ** M;


    int Ebp = 1; // Energy weight of base pair  -2, -1, 0, 1, 2
    int RT = 1; // 'Normalized' temperature 1,2,3,4,5
    float ERT;
    int l = 0; //minimum loop length 0-5
    int delta = 1;  // Base pair weighting  1-5

    char * RNA;  //only ACGU

    int DIM;
    int kind;
    int num_proc = num_threads;


    int paired(int i, int j) {
        char nt1 = RNA[i];
        char nt2 = RNA[j];
        if ((nt1 == 'A' && nt2 == 'U') || (nt1 == 'U' && nt2 == 'A') ||
            (nt1 == 'G' && nt2 == 'C') || (nt1 == 'C' && nt2 == 'G') ||
            (nt1 == 'G' && nt2 == 'U') || (nt1 == 'U' && nt2 == 'G')){

            return 1;}
        else
            return 0;
    }

    void rand_seq(unsigned char*a, int N){
        int i, tmp;
        srand(time(NULL));
        for(i=0; i<N; i++)
        {
            tmp = rand()%4;

            switch(tmp){
                case 0 : a[i] = 'A'; break;
                case 1 : a[i] = 'G'; break;
                case 2 : a[i] = 'C'; break;
                case 3 : a[i] = 'U'; break;
            }

        }

    }


    #include "mem.h"
    #include "mcc.h"
    #include "pb.h"
    #include "pu.h"
    #include "_mea.h"




    void mea() {
        DIM = N+2;

        if(Comp == "Original")
            kind = 1;

        if(Comp == "Dapt")
            kind = 5;

        if(Comp == "Pluto")
            kind = 2;

        if(Comp == "Traco")
            kind = 3;

       mcc();
       pb();
       pu();
       _mea();

    }

}