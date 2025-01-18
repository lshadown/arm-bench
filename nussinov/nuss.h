#include <iostream>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define floord(n,d) floor(((double)(n))/((double)(d)))
#define ceild(n,d) ceil(((double)(n))/((double)(d)))



using namespace std;

namespace Nussinov {
    char *RNA;
    int **S;

#include "library.h"
#include "oryg.h"
#include "tilecorr.h"
#include "dapt.h"
#include "pluto.h"




    void nussinov() {
        S = allocate2DArray(N);
        RNA = new char[N];

        if (Comp == "Pluto")
            plutoom();

        if (Comp == "Traco")
            tilecorr_omp();

        if (Comp == "Dapt")
            daptomp();

        if (Comp == "Original")
            oryg();

        deallocate2DArray(S);
        delete RNA;


    }
}