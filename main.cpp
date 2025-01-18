#include <string>
#include <thread>
#include <chrono>
#include <vector>
#include <algorithm>
#include <omp.h>




// github check
using namespace std;

int num_threads = 10;
int N = 10000;
string Comp;



#include "nussinov/nuss.h"
#include "counting/counting.h"
#include "knuth/knuth.h"
#include "mcc/mcc.h"
#include "triang/triang.h"
#include "nw/nw.h"
#include "zuker/zuker.h"
#include "sw/sw.h"
#include "sw3d/sw3d.h"
#include "mea/mea.h"


int main() {
    double start = omp_get_wtime();

    //Nussinov::nussinov();
    //Counting::counting();
    //Knuth::knuth();
    //Mcc::mcc();
    //Triang::triang();
    //Zuker::zuker();
    //NW::nw();
    //SW::sw();
    //SW3D::sw3d();
    MEA::mea();


    double stop = omp_get_wtime();

    cout <<  "Elapsed time: " + to_string(stop-start);

    return 0;
}
