#include <iostream>
#include <fstream>
#include <string>
#include <thread>
#include <chrono>
#include <vector>
#include <algorithm>
#include <omp.h>

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

void runAndLog(const std::string &algorithmName, const std::function<void()>& algorithm, const std::vector<int>& testSizes) {
    std::ofstream outFile("./result/" + algorithmName + ".txt");

    if (!outFile) {
        std::cerr << "Cannot open file: " << algorithmName + ".txt" << std::endl;
        return;
    }

    std::vector<std::string> implementations = {"Pluto", "Traco", "Dapt", "Original"};

    for (const auto& comp : implementations) {
        Comp = comp;
        outFile << "Implementation: " << comp << "\n";
        std::cout << "Running " << algorithmName << " with Comp = " << comp << std::endl;

        for (int a : testSizes) {
            try {
                N = a;
                double start = omp_get_wtime();
                algorithm();
                double stop = omp_get_wtime();

                outFile << "N: " << a << ", Elapsed time: " << (stop - start) << " seconds\n";
                std::cout << algorithmName << " with N = " << a << ", Comp = " << comp << " completed in "
                          << (stop - start) << " seconds\n";
            } catch (const std::exception& e) {
                outFile << "N: " << a << ", ERROR: " << e.what() << "\n";
                std::cerr << "Error in " << algorithmName << " with N = " << a << ", Comp = " << comp << ": " << e.what() << "\n";
            } catch (...) {
                outFile << "N: " << a << ", ERROR: Unknown exception occurred\n";
                std::cerr << "Unknown error in " << algorithmName << " with N = " << a << ", Comp = " << comp << "\n";
            }
        }
        outFile << "\n";
    }

    outFile.close();
}



int main() {
    std::vector<int> testSizes = {1000, 2500, 5000, 7500, 10000, 15000, 20000, 30000};

    //runAndLog("Nussinov", []() { Nussinov::nussinov(); }, { 2200, 5000, 10000});
    //runAndLog("Counting", []() { Counting::counting(); }, { 2200, 5000, 10000});
    //runAndLog("Knuth", []() { Knuth::knuth(); }, { 2200, 5000});
    //runAndLog("Mcc", []() { Mcc::mcc(); }, { 1000, 2200, 5000});
    //runAndLog("Triang", []() { Triang::triang(); }, { 10000});
    //runAndLog("Zuker", []() { Zuker::zuker(); }, { 5000});
    //runAndLog("NW", []() { NW::nw(); }, { 2200, 5000, 10000});
    //runAndLog("SW", []() { SW::sw(); }, {  10000});
    //runAndLog("SW3D", []() { SW3D::sw3d(); }, { 2200});
    runAndLog("MEA", []() { MEA::mea(); }, { 2200, 5000, 10000});
    
    return 0;
}

