//
// Created by marek on 15.11.2024.
//
#include <iostream>
#include <thread>
#include <vector>
#include <algorithm>

#ifndef MY_APPLICATION_PLUTO_H
#define MY_APPLICATION_PLUTO_H

#endif //MY_APPLICATION_PLUTO_H


void plutoom() {

    int t1, t2, t3, t4, t5, t6, t7, t8, t9, t10;
    int lb, ub, lbp, ubp, lb2, ub2;
    int lbv, ubv;
/* Start of CLooG code */
    if (N >= 2) {
        for (t2 = MAX(-1, ceild(-N - 13, 16)); t2 <= floord(N - 1, 16); t2++) {
            lbp = MAX(0, t2);
            ubp = MIN(floord(N - 1, 16), floord(16 * t2 + N + 13, 16));
#pragma omp parallel for private(lbv, ubv, t5, t6, t7, t8, t9, t10)
            for (t4 = lbp; t4 <= ubp; t4++) {
                for (t5 = MAX(MAX(-N + 2, 16 * t2 - 16 * t4), -16 * t4 - 14);
                     t5 <= MIN(0, 16 * t2 - 16 * t4 + 15); t5++) {
                    for (t7 = MAX(16 * t4, -t5 + 1); t7 <= MIN(N - 1, 16 * t4 + 15); t7++) {
                        for (t9 = 0; t9 <= t5 + t7 - 1; t9++) {
                            S[-t5][t7] = MAX(S[-t5][t9 + -t5] + S[t9 + -t5 + 1][t7], S[-t5][t7]);;
                        }
                        S[-t5][t7] = MAX(S[-t5][t7], S[-t5 + 1][t7 - 1] + can_pair(RNA, -t5, t7));;
                    }
                }
            }
        }
    }
}
/* End of CLooG code */



void pluto() {

    if (N >= 2) {
        for (int t2 = MAX(-1, ceild(-N - 13, 16)); t2 <= floord(N - 1, 16); t2++) {
            int lbp = MAX(0, t2);
            int ubp = MIN(floord(N - 1, 16), floord(16 * t2 + N + 13, 16));

            // Vector to store threads
            std::vector<std::thread> threads;

            for (int t4 = lbp; t4 <= ubp; t4++) {
                threads.emplace_back([&, t2, t4]() {   //lambda
                    for (int t5 = MAX(MAX(-N + 2, 16 * t2 - 16 * t4), -16 * t4 - 14);
                         t5 <= MIN(0, 16 * t2 - 16 * t4 + 15); t5++) {
                        for (int t7 = MAX(16 * t4, -t5 + 1);
                             t7 <= MIN(N - 1, 16 * t4 + 15); t7++) {
                            for (int t9 = 0; t9 <= t5 + t7 - 1; t9++) {
                                {
                                    S[-t5][t7] = MAX(S[-t5][t9 + -t5] + S[t9 + -t5 + 1][t7],
                                                     S[-t5][t7]);
                                }
                            }
                            {
                                S[-t5][t7] = MAX(S[-t5][t7],
                                                 S[-t5 + 1][t7 - 1] + can_pair(RNA, -t5, t7));
                            }
                        }
                    }
                });
            }

            // Join all threads
            for (auto &th: threads) {
                if (th.joinable()) {
                    th.join();
                }
            }
        }
    }
}