//
// Created by marek on 15.11.2024.
//

#ifndef MY_APPLICATION_DAPT_H
#define MY_APPLICATION_DAPT_H

#endif //MY_APPLICATION_DAPT_H



void daptomp(){


    int n = N;
    int c0,c1,c2,c3,c5,c6,c7,c9,c11,c10,c4,c12;
    for (int c0 = floord(-31 * n + 115, 3132) + 2; c0 <= floord(79 * n - 158, 2436) + 2; c0 += 1) {
#pragma omp parallel for
        for (int c1 = MAX(-c0 - (n + 52) / 54 + 2, -((n + 114) / 116)); c1 <= MIN(MIN(-c0 + (n - 2) / 42 + 1, c0 + ((-4 * c0 + 3)/31) - 1), (-21 * c0 + 20)/79); c1 += 1) {
            for (int c2 = MAX(-c0 + c1 + floord(21 * c0 - 17 * c1 - 21, 48) + 1, -c0 - c1 - (n - 42 * c0 - 42 * c1 + 136) / 96 + 1); c2 <= MIN(MIN(-1, -c0 - c1), -((27 * c0 - 31 * c1 + 54) / 69) + 1); c2 += 1) {
                for (int c5 = MAX(27 * c0 - 31 * c1 + 27 * c2 - 83, -42 * c2 - 41); c5 <= MIN(MIN(n + 54 * c0 + 54 * c1 + 54 * c2 - 1, -42 * c2), 54 * c0 - 62 * c1 + 54 * c2); c5 += 1) {
                    for (int c6 = MAX(-54 * c0 - 54 * c1 - 54 * c2, -116 * c1 - 2 * c5 - 114); c6 <= MIN(MIN(-54 * c0 - 54 * c1 - 54 * c2 + 53, n - c5 - 1), -116 * c1 - c5); c6 += 1) {
                        for (int c7 = MAX(-116 * c1 - 115, c5 + c6); c7 <= MIN(MIN(n - 1, -116 * c1), 2 * c5 + c6 - 1); c7 += 1) {
                            if (2 * c5 + c6 >= c7 + 2) {
                                S[c6][c7] = MAX(S[c6][-c5 + c7] + S[-c5 + c7 + 1][c7], S[c6][c7]);
                                if (c7 == c5 + c6) {
                                    S[c6][c5 + c6] = MAX(S[c6][c5 + c6], S[c6 + 1][c5 + c6 - 1] +  can_pair(RNA, c6, c5 + c6));
                                }
                            }
                            S[c6][c7] = MAX(S[c6][c5 + c6 - 1] + S[c5 + c6][c7], S[c6][c7]);
                            if (c7 == c5 + c6) {
                                S[c6][c5 + c6] = MAX(S[c6][c5 + c6], S[c6 + 1][c5 + c6 - 1] + can_pair(RNA, c6, c5 + c6));
                            }
                        }
                    }
                }
            }
        }
    }


}

void dapt() {
    int n = N;
    for (int c0 = floord(-31 * n + 115, 3132) + 2; c0 <= floord(79 * n - 158, 2436) + 2; c0 += 1) {
        std::vector<std::thread> threads;

        for (int c1 = MAX(-c0 - (n + 52) / 54 + 2, -((n + 114) / 116));
             c1 <= MIN(MIN(-c0 + (n - 2) / 42 + 1, c0 + ((-4 * c0 + 3) / 31) - 1), (-21 * c0 + 20) / 79);
             c1 += 1) {
                threads.emplace_back([&, c0, c1]() {
                for (int c2 = MAX(-c0 + c1 + floord(21 * c0 - 17 * c1 - 21, 48) + 1,
                                       -c0 - c1 - (n - 42 * c0 - 42 * c1 + 136) / 96 + 1);
                     c2 <= MIN(MIN(-1, -c0 - c1), -((27 * c0 - 31 * c1 + 54) / 69) + 1);
                     c2 += 1) {
                    for (int c5 = MAX(27 * c0 - 31 * c1 + 27 * c2 - 83, -42 * c2 - 41);
                         c5 <= MIN(MIN(n + 54 * c0 + 54 * c1 + 54 * c2 - 1, -42 * c2), 54 * c0 - 62 * c1 + 54 * c2);
                         c5 += 1) {
                        for (int c6 = MAX(-54 * c0 - 54 * c1 - 54 * c2, -116 * c1 - 2 * c5 - 114);
                             c6 <= MIN(MIN(-54 * c0 - 54 * c1 - 54 * c2 + 53, n - c5 - 1), -116 * c1 - c5);
                             c6 += 1) {
                            for (int c7 = MAX(-116 * c1 - 115, c5 + c6);
                                 c7 <= MIN(MIN(n - 1, -116 * c1), 2 * c5 + c6 - 1);
                                 c7 += 1) {
                                if (2 * c5 + c6 >= c7 + 2) {
                                    S[c6][c7] = MAX(S[c6][-c5 + c7] + S[-c5 + c7 + 1][c7], S[c6][c7]);
                                    if (c7 == c5 + c6) {
                                        S[c6][c5 + c6] = MAX(S[c6][c5 + c6],
                                                             S[c6 + 1][c5 + c6 - 1] + can_pair(RNA, c6, c5 + c6));
                                    }
                                }
                                S[c6][c7] = MAX(S[c6][c5 + c6 - 1] + S[c5 + c6][c7], S[c6][c7]);
                                if (c7 == c5 + c6) {
                                    S[c6][c5 + c6] = MAX(S[c6][c5 + c6],
                                                         S[c6 + 1][c5 + c6 - 1] + can_pair(RNA, c6, c5 + c6));
                                }
                            }
                        }
                    }
                }
            });
        }

        // Join all threads
        for (auto& th : threads) {
            if (th.joinable()) {
                th.join();
            }
        }
    }
}

void dapt_extra()
{

        int n = N;

        // Worker function for threads
        auto thread_worker = [&](int start_c1, int end_c1) {
            for (int c1 = start_c1; c1 <= end_c1; ++c1) {
                for (int c2 = std::max(-1, -c1 - (n + 52) / 54 + 2);
                     c2 <= std::min(-1, -c1);
                     ++c2) {
                    for (int c5 = std::max(27 * c2 - 83, -42 * c2 - 41);
                         c5 <= std::min(n + 54 * c2 - 1, -42 * c2);
                         ++c5) {
                        for (int c6 = std::max(-54 * c2, -116 * c1 - 2 * c5 - 114);
                             c6 <= std::min(-54 * c2 + 53, n - c5 - 1);
                             ++c6) {
                            for (int c7 = std::max(c5 + c6, -116 * c1 - 115);
                                 c7 <= std::min(n - 1, 2 * c5 + c6 - 1);
                                 ++c7) {
                                if (2 * c5 + c6 >= c7 + 2) {
                                    S[c6][c7] = MAX(S[c6][-c5 + c7] + S[-c5 + c7 + 1][c7], S[c6][c7]);
                                    if (c7 == c5 + c6) {
                                        S[c6][c5 + c6] = MAX(S[c6][c5 + c6], S[c6 + 1][c5 + c6 - 1] + can_pair(RNA, c6, c5 + c6));
                                    }
                                }
                                S[c6][c7] = MAX(S[c6][c5 + c6 - 1] + S[c5 + c6][c7], S[c6][c7]);
                                if (c7 == c5 + c6) {
                                    S[c6][c5 + c6] = MAX(S[c6][c5 + c6], S[c6 + 1][c5 + c6 - 1] + can_pair(RNA, c6, c5 + c6));
                                }
                            }
                        }
                    }
                }
            }
        };

        // Range for c1
        int start_c1 = -((n + 114) / 116);
        int end_c1 = (n - 2) / 42 + 1;

        int total_c1 = end_c1 - start_c1 + 1;
        int chunk_size = (total_c1 + num_threads - 1) / num_threads;

        // Spawn threads
        std::vector<std::thread> threads;
        for (int i = 0; i < num_threads; ++i) {
            int thread_start_c1 = start_c1 + i * chunk_size;
            int thread_end_c1 = std::min(start_c1 + (i + 1) * chunk_size - 1, end_c1);

            if (thread_start_c1 <= thread_end_c1) {
                threads.emplace_back(thread_worker, thread_start_c1, thread_end_c1);
            }
        }

        // Join threads
        for (auto& th : threads) {
            if (th.joinable()) {
                th.join();
            }
        }

}