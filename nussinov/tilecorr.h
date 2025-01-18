//
// Created by marek on 15.11.2024.
//

#ifndef MY_APPLICATION_TILECORR_H
#define MY_APPLICATION_TILECORR_H

#endif //MY_APPLICATION_TILECORR_H





void tilecorr() {

    for (int c1 = 1; c1 < N + floord(N - 2, 128); c1 += 1) {
        // Vector to store threads
        std::vector<std::thread> threads;

        for (int c3 = MAX(0, -N + c1 + 1); c3 <= (c1 - 1) / 129; c3 += 1) {
            threads.emplace_back([&, c1, c3]() {
                for (int c4 = 0; c4 <= 1; c4 += 1) {
                    if (c4 == 1) {
                        for (int c9 = N - c1 + 129 * c3; c9 <= MIN(N - 1, N - c1 + 129 * c3 + 127); c9 += 1) {
                            for (int c10 = MAX(0, N - c1 + 129 * c3 - c9 + 1); c10 <= 1; c10 += 1) {
                                if (c10 == 1) {
                                    S[(N - c1 + c3 - 1)][c9] = MAX(S[(N - c1 + c3 - 1)][c9],
                                                                   S[(N - c1 + c3 - 1) + 1][c9 - 1] +
                                                                   can_pair(RNA, (N - c1 + c3 - 1), c9));
                                } else {
                                    for (int c11 = 128 * c3 + 1; c11 <= -N + c1 - c3 + c9; c11 += 1) {
                                        S[(N - c1 + c3 - 1)][c9] = MAX(S[(N - c1 + c3 - 1)][c11 + (N - c1 + c3 - 1)] +
                                                                       S[c11 + (N - c1 + c3 - 1) + 1][c9],
                                                                       S[(N - c1 + c3 - 1)][c9]);
                                    }
                                }
                            }
                        }
                    } else {
                        for (int c5 = 0; c5 <= 8 * c3; c5 += 1) {
                            for (int c9 = N - c1 + 129 * c3; c9 <= MIN(N - 1, N - c1 + 129 * c3 + 127); c9 += 1) {
                                for (int c11 = 16 * c5; c11 <= MIN(128 * c3, 16 * c5 + 15); c11 += 1) {
                                    S[(N - c1 + c3 - 1)][c9] = MAX(S[(N - c1 + c3 - 1)][c11 + (N - c1 + c3 - 1)] +
                                                                   S[c11 + (N - c1 + c3 - 1) + 1][c9],
                                                                   S[(N - c1 + c3 - 1)][c9]);
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

void tilecorr_omp()
{

    int c0,c1,c2,c3,c5,c6,c7,c9,c11,c10,c4,c12;

    for( c1 = 1; c1 < N + floord(N - 2, 128); c1 += 1)
#pragma omp parallel for schedule(dynamic, 1) shared(c1) private(c3,c4,c9,c10,c11,c5)
            for( c3 = MAX(0, -N + c1 + 1); c3 <= (c1 - 1) / 129; c3 += 1)
                for( c4 = 0; c4 <= 1; c4 += 1) {
                    if (c4 == 1) {
                        for( c9 = N - c1 + 129 * c3; c9 <= MIN(N - 1, N - c1 + 129 * c3 + 127); c9 += 1)
                            for( c10 = MAX(0, N - c1 + 129 * c3 - c9 + 1); c10 <= 1; c10 += 1) {
                                if (c10 == 1) {
                                    S[(N-c1+c3-1)][c9] = MAX(S[(N-c1+c3-1)][c9], S[(N-c1+c3-1)+1][c9-1] + can_pair(RNA, (N-c1+c3-1), c9));
                                } else
                                    for( c11 = 128 * c3 + 1; c11 <= -N + c1 - c3 + c9; c11 += 1)
                                        S[(N-c1+c3-1)][c9] = MAX(S[(N-c1+c3-1)][c11+(N-c1+c3-1)] + S[c11+(N-c1+c3-1)+1][c9], S[(N-c1+c3-1)][c9]);
                            }
                    } else
                        for( c5 = 0; c5 <= 8 * c3; c5 += 1)
                            for( c9 = N - c1 + 129 * c3; c9 <= MIN(N - 1, N - c1 + 129 * c3 + 127); c9 += 1)
                                for( c11 = 16 * c5; c11 <= MIN(128 * c3, 16 * c5 + 15); c11 += 1)
                                    S[(N-c1+c3-1)][c9] = MAX(S[(N-c1+c3-1)][c11+(N-c1+c3-1)] + S[c11+(N-c1+c3-1)+1][c9], S[(N-c1+c3-1)][c9]);
                }

}
