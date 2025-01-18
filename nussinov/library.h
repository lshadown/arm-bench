//
// Created by marek on 15.11.2024.
//

#ifndef MY_APPLICATION_LIBRARY_H
#define MY_APPLICATION_LIBRARY_H

#endif //MY_APPLICATION_LIBRARY_H

int paired(char a1, char a2) {
    if(a1 == 'A' && a2 == 'U') return 1;
    if(a1 == 'U' && a2 == 'A') return 1;
    if(a1 == 'G' && a2 == 'C') return 1;
    if(a1 == 'C' && a2 == 'G') return 1;
    return 0;
}

int can_pair(char *RNA, int i, int j)
{
    return paired(RNA[i], RNA[j]);
}



int** allocate2DArray(int n) {
    // Alokujemy ciągły blok pamięci dla n*n elementów typu int
    int* data = new int[n * n];

    // Alokujemy tablicę wskaźników do wierszy
    int** array = new int*[n];

    // Inicjalizujemy wskaźniki tak, aby każdy wskazywał na odpowiedni segment w ciągłym bloku pamięci
    for (int i = 0; i < n; ++i) {
        array[i] = data + i * n;
    }

    return array;
}

void deallocate2DArray(int** array) {
    if (array != nullptr) {
        delete[] array[0];  // Zwalniamy ciągły blok pamięci
        delete[] array;     // Zwalniamy tablicę wskaźników
    }
}
