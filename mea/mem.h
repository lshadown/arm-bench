#include <stdio.h>
#include <stdlib.h>
#include <time.h>


int **mem()
{
    int i;
    int **S;
    S = (int **) malloc(
            DIM * sizeof(int*));

    for (i=0; i<DIM; i++)
        S[i] = (int*)malloc(DIM * sizeof(int));

    return S;
}

float **memd()
{
    int i;
    float **S = (float **) malloc(DIM * sizeof(float*));

    for (i=0; i<DIM; i++)
        S[i] = (float*)malloc(DIM * sizeof(float));




    return S;
}




void rna_array_init(float **S, float def, float def2){

    int i,j;

    for(i=0; i<=N; i++)
        for(j=0; j<=N; j++)
            if(i==j || i==0)
                S[i][j] = def;
            else
                S[i][j] = def2;



}

void rna_array_print(float **S){

    int i,j;

    for(i=0; i<N; i++){
        for(j=0; j<N; j++)
            if(i>j)
                printf("       ");
            else
                printf(" %5.3f ", S[i][j]);
        printf("\n");
    }
    printf("\n");
}