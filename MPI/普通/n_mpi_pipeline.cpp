#include<iostream>
#include <stdio.h>
#include<cstring>
#include<typeinfo>
#include <stdlib.h>
#include<cmath>
#include<mpi.h>
#include<windows.h>
#include<omp.h>
#include<tmmintrin.h>
#include<xmmintrin.h>
#include<emmintrin.h>
#include<pmmintrin.h>
#include<smmintrin.h>
#include<nmmintrin.h>
#include<immintrin.h>
#include <iostream>
#include <cstring>
#include <cstdlib>
#include <chrono>
#include <cstring>

long long head, tail, freq;
const int n = 1000;
float a[n][n];

void init() {
    memset(a, 0, sizeof(a));
    for (int i = 0; i < n; i++) {
        a[i][i] = 1;
        for (int j = i + 1; j < n; j++)
            a[i][j] = rand();
    }
    for (int k = 0; k < n; k++)
        for (int i = k + 1; i < n; i++)
            for (int j = 0; j < n; j++)
                a[i][j] += a[k][j];
}

double gauss(int argc, char* argv[]) {  //块划分
    double start_time = 0;
    double end_time = 0;
    MPI_Init(&argc, &argv);
    int total = 0;
    int rank = 0;
    int i = 0;
    int j = 0;
    int k = 0;
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &total);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int begin = n / total * rank;
    int end = (rank == total - 1) ? n : n / total * (rank + 1);
    if (rank == 0) {  //0号进程初始化矩阵
        init();

        for (j = 1; j < total; j++) {
            int b = j * (n / total), e = (j == total - 1) ? n : (j + 1) * (n / total);
            for (i = b; i < e; i++)
                MPI_Send(&a[i][0], n, MPI_FLOAT, j, 1, MPI_COMM_WORLD);//1是初始矩阵信息，向每个进程发送数据
        }

    }
    else {
        for (i = begin; i < end; i++) {
            MPI_Recv(&a[i][0], n, MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);
        }

    }

    MPI_Barrier(MPI_COMM_WORLD);  //此时每个进程都拿到了数据
    start_time = MPI_Wtime();
    for (k = 0; k < n; k++) {
        if ((begin <= k && k < end)) {
            for (j = k + 1; j < n; j++)
                a[k][j] = a[k][j] / a[k][k];
            a[k][k] = 1.0;
        }
        else if(k<begin) MPI_Recv(&a[k][0], n, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &status);
        for (i = max(begin, k + 1); i < end; i++) {
            for (j = k + 1; j < n; j++)
                a[i][j] = a[i][j] - a[i][k] * a[k][j];
            a[i][k] = 0;
        }
        if (rank != total - 1 && k < end) {
            MPI_Send(&a[k][0],n,MPI_FLOAT,rank+1,0,MPI_COMM_WORLD);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);	//各进程同步
    if (rank == 0) {
        end_time = MPI_Wtime();
        printf("%.4lf ms\n", 1000 * (end_time - start_time));
    }
    MPI_Finalize();
    return end_time - start_time;
}
int main(int argc, char* argv[]) {
    gauss(argc, argv);
    return 0;
}
