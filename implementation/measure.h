#ifndef MEASURE_H
#define MEASURE_H

#include <mpi.h>
#include <string>
#include <fstream>
#include <cstdio>
#include <ctime>

using namespace std;

template<int N, int D>
void __safe_local_result(int rank, complex<double> **imageMatrix){
    std::ofstream ofs;
    ofs.open("results/doublestep_test_result_rank_" + std::to_string(rank) + ".txt");
    for(int row; row < D; ++row){
        for(int column = 0; column < N; ++column){
            ofs << imageMatrix[row][column] << ", ";
        }
        ofs << endl;
    }
    ofs.close();
}

void __safe_local_timer(int rank, double &start, string postfix){
    double dif = (omp_get_wtime() - start );
    std::ofstream ofs;
    ofs.open("results/doublestep_test_times_rank_" + std::to_string(rank)+"_" + postfix + ".txt", std::ios_base::app);
    ofs << dif << endl;
    ofs.close();
}

#endif
