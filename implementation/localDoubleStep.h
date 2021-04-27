#ifndef LOCAL_DOUBLE_H
#define LOCAL_DOUBLE_H

#include <vector>
#include <iostream>
#include <mpi.h>
#include "magickNumbers.h"
#include "perpendicularDoubleStepSteps.h"

template<int N, int lp, int D>
void __localDoubleStep(int rank, complex<double> **imageMatrix, vector<complex<double>> &W);

template<int N>
void local_double_step(int rank){

    std::cout << "local process " << endl;

    const int D = mn::N/mn::size;
    //create image matrix
    complex<double> **imageMatrix = new complex<double>*[D];
    for( int i = 0; i < D; ++i){
        imageMatrix[i] = new complex<double>[N];
        for( int j = 0; j < N; ++j){
            imageMatrix[i][j] = complex<double>(0,0);
        }
    }


    //calculating W
    vector<complex<double>> W (N);
    complex<double> W_0 = exp(complex<double>(0,2*M_PI/N));
    for( int i = 0; i < N; ++i){
        W[i] = pow(W_0,i);
    }

    while(true){
        __localDoubleStep<mn::N, mn::lp, D>(rank, imageMatrix,W);
    }

    for(int i = 0; i < D; ++i){
        delete imageMatrix[i];
    }
    delete imageMatrix;
}

template<int N, int lp, int D>
void __localDoubleStep(int rank, complex<double> **imageMatrix, vector<complex<double>> &W){

    int buf[ mn::visibilities_per_message*6 ]; //note here v, u, vis_i, vis_r, isCS, shift_index
    int res = MPI_Recv( &buf, mn::visibilities_per_message*6, MPI_INT, mn::MPI_host, mn::performDoubleStep, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    vector<vector<visibility<N>>> parallelVisibilities(N);
    vector<vector<visibility<N>>> perpendicularVisibilities(N);

    for( int i = 0; i < N; i += 6){
        double vis_r = *(float*)&buf[i+mn::VIS_ENCODING::vis_r];
        double vis_i = *(float*)&buf[i+mn::VIS_ENCODING::vis_i];
        int u = buf[i+mn::VIS_ENCODING::u];
        int v = buf[i+mn::VIS_ENCODING::v];
        int shift_index = buf[i+mn::VIS_ENCODING::shift_index];
        bool isCS = buf[i+mn::VIS_ENCODING::isCS];
        visibility<N> vis(complex<double>(vis_r, vis_i), u, v, isCS, shift_index);
        if(isCS){
            perpendicularVisibilities[shift_index].push_back(vis);
        } else {
            parallelVisibilities[shift_index].push_back(vis);
        }
    }

    //print out visibilities
    for( int j = 0; j < N; ++j){
    for(int i = 0; i < perpendicularVisibilities[j].size(); ++i){
        visibility<N> v = perpendicularVisibilities[j][i];
        std::cout << "per [" << v.vis << ", " << v.u <<  ", "<< v.v << "]";
    }
    for(int i = 0; i < parallelVisibilities[j].size(); ++i){
        visibility<N> v = parallelVisibilities[j][i];
        std::cout << "par [" << v.vis << ", " << v.u <<  ", "<< v.v << "]" << "shift-index: " << v.shift_index << endl;
    }
    }

    //perform SPIFT
    perd::perpendicularDoubleStep<N,lp,D>(perpendicularVisibilities,W,imageMatrix,rank);
    pard::parallelDoubleStep<N,lp,D>(parallelVisibilities,W,imageMatrix,rank);
}

#endif
