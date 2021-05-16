#ifndef LOCAL_DOUBLE_H
#define LOCAL_DOUBLE_H

#include <vector>
#include <iostream>
#include <mpi.h>
#include <string>
#include <fstream>
#include <cstdio>
#include <ctime>
#include "magickNumbers.h"
#include "perpendicularDoubleStepSteps.h"
#include "measure.h"

template<int N, int lp, int D>
bool __localDoubleStep(int rank, complex<double> **imageMatrix, vector<complex<double>> &W, vector<vector<perd::line>> &required_lines_for_perpendicular ,double &t);


template<int N>
void local_double_step(int rank){
    double t;
    const int D = mn::N/mn::size;
    //create image matrix
    complex<double> **imageMatrix = new complex<double>*[D];
    for( int i = 0; i < D; ++i){
        imageMatrix[i] = new complex<double>[N];
        for( int j = 0; j < N; ++j){
            imageMatrix[i][j] = complex<double>(0,0);
        }
    }

    //calculating lines for perpendicular doublestep
    vector<vector<perd::line>> required_lines_for_perpendicular(log2(N)+1);
    vector<perd::line> lines;
    for(int i = D*rank; i < D*rank+D; ++i){
        lines.push_back(perd::line(i,0));
    }
    perd::calc_required_lines<N>(lines, 0, N, required_lines_for_perpendicular);


    //calculating W
    vector<complex<double>> W (N);
    complex<double> W_0 = exp(complex<double>(0,2*M_PI/N));
    for( int i = 0; i < N; ++i){
        W[i] = pow(W_0,i);
    }

    bool ceep_running = true;
    while(ceep_running){
        ceep_running = __localDoubleStep<mn::N, mn::lp, D>(rank, imageMatrix,W, required_lines_for_perpendicular, t);
    }

    for(int i = 0; i < D; ++i){
        delete imageMatrix[i];
    }
    delete imageMatrix;
}

template<int N, int lp, int D>
bool __localDoubleStep(int rank, complex<double> **imageMatrix, vector<complex<double>> &W, vector<vector<perd::line>> &required_lines_for_perpendicular ,double &t){

    //int buf[ mn::visibilities_per_message*6+1 ]; //note here v, u, vis_i, vis_r, isCS, shift_index
    int *buf = new int[mn::visibilities_per_message*6+1];


    int res = MPI_Recv( buf, mn::visibilities_per_message*6+1, MPI_INT, mn::MPI_host, mn::performDoubleStep, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (buf[mn::visibilities_per_message*6] == mn::COMMANDS::stop){
        return false;
    } else if (buf[mn::visibilities_per_message*6] == mn::COMMANDS::safe_results){
        __safe_local_result<N,D>(rank, imageMatrix);
        return true;
    } else if(buf[mn::visibilities_per_message*6] == mn::COMMANDS::safe_timer){
        __safe_local_timer(rank, t, "time per run" );
        return true;
    } else if(buf[mn::visibilities_per_message*6] == mn::COMMANDS::start_timer){
        t = omp_get_wtime();
        return true;
    }


    vector<vector<visibility<N>>> *parallelVisibilities= new vector<vector<visibility<N>>>(N);
    vector<vector<visibility<N>>> *perpendicularVisibilities = new vector<vector<visibility<N>>>(N);



    for( int i = 0; i < mn::visibilities_per_message*6; i += 6){
        double vis_r = *(float*)&buf[i+mn::VIS_ENCODING::vis_r];
        double vis_i = *(float*)&buf[i+mn::VIS_ENCODING::vis_i];
        int u = buf[i+mn::VIS_ENCODING::u];
        int v = buf[i+mn::VIS_ENCODING::v];
        int shift_index = buf[i+mn::VIS_ENCODING::shift_index];
        bool isCS = buf[i+mn::VIS_ENCODING::isCS];
        visibility<N> vis(complex<double>(vis_r, vis_i), u, v, isCS, shift_index);
        if(isCS){
            (*perpendicularVisibilities)[shift_index].push_back(vis);
        } else {
            (*parallelVisibilities)[shift_index].push_back(vis);
        }
    }

    //perform SPIFT
    pard::parallelDoubleStep<N,lp,D>(*parallelVisibilities,W,imageMatrix,rank);
    perd::perpendicularDoubleStep<N,lp,D>(*perpendicularVisibilities,W,imageMatrix,rank, required_lines_for_perpendicular);

    delete buf;

    return true;
}


#endif
