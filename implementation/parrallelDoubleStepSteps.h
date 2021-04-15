#ifndef PARA_H
#define PARA_H

#include<vector>
#include<complex>
#include"visibilities.h"
#include<math.h>
#include<iostream>
#include<omp.h>

using namespace std;

namespace pard{

template<int N, int lp, int D>
void __createShiftVectors(complex<double> **mem1, vector<vector<visibility<N>>> parallelVisibilities, vector<complex<double>> W);

template<int N, int lp, int D>
void __cacDoubleStepPara(complex<double> **mem1, complex<double> **mem2);

template<int N, int lp, int D>
void __calcFinalImageMatrix(complex<double> **mem1, complex<double> **mem2, complex<double> **imageMatrix, int machine_index);

template<int N, int lp, int D>
void parallelDoubleStep( vector<vector<visibility<N>>> parallelVisibilities, vector<complex<double>> W, complex<double> **imageMatrix, int machine_index ){
    cout << " is in parallel " << endl;

    //creating memory
    complex<double> **mem1 = new complex<double>*[N];
    complex<double> **mem2 = new complex<double>*[N];
    for(int i = 0; i < N; ++i){
            mem1[i] = new complex<double>[N];
            mem2[i] = new complex<double>[N];
    }

    __createShiftVectors<N,lp,D>(mem1,parallelVisibilities,W);

    cout << " created shift-vectors " << endl;

    __cacDoubleStepPara<N,lp,D>(mem1, mem2);

    cout << " performed steps of doublestep " << endl;

    __calcFinalImageMatrix<N,lp,D>(mem1,mem2,imageMatrix,machine_index);

    cout << " done with parallel " << endl;

}


template<int N, int lp, int D>
void __createShiftVectors(complex<double> **mem1,vector<vector<visibility<N>>> parallelVisibilities, vector<complex<double>> W){
    int r = N/lp;

    //calculating shift-vectors and storing it in memory 1
    #pragma omp parallel for
    for( int k = 0; k < lp; ++k){

        for( int shift_index = k*r; shift_index < (k+1)*r; ++shift_index){

            //set up with initial values (0 if no visibilities provided
            if(parallelVisibilities[shift_index].size() == 0){
                for(int i = 0; i < N; ++i){
                    mem1[shift_index][i] = complex<double>(0,0);
                }
            } else {
                visibility<N> v = parallelVisibilities[shift_index][0];
                for( int i = 0; i < N; ++i){
                    mem1[shift_index][i] = v.vis*W[(i*v.v)%N];
                }
                for(unsigned int j = 1; j < parallelVisibilities[shift_index].size(); ++ j){
                    visibility<N> vis = parallelVisibilities[shift_index][j];
                    for( int i = 0; i < N; ++i){
                        mem1[shift_index][i] += vis.vis*W[(i*v.v)%N];
                    }
                }
            }
        }
    }

}

template<int N, int lp, int D>
void __cacDoubleStepPara(complex<double> **mem1, complex<double> **mem2){
    //calculating doublestep
    //for small j
    for( int j = 1; j <= log2(D) && j <= log2(N/lp); ++j){

        #pragma omp parallel for
        for( int k = 0; k < lp; ++k){
            int D_j = pow(2,j);
            int r = N/D_j;
            int starting = N/lp*k;
            int ending = N/lp*(k+1);

            for(int shift_index = starting; shift_index<ending; shift_index+=D_j){
                //calculate it
                for(int row = 0; row < D_j/2; ++row){
                    for(int column = 0; column < N; ++column){
                        mem2[shift_index+row][column] = mem1[shift_index/2+row][column] + mem1[shift_index/2 + N/2+row][column];
                        mem2[shift_index+D_j/2+row][column] = mem1[shift_index/2+row][(column+shift_index/2)%N] + mem1[shift_index/2+N/2+row][(column+(shift_index/2+N/2))%N];
                    }
                }
            }

        }
        complex<double> **prov= mem1;
        mem1 = mem2;
        mem2 = prov;
    }

    //for large J
    for( int j = log2(N/lp)+1; j <= log2(D); ++j){
        int D_j = pow(2,j);
        #pragma omp parallel for
        for( int k = 0; k < lp; ++k){
            double r = k*N/(lp*D_j);
            int shiftIndex = r;
            for(int row = D_j*(r-shiftIndex); row < D_j*((k+1)*N/(lp*D_j)-shiftIndex);++row){
                if(row < D_j/2){
                    for(int column = 0; column < N; ++column){
                        mem2[shiftIndex+row][column] = mem1[shiftIndex/2+row][column] + mem1[shiftIndex/2+N/2+row][column];
                    }
                }
                else{
                    for(int column = 0; column < N; ++column){
                        mem2[shiftIndex+row][column] = mem1[shiftIndex/2+row-D_j/2][(column+shiftIndex/2)%N] + mem1[shiftIndex/2+N/2+row-D_j/2][(column+(shiftIndex+N)/2)%N];
                    }
                }
            }
        }
        complex<double> **prov = mem1;
        mem1 = mem2;
        mem2 = prov;
    }

}

template<int N, int lp, int D>
void __calcFinalImageMatrix(complex<double> **mem1, complex<double> **mem2, complex<double> **imageMatrix, int machine_index){
    #pragma omp parallel for
    for(int k = 0; k < lp; ++k){
        int r = N/lp;
        for(int column = k*r; column < (k+1)*r; ++column){
            for(int row = 0; row < D; ++row){
                for(int shift_index = 0; shift_index < N; shift_index+=D){
                    imageMatrix[row][column] += mem1[shift_index+row][(shift_index*machine_index+column)%N];
                }
            }
        }
    }
}

}

#endif // PARA_H
