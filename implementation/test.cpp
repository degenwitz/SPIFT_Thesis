#include<vector>
#include<complex>
#include"visibilities.h"
#include<math.h>
#include<iostream>
#include<omp.h>

using namespace std;

template<int N, int lp, int D>
void parallelDoubleStep( vector<vector<visibility<N>>> parallelVisibilities, vector<complex<double>> W, complex<double> **imageMatrix, int machine_index ){
    int r = N/lp;

    complex<double> **mem1 = new complex<double>*[N];
    complex<double> **mem2 = new complex<double>*[N];

    //calculating shift-vectors and storing it in memory 1
    #pragma omp parallel for
    for( int k = 0; k < lp; ++k){

        for( int shift_index = k*r; shift_index < (k+1)*r; ++shift_index){
            mem1[shift_index] = new complex<double>[N];
            mem2[shift_index] = new complex<double>[N];

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

    //calculating doublestep
    //for small j
    for( int j = 1; j <= log2(D) && j <= log2(N/lp); ++j){

        #pragma omp parallel for
        for( int k = 0; k < lp; ++k){
            int D_j = pow(2,j);
            int r = N/D_j;

            for(int log_si = k*r; log_si<(k+1)*r; ++log_si){
                int shift_index = pow(2, log_si);

                //calculate it
                for(int row = 0; row < D_j/2; ++row){
                    for(int column = 0; column < N; ++column){
                        mem2[shift_index+row][column] = mem1[shift_index/2+row][column] + mem1[shift_index/2 + N/2+row][column];
                        mem2[shift_index+D_j/2+row][column] = mem1[shift_index/2+row][(column+shift_index/2)%N] + mem1[shift_index/2+row][column+((shift_index/2+N/2))%N];
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

    //adding it all to the image matrix
    #pragma omp parallel for
    for(int k = 0; k < lp; ++k){
        int r = N/lp;
        for(int column = k*r; column < (k+1)*r; ++column){
            for(int row = 0; row < N; ++row){
                for(int shift_index = 0; shift_index < N; shift_index+D){
                    imageMatrix[row][column] += mem1[shift_index+row][(shift_index*machine_index+column)%N;
                }
            }
        }
    }
}


int main(){

    const int N = 32;
    //calculating W
    vector<complex<double>> W (N);
    complex<double> W_0 = exp(complex<double>(0,2*M_PI/N));
    for( int i = 0; i < N; ++i){
        W[i] = pow(W_0,i);
    }


    //test
    visibility<N> a(complex<double>(4,3), 4,5);
    a.calc_isCS();
    a.calc_shift_index();
    vector<vector<visibility<N>>> b (N,vector<visibility<N>>(N,a));
    parallelDoubleStep<N, 4, 8>(b, W);
}
