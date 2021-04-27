#ifndef LOCAL_Spift_Node
#define LOCAL_Spift_Node

#include"magickNumbers.h"
#include"visibilities.h"
#include <mpi.h>
#include<vector>
#include<complex>
#include<omp.h>

using namespace std;

template<int N, int D>
void printImageMatrix( int rank, complex<double>** image){
    for( int row = 0; row < D; ++row){
        cout << "\n[ ";
        for(int column = 0; column < N; ++column){
            cout << image[row][column] << ", ";
        }
        cout << "]";
    }
    cout << "-> " << rank << endl;
}


template<int N, int D>
void __local_spift_node(int rank, complex<double> **imageMatrix, vector<vector<visibility<N>>> &rowVisibilities, vector<vector<visibility<N>>> &colVisibilities, vector<complex<double>> &W){
    complex<double> shift_vectors_col [N][N];
    complex<double> shift_vectors_row [N][N];


    //create row shift vectors
    for( int shift_index = 0; shift_index < N; ++shift_index){
        //set up with initial values (0 if no visibilities provided
        if(rowVisibilities[shift_index].size() == 0){
            for(int i = 0; i < N; ++i){
                shift_vectors_row[shift_index][i] = complex<double>(0,0);
            }
        } else {
            visibility<N> v = rowVisibilities[shift_index][0];
            for( int i = 0; i < N; ++i){
                shift_vectors_row[shift_index][i] = v.vis*W[(i*v.v)%N];
            }
            for(unsigned int j = 1; j < rowVisibilities[shift_index].size(); ++ j){
                visibility<N> vis = rowVisibilities[shift_index][j];
                for( int i = 0; i < N; ++i){
                    shift_vectors_row[shift_index][i] += vis.vis*W[(i*v.v)%N];
                }
            }
        }
    }

    //create col shift vectors
    for( int shift_index = 0; shift_index < N; ++shift_index){
        //set up with initial values (0 if no visibilities provided
        if(colVisibilities[shift_index].size() == 0){
            for(int i = 0; i < N; ++i){
                shift_vectors_col[shift_index][i] = complex<double>(0,0);
            }
        } else {
            visibility<N> v = colVisibilities[shift_index][0];
            for( int i = 0; i < N; ++i){
                shift_vectors_col[shift_index][i] = v.vis*W[(i*v.v)%N];
            }
            for(unsigned int j = 1; j < colVisibilities[shift_index].size(); ++ j){
                visibility<N> vis = colVisibilities[shift_index][j];
                for( int i = 0; i < N; ++i){
                    shift_vectors_col[shift_index][i] += vis.vis*W[(i*v.v)%N];
                }
            }
        }
    }

    //add row-shifts
    for( int shift_index = 0; shift_index < N; ++ shift_index){
        for(int j = 0; j < D; ++j){
            int startindex = (shift_index*(j+(rank*D)))%N;
            for( int k = 0; k < N; ++k){
                int idx = (startindex+k)%N;
                imageMatrix[j][k] += shift_vectors_row[shift_index][idx];
            }
        }
    }

    //add column-shift
    for( int shift_index = 0; shift_index < N; ++shift_index){
        for( int k = 0; k < N; ++k){
            int startindex = ((rank*D)+(shift_index*k))%N;
            for(int j = 0; j < D; ++j){
                int idx = (startindex+j)%N;
                imageMatrix[j][k] += shift_vectors_col[shift_index][idx];
            }
        }
    }
}

template<int N, int D>
void local_spfit_machine(int rank){

    //create image matrices
    vector<complex<double>**> imageMatrices(mn::lp);
    for(int machine = 0; machine < mn::lp; ++machine){
        complex<double> **imageMatrix = new complex<double>*[D];
        for( int i = 0; i < D; ++i){
            imageMatrix[i] = new complex<double>[N];
            for( int j = 0; j < N; ++j){
                imageMatrix[i][j] = complex<double>(0,0);
            }
        }
        imageMatrices[machine] = imageMatrix;
    }


    //calculating W
    vector<complex<double>> W (N);
    complex<double> W_0 = exp(complex<double>(0,2*M_PI/N));
    for( int i = 0; i < N; ++i){
        W[i] = pow(W_0,i);
    }


    while(true){
        int buf[ mn::visibilities_per_message*6 ]; //note here v, u, vis_i, vis_r, isCS, shift_index
        int res = MPI_Recv( &buf, mn::visibilities_per_message*6, MPI_INT, mn::MPI_host, mn::performDoubleStep, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        vector<vector<visibility<N>>> columShifts(N);
        vector<vector<visibility<N>>> rowShifts(N);

        for( int i = 0; i < N; i += 6){
            double vis_r = *(float*)&buf[i+mn::VIS_ENCODING::vis_r];
            double vis_i = *(float*)&buf[i+mn::VIS_ENCODING::vis_i];
            int u = buf[i+mn::VIS_ENCODING::u];
            int v = buf[i+mn::VIS_ENCODING::v];
            int shift_index = buf[i+mn::VIS_ENCODING::shift_index];
            bool isCS = buf[i+mn::VIS_ENCODING::isCS];
            visibility<N> vis(complex<double>(vis_r, vis_i), u, v, isCS, shift_index);
            if(isCS){
                columShifts[shift_index].push_back(vis);
            } else {
                rowShifts[shift_index].push_back(vis);
            }
        }

        #pragma omp parallel for
        for( int k = 0; k < mn::lp; ++k){
            __local_spift_node<N,D>(rank*mn::lp+k, imageMatrices[k], rowShifts, columShifts, W);
        }
    }


    //deleting memory
    for(int machine = 0; machine < mn::lp; ++machine){
        complex<double> **imageMatrix = new complex<double>*[D];
        for( int i = 0; i < D; ++i){
            delete imageMatrix[i];
        }
        delete imageMatrices[machine];
    }


}

#endif // LOCAL_Spift_Node

