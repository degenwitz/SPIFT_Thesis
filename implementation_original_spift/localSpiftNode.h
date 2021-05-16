#ifndef LOCAL_Spift_Node
#define LOCAL_Spift_Node

#include"magickNumbers.h"
#include"visibilities.h"
#include <mpi.h>
#include<vector>
#include<complex>
#include<omp.h>
#include <fstream>
#include "measure.h"
#include "shiftDictionary.h"

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
void __local_spift_node_with_dic(int rank, complex<double> **imageMatrix,
    complex<double> **shift_vectors_row,
    complex<double> **shift_vectors_col){

    double t_calculateImageMatrix = omp_get_wtime();

    //add row-shifts
    for( int shift_index = 0; shift_index < N; ++shift_index){
        for(int j = 0; j < D; ++j){
            int startindex = (shift_index*(j+(rank*D)))%N;
            for( int k = 0; k < N; ++k){
                int idx = (startindex+k)%N;
                imageMatrix[j][k] += shift_vectors_row[shift_index][idx];
            }
        }
    }

    __safe_local_timer(rank, t_calculateImageMatrix, "calculated_spift_row_with_dic");
    t_calculateImageMatrix = omp_get_wtime();
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

    __safe_local_timer(rank, t_calculateImageMatrix, "calculated_spift_col_with_dic");

}


template<int N, int D>
void __local_spift_node_without_dic(int rank, complex<double> **imageMatrix,
vector<complex<double>> W,
    visibility<N> vis){

    //add row-shifts
    if(!vis.isCS){
        //calc
        int shift_index = vis.shift_index;
        vector<complex<double>> q(N);
        for( int i = 0; i < N; ++i){
            q[i] += vis.vis*W[(i*vis.v)%N];
        }


        for(int j = 0; j < D; ++j){
            int startindex = (shift_index*(j+(rank*D)))%N;
            for( int k = 0; k < N; ++k){
                int idx = (startindex+k)%N;
                imageMatrix[j][k] += q[idx];
            }
        }
    } else {
        //calc
        int shift_index = vis.shift_index;
        vector<complex<double>> q(N);
        for( int i = 0; i < N; ++i){
            q[i] += vis.vis*W[(i*vis.u)%N];
        }

        for( int k = 0; k < N; ++k){
            int startindex = ((rank*D)+(shift_index*k))%N;
            for(int j = 0; j < D; ++j){
                int idx = (startindex+j)%N;
                imageMatrix[j][k] += q[idx];
            }
        }
    }
}

template<int N, int D>
void local_spfit_machine_with_dic(int rank){

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

    //measuring
    double t;

    cout << " done with setup " << endl;

    complex<double> ** shift_vectors_col = new complex<double>*[N];
    complex<double> ** shift_vectors_row = new complex<double>*[N];

    for(int i = 0; i< N; ++i){
        shift_vectors_col[i] = new complex<double>[N];
        shift_vectors_row[i] = new complex<double>[N];
    }


    while(true){
        int buf[ mn::visibilities_per_message*6+1 ]; //note here v, u, vis_i, vis_r, isCS, shift_index
        int res = MPI_Recv( &buf, mn::visibilities_per_message*6+1, MPI_INT, mn::MPI_host, mn::performDoubleStep, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (buf[mn::visibilities_per_message*6] == mn::COMMANDS::stop){
            cout << "stopping" << endl;
            break;
        } else if (buf[mn::visibilities_per_message*6] == mn::COMMANDS::safe_results){
            __safe_local_result<N,D>(rank, imageMatrices);
            continue;
        } else if(buf[mn::visibilities_per_message*6] == mn::COMMANDS::safe_timer){
            __safe_local_timer(rank, t, "time_per_run_with_dic");
            cout << " safed time " << endl;
            continue;
        } else if(buf[mn::visibilities_per_message*6] == mn::COMMANDS::start_timer){
            t = omp_get_wtime();
            continue;
        }


        vector<vector<visibility<N>>> columShifts(N);
        vector<vector<visibility<N>>> rowShifts(N);

        for( int i = 0; i < mn::visibilities_per_message*6; i += 6){
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

        double t_createShiftVectors = omp_get_wtime();
        calculate_col_shift_dictionary(shift_vectors_col, columShifts, W);
        __safe_local_timer(rank,t_createShiftVectors,"create_col_shift_dic_for_spift");

        t_createShiftVectors = omp_get_wtime();
        calculate_row_shift_dictionary(shift_vectors_row, rowShifts, W);
        __safe_local_timer(rank,t_createShiftVectors,"create_row_shift_dic_for_spift");

        omp_set_dynamic(0);
        omp_set_num_threads(mn::lp);
        #pragma omp parallel for
        for( int k = 0; k < mn::lp; ++k){
            __local_spift_node_with_dic<N,D>(rank*mn::lp+k, imageMatrices[k], shift_vectors_row, shift_vectors_col);
        }
    }


    //deleting memory
    for(int machine = 0; machine < mn::lp; ++machine){
        for( int i = 0; i < D; ++i){
            delete imageMatrices[machine][i];
        }
        delete imageMatrices[machine];
    }
    for(int i = 0; i< N; ++i){
        delete shift_vectors_col[i];
        delete shift_vectors_row[i];
    }
    delete shift_vectors_col;
    delete shift_vectors_row;


}

template<int N, int D>
void local_spfit_machine_without_dic(int rank){

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

    //measuring
    double t;

    cout << " done with setup " << endl;


    while(true){
        int buf[ mn::visibilities_per_message*6+1 ]; //note here v, u, vis_i, vis_r, isCS, shift_index
        int res = MPI_Recv( &buf, mn::visibilities_per_message*6+1, MPI_INT, mn::MPI_host, mn::performDoubleStep, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        if (buf[mn::visibilities_per_message*6] == mn::COMMANDS::stop){
            cout << "stopping" << endl;
            break;
        } else if (buf[mn::visibilities_per_message*6] == mn::COMMANDS::safe_results){
            __safe_local_result<N,D>(rank, imageMatrices);
            continue;
        } else if(buf[mn::visibilities_per_message*6] == mn::COMMANDS::safe_timer){
            __safe_local_timer(rank, t, "time_per_run_without_dic");
            cout << " safed time " << endl;
            continue;
        } else if(buf[mn::visibilities_per_message*6] == mn::COMMANDS::start_timer){
            t = omp_get_wtime();
            continue;
        }


        vector<visibility<N>> visibilities();



        double t_createShiftVectors = omp_get_wtime();

        for( int i = 0; i < mn::visibilities_per_message*6; i += 6){
            double vis_r = *(float*)&buf[i+mn::VIS_ENCODING::vis_r];
            double vis_i = *(float*)&buf[i+mn::VIS_ENCODING::vis_i];
            int u = buf[i+mn::VIS_ENCODING::u];
            int v = buf[i+mn::VIS_ENCODING::v];
            int shift_index = buf[i+mn::VIS_ENCODING::shift_index];
            bool isCS = buf[i+mn::VIS_ENCODING::isCS];
            visibility<N> vis(complex<double>(vis_r, vis_i), u, v, isCS, shift_index);

            omp_set_dynamic(0);
            omp_set_num_threads(mn::lp);
            #pragma omp parallel for
            for( int k = 0; k < mn::lp; ++k){
                __local_spift_node_without_dic<N,D>(rank*mn::lp+k, imageMatrices[k], W,vis);
            }
        }

        __safe_local_timer(rank,t_createShiftVectors,"spift_without_dic");
    }


    //deleting memory
    for(int machine = 0; machine < mn::lp; ++machine){
        for( int i = 0; i < D; ++i){
            delete imageMatrices[machine][i];
        }
        delete imageMatrices[machine];
    }


}

#endif // LOCAL_Spift_Node

