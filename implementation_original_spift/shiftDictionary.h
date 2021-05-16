#ifndef SHIFT_DIC_H
#define SHIFT_DIC_H

#include<complex>
#include<vector>
#include"magickNumbers.h"

using namespace std;

template<int N>
void calculate_row_shift_dictionary(complex<double> **row_dic, vector<vector<visibility<N>>> &rowVisibilities, vector<complex<double>> &W){

    int r = N/mn::lp;
    //create row shift vectors
    omp_set_dynamic(0);
    omp_set_num_threads(mn::lp);
    #pragma omp parallel for shared(rowVisibilities)
    for( int k = 0; k < mn::lp; ++k){

        for( int shift_index = k*r; shift_index < (k+1)*r; ++shift_index){

            //set up with initial values (0 if no visibilities provided
            if(rowVisibilities[shift_index].size() == 0){
                for(int i = 0; i < N; ++i){
                    row_dic[shift_index][i] = complex<double>(0,0);
                }
            } else {
                visibility<N> v = rowVisibilities[shift_index][0];
                for( int i = 0; i < N; ++i){
                    row_dic[shift_index][i] = v.vis*W[(i*v.v)%N];
                }
                for(unsigned int j = 1; j < rowVisibilities[shift_index].size(); ++ j){
                    visibility<N> vis = rowVisibilities[shift_index][j];
                    for( int i = 0; i < N; ++i){
                        row_dic[shift_index][i] += vis.vis*W[(i*v.v)%N];
                    }
                }
            }
        }
    }
}

template<int N>
void calculate_col_shift_dictionary(complex<double> **col_dic, vector<vector<visibility<N>>> &colVisibilities, vector<complex<double>> &W){

    int r = N/mn::lp;
    //create row shift vectors
    omp_set_dynamic(0);
    omp_set_num_threads(mn::lp);
    #pragma omp parallel for shared(colVisibilities)
    for( int k = 0; k < mn::lp; ++k){

        for( int shift_index = k*r; shift_index < (k+1)*r; ++shift_index){

            //set up with initial values (0 if no visibilities provided
            if(colVisibilities[shift_index].size() == 0){
                for(int i = 0; i < N; ++i){
                    col_dic[shift_index][i] = complex<double>(0,0);
                }
            } else {
                visibility<N> v = colVisibilities[shift_index][0];
                for( int i = 0; i < N; ++i){
                    col_dic[shift_index][i] = v.vis*W[(i*v.u)%N];
                }
                for(unsigned int j = 1; j < colVisibilities[shift_index].size(); ++ j){
                    visibility<N> vis = colVisibilities[shift_index][j];
                    for( int i = 0; i < N; ++i){
                        col_dic[shift_index][i] += vis.vis*W[(i*v.u)%N];
                    }
                }
            }
        }
    }

}


#endif //SHIFT_DIC
