#include"visibilities.h"
#include"magickNumbers.h"
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <mpi.h>
#include<iostream>

template<int N>
std::vector<visibility<N>> get_random_datapoints(){
        std::vector<visibility<N>> list_of_visibilities;

        //get Data
        for(int i = 0; i < mn::visibilities_per_message; ++i){
            int u = rand()%mn::N;
            int v = rand()%mn::N;
            complex<double> val(rand()%40,rand()%40);
            visibility<N> vis(val,u,v);
            list_of_visibilities.push_back(vis);
        }

        return list_of_visibilities;
}

template<int N>
std::vector<visibility<N>> get_first_row_ones(){
    std::vector<visibility<N>> lst;
    lst.push_back(visibility<N>(complex<double>(1,0),0,0));
    lst.push_back(visibility<N>(complex<double>(1,0),1,1));
    lst.push_back(visibility<N>(complex<double>(1,0),1,1));
    lst.push_back(visibility<N>(complex<double>(1,0),1,1));
    return lst;
}

template<int N>
void spift_coordinator(){

    srand( (unsigned)time(NULL) );

    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    while(true){

        std::vector<visibility<N>> list_of_visibilities = get_first_row_ones<N>();

        //prepare data and insert into buffer
        int buf[mn::visibilities_per_message*6];
        for(int i = 0; i < list_of_visibilities.size() && i < mn::visibilities_per_message; ++i){
            int buf_i = i*6;
            list_of_visibilities[i].calc_isCS();
            list_of_visibilities[i].calc_shift_index();
            buf[buf_i+mn::VIS_ENCODING::u] = list_of_visibilities[i].u;
            buf[buf_i+mn::VIS_ENCODING::v] = list_of_visibilities[i].v;
            buf[buf_i+mn::VIS_ENCODING::isCS] = list_of_visibilities[i].isCS;
            buf[buf_i+mn::VIS_ENCODING::shift_index] = list_of_visibilities[i].shift_index;
            float vis_r_double = list_of_visibilities[i].vis.real();
            float vis_i_double = list_of_visibilities[i].vis.imag();
            buf[buf_i+mn::VIS_ENCODING::vis_r] = *(int*)&vis_r_double;
            buf[buf_i+mn::VIS_ENCODING::vis_i] = *(int*)&vis_i_double;
        }

        visibility<N> v = list_of_visibilities[0];


        //send
        for(int rank = 1; rank < size; ++rank){
            MPI_Send(buf,mn::visibilities_per_message*6, MPI_INT, rank, mn::performDoubleStep, MPI_COMM_WORLD);
        }
    }
}
