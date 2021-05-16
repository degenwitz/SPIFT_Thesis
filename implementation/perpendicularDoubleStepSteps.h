#ifndef PERPEN_H
#define PERPEN_H

#include<vector>
#include<complex>
#include"visibilities.h"
#include<math.h>
#include<iostream>
#include<omp.h>

using namespace std;

namespace perd{

struct line{
    int ln;
    int shiftIndex;
    line(int l, int sI): ln(l), shiftIndex(sI) {}
    line(): ln(-1), shiftIndex(-1) {}
};

template <int N, int lp, int D>
complex<double>** __createShiftVectors(complex<double> mem, vector<vector<visibility<N>>> &perpendicularVisibilities, vector<complex<double>> &W, vector<line> required_lines);

template <int N, int lp, int D>
complex<double>** __calculateDoubleStep(complex<double> **mem1, complex<double> **mem2, vector<line> required_lines);

template<int N, int lp, int D>
void __calcFinalImageMatrix(complex<double> **mem1, complex<double> **mem2, complex<double> **imageMatrix, int machine_index);

template<int N>
void old_calc_required_lines(vector<line> &lines,int shiftIndex, int D_i, vector<vector<line>> &requiLines){
    requiLines[int(log2(D_i))].insert( requiLines[int(log2(D_i))].end(), lines.begin(), lines.end());
    if(D_i == 1){
        return;
    }

    vector<line> firstHalf;
    vector<line> secondHalf;
    for(int i = 0; i < lines.size(); ++i){
        int l = lines[i].ln;
        firstHalf.push_back(line(l,shiftIndex/2));
        if(shiftIndex != 0){
            firstHalf.push_back(line((l+shiftIndex/2)%N,shiftIndex/2));
        }
        secondHalf.push_back(line(l,(shiftIndex/2+N/2)%N));
        secondHalf.push_back(line((l+(shiftIndex+N)/2)%N,(shiftIndex/2+N/2)%N));
    }
    old_calc_required_lines<N>(firstHalf,shiftIndex/2, D_i/2, requiLines);
    old_calc_required_lines<N>(secondHalf,(shiftIndex+N)/2%N, D_i/2, requiLines);
}

template<int N>
void calc_required_lines(vector<line> &lines,int shiftIndex, int D, vector<vector<line>> &requiLines){
    int final_lines_index = int(log2(D));
    bool ***lines_to_calc = new bool**[final_lines_index+1];

    //set up first line
    lines_to_calc[final_lines_index] = new bool*[1];
    lines_to_calc[final_lines_index][0] = new bool[N];
    for(int i = 0; i < N; ++i){
        lines_to_calc[final_lines_index][0][i] = false;
    }
    for(int i = 0; i < lines.size(); ++i){
        lines_to_calc[final_lines_index][0][lines[i].ln] = true;
    }


    //calc all other lines
    for(int i = final_lines_index; i>0; --i){
        //setup
        lines_to_calc[i-1] = new bool*[int(pow(2,final_lines_index-(i-1)))];
        for(int j = 0; j < int(pow(2,final_lines_index-(i-1))); ++j){
            lines_to_calc[i-1][j] = new bool[N];
            for(int k = 0; k < N; ++k){
                lines_to_calc[i-1][j][k] = false;
            }
        }
        //calculate next ones
        for(int j = 0; j < int(pow(2,final_lines_index-(i))); ++j){
            int shiftIndex = pow(2,i)*j;
            int D_i_1 = pow(2,i-1);
            for(int k = 0; k < N; ++k){
                if(lines_to_calc[i][j][k] == true){
                    int l = k;
                    lines_to_calc[i-1][int(shiftIndex/2/D_i_1)][l] = true;
                    lines_to_calc[i-1][int(shiftIndex/2/D_i_1)][(l+shiftIndex/2)%N] = true;
                    lines_to_calc[i-1][int((((shiftIndex+N)/2)%N)/D_i_1)][l] = true;
                    lines_to_calc[i-1][int((((shiftIndex+N)/2)%N)/D_i_1)][(l+(shiftIndex+N)/2)%N] = true;
                }
            }
        }
    }
    //update the requiLines vector
    for(int i = 0; i < final_lines_index+1; ++i){
        int D_i = pow(2,i);
        for(int j = 0; j < int(pow(2,final_lines_index-(i))); ++j){
            int shiftIndex = D_i*j;
            for(int k = 0; k < N; ++k){
                if(lines_to_calc[i][j][k] == true){
                    requiLines[i].push_back(line(k,shiftIndex));
                }
            }
        }
    }
}

template<int N, int lp, int D>
void perpendicularDoubleStep( vector<vector<visibility<N>>> &perpendicularVisibilities, vector<complex<double>> &W, complex<double> **imageMatrix, int machine_index, vector<vector<line>> &required_lines ){

    cout << "starting with perpendicular" << endl;

    complex<double> **mem1 = new complex<double>*[N];
    complex<double> **mem2 = new complex<double>*[N];
    for( int i = 0; i < N; ++i){
        mem1[i] = new complex<double>[N];
        mem2[i] = new complex<double>[N];
    }

    double t_createShiftVectors = omp_get_wtime();
    //calculating shift-vectors
    complex<double>** mem_new = __createShiftVectors<N, lp, D>(mem2, perpendicularVisibilities,W,required_lines[0]);

    __safe_local_timer(machine_index,t_createShiftVectors,"create_col_shift_dic");


    if(mem2 == mem_new){
        mem2 = mem1;
        mem1 = mem_new;
    }

    double t_calculateImageMatrix = omp_get_wtime();
    //doublestep
    mem_new = __calculateDoubleStep<N, lp, D>(mem1, mem2, required_lines);

    cout << "calculated perpendicular doublestel" << endl;

    if(mem2 == mem_new){
        mem2 = mem1;
        mem1 = mem_new;
    }

    //update image-matrix
    __calcFinalImageMatrix<N, lp, D>(mem1, mem2, imageMatrix, machine_index);

    __safe_local_timer(machine_index, t_calculateImageMatrix, "calculated_col_doublestep");
    cout << "updated parallel image matrix" << endl;


    //deleting memory
    for(int i = 0; i < N; ++i){
        delete mem1[i];
        delete mem2[i];
    }
    delete mem1;
    delete mem2;
}

template <int N, int lp, int D>
complex<double>** __createShiftVectors(complex<double>** mem, vector<vector<visibility<N>>> &perpendicularVisibilities, vector<complex<double>> &W, vector<line> required_lines){
    int allLines = required_lines.size();
    int r = allLines/lp;

    omp_set_dynamic(0);
    omp_set_num_threads(mn::lp);
    #pragma omp parallel for shared(perpendicularVisibilities, required_lines)
    for( int k = 0; k < lp; ++k){
        for(int i = k*r; i < (k+1)*r && i < allLines; ++i){
            line l = required_lines[i];
            int shift_index = l.shiftIndex;
            int row = l.ln;
            if(perpendicularVisibilities[shift_index].size()<1){
                mem[row][shift_index] = complex<double>(0,0);
            }
            else{
                visibility<N> vis = perpendicularVisibilities[shift_index][0];
                mem[row][shift_index] = vis.vis * W[(vis.u*row)%N];
                for(int v = 1; v < perpendicularVisibilities[shift_index].size(); ++ v){
                    vis = perpendicularVisibilities[shift_index][v];
                    mem[row][shift_index] += vis.vis * W[(vis.u*row)%N];
                }
            }
        }
    }
    return mem;
}


template <int N, int lp, int D>
complex<double>** __calculateDoubleStep(complex<double> **mem1, complex<double> **mem2, vector<vector<line>> required_lines){

    for( int j = 1; j < log(N)+1; ++j){
        int allLines = required_lines[j].size();//<---------------------------------------------------------geht nicht wenn lp zu gross ist
        int D_i = int(pow(2,j));

        if(allLines >= lp){
            int r = allLines/lp;

            omp_set_dynamic(0);
            omp_set_num_threads(mn::lp);
            #pragma omp parallel for
            for( int k = 0; k < lp; ++k){
                for(int i = k*r; i < (k+1)*r && i < allLines; ++i){
                    line l = required_lines[j][i];
                    int shift_index = l.shiftIndex;
                    int row = l.ln;
                    for(int column = 0; column < D_i/2; ++column){
                        mem2[row][shift_index+column] = mem1[row][shift_index/2+column] + mem1[row][(shift_index+N)/2+column];
                    }
                    for(int column = D_i/2; column < D_i; ++column){
                        int col = column-D_i/2;
                        mem2[row][shift_index+column] = mem1[(row+shift_index/2)%N][shift_index/2+col] + mem1[((shift_index+N)/2+row)%N][((shift_index+N)/2+col)%N];
                    }

                }
            }
        } else {
            double r = allLines*1./lp;

            omp_set_dynamic(0);
            omp_set_num_threads(mn::lp);
            #pragma omp parallel for
            for( int k = 0; k < lp; ++k){
                line l = required_lines[j][int(r*k)];
                int firstCol = (r*k-int(r*k)) * D_i + 0.5;
                int lastCol = (r*(k+1)-int(r*k)) * D_i + 0.5;
                int shift_index = l.shiftIndex;
                int row = l.ln;
                for(int column = firstCol; column < D_i/2 && column < lastCol; ++column){
                    mem2[row][shift_index+column] = mem1[row][shift_index/2+column] + mem1[row][(shift_index+N)/2+column];
                }
                for(int column = (D_i/2<firstCol)?firstCol:D_i/2; column < D_i && column < lastCol; ++column){
                    int col = column-D_i/2;
                    mem2[row][shift_index+column] = mem1[(row+shift_index/2)%N][shift_index/2+col] + mem1[((shift_index+N)/2+row)%N][((shift_index+N)/2+col)%N];
                }

            }
        }
        complex<double> **prov = mem1;
        mem1 = mem2;
        mem2 = prov;
    }

    return mem1;
}

template<int N, int lp, int D>
void __calcFinalImageMatrix(complex<double> **mem1, complex<double> **mem2, complex<double> **imageMatrix, int machine_index){
    int r = N/lp;

    omp_set_dynamic(0);
    omp_set_num_threads(mn::lp);
    #pragma omp parallel for
    for( int k = 0; k < lp; ++k){
        for( int column = r*k; column < r*(k+1); ++column){
            for(int row = 0; row < D; ++row){
                imageMatrix[row][column] += mem1[D*machine_index+row][column];
            }
        }
    }
}
}

#endif
