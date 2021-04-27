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
void calc_required_lines(vector<line> &lines,int shiftIndex, int D_i, vector<vector<line>> &requiLines){
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
    calc_required_lines<N>(firstHalf,shiftIndex/2, D_i/2, requiLines);
    calc_required_lines<N>(secondHalf,(shiftIndex+N)/2%N, D_i/2, requiLines);
}

template<int N, int lp, int D>
void perpendicularDoubleStep( vector<vector<visibility<N>>> &perpendicularVisibilities, vector<complex<double>> &W, complex<double> **imageMatrix, int machine_index ){

    cout << " starting perpendicular " << endl;

    complex<double> **mem1 = new complex<double>*[N];
    complex<double> **mem2 = new complex<double>*[N];
    for( int i = 0; i < N; ++i){
        mem1[i] = new complex<double>[N];
        mem2[i] = new complex<double>[N];
    }


    //calculating required lines for each shift-matrix
    vector<vector<line>> required_lines(log2(N)+1);
    vector<line> lines;
    for(int i = D*machine_index; i < D*machine_index+D; ++i){
        lines.push_back(line(i,0));
    }
    calc_required_lines<N>(lines, 0, N, required_lines);

    //calculating shift-vectors
    complex<double>** mem_new = __createShiftVectors<N, lp, D>(mem2, perpendicularVisibilities,W,required_lines[0]);

    if(mem2 == mem_new){
        mem2 = mem1;
        mem1 = mem_new;
    }

    //doublestep
    mem_new = __calculateDoubleStep<N, lp, D>(mem1, mem2, required_lines);

    if(mem2 == mem_new){
        mem2 = mem1;
        mem1 = mem_new;
    }

    //update image-matrix
    __calcFinalImageMatrix<N, lp, D>(mem1, mem2, imageMatrix, machine_index);

    cout << " done with perpendicular " << endl;

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
    #pragma omp parallel for
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
