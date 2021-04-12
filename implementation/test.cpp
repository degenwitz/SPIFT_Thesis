#include<vector>
#include<complex>
#include"visibilities.h"
#include<math.h>
#include<iostream>
#include<omp.h>

using namespace std;

template<int N, int D>
void printImage(complex<double> **imageMatrix){
    for(int row = 0; row < D; ++row){
        cout << "[";
        for(int column = 0; column < N; ++column){
            cout << imageMatrix[row][column] << ", ";
        }
        cout << "]" << endl;
    }
}

template<int N, int lp, int D>
void parallelDoubleStep( vector<vector<visibility<N>>> parallelVisibilities, vector<complex<double>> W, complex<double> **imageMatrix, int machine_index ){
    cout << " is in parallel " << endl;
    int r = N/lp;

    complex<double> **mem1 = new complex<double>*[N];
    complex<double> **mem2 = new complex<double>*[N];

    cout << " created memoriy-slots " << endl;

    //calculating shift-vectors and storing it in memory 1
    #pragma omp parallel for
    for( int k = 0; k < lp; ++k){

        for( int shift_index = k*r; shift_index < (k+1)*r; ++shift_index){
            cout<< "one time " << endl;
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
    cout << "calculated parallel-shift vectors: " << endl;
    printImage<N,N>(mem1);

    //calculating doublestep
    //for small j
    for( int j = 1; j <= log2(D) && j <= log2(N/lp); ++j){

        //#pragma omp parallel for
        for( int k = 0; k < lp; ++k){
            int D_j = pow(2,j);
            int r = N/D_j;

            for(int shift_index = k*r; shift_index<(k+1)*r; shift_index+=D_j){


                //calculate it
                for(int row = 0; row < D_j/2; ++row){
                    for(int column = 0; column < N; ++column){
                        //cout << "row: "<< row <<"shift index: " << shift_index <<"[ ]"<< shift_index+row << endl;
                        mem2[shift_index+row][column] = mem1[shift_index/2+row][column] + mem1[shift_index/2 + N/2+row][column];
                        mem2[shift_index+D_j/2+row][column] = mem1[shift_index/2+row][(column+shift_index/2)%N] + mem1[shift_index/2+N/2+row][(column+(shift_index/2+N/2))%N];
                    }
                }
            }

        }
        complex<double> **prov= mem1;
        mem1 = mem2;
        mem2 = prov;
        cout << "parallel step j: " << j << endl;
        printImage<N,N>(mem1);
    }

    cout << "done with small j's" << endl;

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

    cout << "done with big j's" << endl;

    //adding it all to the image matrix
    #pragma omp parallel for
    for(int k = 0; k < lp; ++k){
        cout << "-----------------k= " << k << " -------------------" << endl;
        int r = N/lp;
        for(int column = k*r; column < (k+1)*r; ++column){
            for(int row = 0; row < D; ++row){
                for(int shift_index = 0; shift_index < N; shift_index+=D){
                    imageMatrix[row][column] += mem1[shift_index+row][(shift_index*machine_index+column)%N];
                }
            }
        }
    }

    cout << "done with parallel shift-indexes" << endl;
}


struct line{
    int ln;
    int shiftIndex;
    line(int l, int sI): ln(l), shiftIndex(sI) {}
    line(): ln(-1), shiftIndex(-1) {}
};

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
void perpendicularDoubleStep( vector<vector<visibility<N>>> perpendicularVisibilities, vector<complex<double>> W, complex<double> **imageMatrix, int machine_index ){
    complex<double> **mem1 = new complex<double>*[N];
    complex<double> **mem2 = new complex<double>*[N];
    for( int i = 0; i < N; ++i){
        mem1[i] = new complex<double>[N];
        mem2[i] = new complex<double>[N];
    }

    cout << "setup memory" << endl;

    //calculating required lines for each shift-matrix
    vector<vector<line>> required_lines(log2(N)+1);
    vector<line> lines;
    for(int i = D*machine_index; i < D*machine_index+D; ++i){
        lines.push_back(line(i,0));
    }
    calc_required_lines<N>(lines, 0, N, required_lines);

    cout << "found required lines" << endl;

    //calculating shift-vectors
    int allLines = required_lines[0].size();
    int r = allLines/lp;
    #pragma omp parallel for
    for( int k = 0; k < lp; ++k){
        for(int i = k*r; i < (k+1)*r && i < allLines; ++i){
            line l = required_lines[0][i];
            int shift_index = l.shiftIndex;
            int row = l.ln;
            cout << "row: " << row << " shift_index: " << shift_index << endl;
            if(perpendicularVisibilities[shift_index].size()<1){
                mem1[row][shift_index] = complex<double>(0,0);
            }
            else{
                visibility<N> vis = perpendicularVisibilities[shift_index][0];
                mem1[row][shift_index] = vis.vis * W[(vis.u*row)%N];
                for(int v = 1; v < perpendicularVisibilities[shift_index].size(); ++ v){
                    vis = perpendicularVisibilities[shift_index][v];
                    mem1[row][shift_index] += vis.vis * W[(vis.u*row)%N];
                }
            }
        }
    }
    complex<double> **prov = mem1;
    mem1 = mem2;
    mem2 = prov;
    cout << "calculated shift-vectors for perpendicular: " << endl;
    printImage<N,N>(mem2);



    //doublestep
    for( int j = 1; j < log(N)+1; ++j){

        int allLines = required_lines[j].size();
        int r = allLines/lp;
        int D_i = int(pow(2,j));

        #pragma omp parallel for
        for( int k = 0; k < lp; ++k){
            for(int i = k*r; i < (k+1)*r && i < allLines; ++i){
                line l = required_lines[0][i];
                int shift_index = l.shiftIndex;
                int row = l.ln;
                for(int column = 0; column < D_i/2; ++column){
                    mem1[row][shift_index+column] = mem2[row][shift_index/2+column] + mem2[row][(shift_index+N)/2+column];
                }
                for(int column = D_i/2; column < D_i; ++column){
                    int col = column-D_i/2;
                    mem1[row][shift_index+column] = mem2[(row+shift_index/2)%N][shift_index/2+col] + mem2[((shift_index+N)/2+row)%N][((shift_index+N)/2+col)%N];
                }

            }
        }
        complex<double> **prov = mem1;
        mem1 = mem2;
        mem2 = prov;
        cout << "after step j: " << j << endl;
        printImage<N,N>(mem2);
    }

    //update image-matrix
    r = N/lp;
    #pragma omp parallel for
    for( int k = 0; k < lp; ++k){
        for( int column = r*k; column < r*(k+1); ++column){
            for(int row = 0; row < D; ++row){
                imageMatrix[row][column] += mem2[D*machine_index+row][column];
            }
        }
    }

}

template<int N>
void addValues(visibility<N> &v,vector<vector<visibility<N>>> &parrallelShiftVectors, vector<vector<visibility<N>>> &perpendicularShiftVectors ){
v.calc_isCS();
v.calc_shift_index();
    if(v.isCS){
        perpendicularShiftVectors[v.shift_index] = vector<visibility<N>> (1, v);
    } else {
        parrallelShiftVectors[v.shift_index] = vector<visibility<N>> (1, v);
    }
}


int main(){
    const int N = 4;
    const int D = 4;
    //calculating W
    vector<complex<double>> W (N);
    complex<double> W_0 = exp(complex<double>(0,2*M_PI/N));
    for( int i = 0; i < N; ++i){
        W[i] = pow(W_0,i);
    }
    //create image matrix
    complex<double> **imageMatrix = new complex<double>*[D];
    for( int i = 0; i < D; ++i){
        imageMatrix[i] = new complex<double>[N];
        for( int j = 0; j < N; ++j){
            imageMatrix[i][j] = complex<double>(0,0);
        }
    }

    //test
    vector<vector<visibility<N>>> parrallelShiftVectors (N);
    vector<vector<visibility<N>>> perpendicularShiftVectors (N);
    visibility<N> a1(complex<double>(1,0), 2,2);
    addValues(a1, parrallelShiftVectors, perpendicularShiftVectors);
    visibility<N> a2(complex<double>(1,0), 1,3);
    addValues(a2, parrallelShiftVectors, perpendicularShiftVectors);

    parallelDoubleStep<N, 2, D>(parrallelShiftVectors, W, imageMatrix, 0);
    perpendicularDoubleStep<N, 2, D>(perpendicularShiftVectors, W, imageMatrix, 0);
    cout << "image-matrix: " << endl;
    printImage<N,D>(imageMatrix);
}
