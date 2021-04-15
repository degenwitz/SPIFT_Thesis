#ifndef VISIBILITY_H
#define VISIBILITY_H

#include<complex>
#include<algorithm>

template<int N>
struct visibility{
public:
    int u,v;
    std::complex<double> vis;
    bool isCS;
    int shift_index;

    visibility(std::complex<double> vis,int u,int v): u(u), v(v), vis(vis) { }

    visibility(std::complex<double> vis,int u,int v, bool isCS, int shift_index): u(u), v(v), vis(vis), isCS(isCS), shift_index(shift_index) { }

    void calc_isCS(){
        isCS = (v == 0) || (u%2==1 && v%2 == 0) || (v%2 == 0 && std::__gcd(u,N) < std::__gcd(v,N) );
    }

    void calc_shift_index(){
        if( u == 0 || v == 0) {shift_index = 0; return; }
        if(isCS){
            for( int j = 0; j < N; ++j){
                if(v==(j*u)%N){ shift_index = j; return; }
            }
        } else {
            for( int k = 0; k < N; ++k){
                if(u == (k*v)%N){ shift_index = k; return; }
            }
        }
    }
};
#endif // VISIBILITY_H
