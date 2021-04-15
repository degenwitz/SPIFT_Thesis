#include<vector>
#include<complex>
#include<math.h>
#include<iostream>
#include<omp.h>
#include"visibilities.h"
#include"parrallelDoubleStepSteps.h"
#include"perpendicularDoubleStepSteps.h"
#include"localDoubleStep.h"
#include"doubleStep.h"

using namespace std;


int main(int argc, char** argv){
    double_step(argc, argv);
}


