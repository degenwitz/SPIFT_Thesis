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
#include"magickNumbers.h"

using namespace std;


int main(int argc, char** argv){
    const int plan_size = mn::number_of_runs*3;
    int plan [plan_size];
    for(int i = 0; i < plan_size; i+= 3){
        plan[i] = mn::COMMANDS::start_timer;
        plan[i+1] = mn::COMMANDS::calc;
        plan[i+2] = mn::COMMANDS::safe_timer;
    }
    double_step_planned(argc, argv, plan, plan_size);
    //double_step(argc, argv);
}


