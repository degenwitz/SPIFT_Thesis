#include"spift.h"


int main(int argc, char** argv){
    const int plan_size = mn::number_of_runs*3;
    int plan [plan_size];
    for(int i = 0; i < plan_size; i+= 3){
        plan[i] = mn::COMMANDS::start_timer;
        plan[i+1] = mn::COMMANDS::calc;
        plan[i+2] = mn::COMMANDS::safe_timer;
    }
    spift(argc, argv, plan, plan_size, true);
}


