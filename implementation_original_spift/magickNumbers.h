#ifndef MAGICKNUMBERS_H
#define MAGICKNUMBERS_H

namespace mn{
const int number_of_runs=50;
const int N = 8192;
const int lp = 32;
const bool cheating = true;
const int cheating_number = 2;
const int size = 2;
const int MPI_host = 0;
const int visibilities_per_message = 666*8*8;
enum MPI_COMMUNICATION_{performDoubleStep=0};
namespace VIS_ENCODING{
    enum VIS_ENCODING{u=0, v=1, vis_r =2, vis_i=3, isCS=4, shift_index=5};
}
namespace COMMANDS{
    enum COMMANDS{calc=0, stop=1, safe_results=2, start_timer=3, safe_timer=4};
}
}
#endif
