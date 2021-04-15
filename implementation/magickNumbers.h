#ifndef MAGICKNUMBERS_H
#define MAGICKNUMBERS_H

namespace mn{

const int N = 8;
const int lp = 4;
const int size = 2;
const int MPI_host = 0;
const int visibilities_per_message = N*3/2;
enum MPI_COMMUNICATION_{performDoubleStep=0};
namespace VIS_ENCODING{
    enum VIS_ENCODING{u=0, v=1, vis_r =2, vis_i=3, isCS=4, shift_index=5};
}
}
#endif
