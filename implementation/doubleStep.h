#include"magickNumbers.h"
#include"localDoubleStep.h"
#include"doubleStepCoordinator.h"

void double_step(int argc, char** argv){
  // Initialisation
  MPI_Init(&argc, &argv);

  // Reading size and rank
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank == 0){
    double_step_coordinator<mn::N>();
  } else {
    local_double_step<mn::N>(rank-1);
  }

  // Finalisation

  MPI_Finalize();
}
