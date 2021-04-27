#include"localSpiftNode.h"
#include"spiftCoordinator.h"

void spift(int argc, char** argv){
  // Initialisation
  MPI_Init(&argc, &argv);

  // Reading size and rank
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank == 0){
    spift_coordinator<mn::N>();
  } else {
    local_spfit_machine<mn::N, mn::N/mn::lp/mn::size>(rank-1);
  }

  // Finalisation

  MPI_Finalize();
}
