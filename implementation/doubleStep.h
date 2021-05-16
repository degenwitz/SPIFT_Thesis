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

void double_step_planned(int argc, char** argv, int *plan, int plan_size){
  // Initialisation
  MPI_Init(&argc, &argv);

  // Reading size and rank
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank == mn::MPI_host){
    std::cout << "one coordinator" << endl;
    double_step_coordinator_planned<mn::N>(plan, plan_size);
  } else if(!mn::cheating) {
    local_double_step<mn::N>(rank-1);
  } else if(rank ==mn::cheating_number){
    local_double_step<mn::N>(rank-1);
  }

  // Finalisation

  MPI_Finalize();
}
