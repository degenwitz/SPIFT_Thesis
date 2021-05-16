#include"localSpiftNode.h"
#include"spiftCoordinator.h"

void spift(int argc, char** argv, int *plan, int plan_size, bool with_dic){
  // Initialisation
  MPI_Init(&argc, &argv);

  // Reading size and rank
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if(rank == 0){
    spift_coordinator<mn::N>(plan, plan_size);
  } else {
    if(with_dic){
        if(!mn::cheating){
            local_spfit_machine_with_dic<mn::N, mn::N/mn::lp/mn::size>(rank-1);
        } else if(rank ==mn::cheating_number){
            local_spfit_machine_with_dic<mn::N, mn::N/mn::lp/mn::size>(rank-1);
        }
    } else {
        if(!mn::cheating){
            local_spfit_machine_without_dic<mn::N, mn::N/mn::lp/mn::size>(rank-1);
        } else if(rank==mn::cheating_number){
            local_spfit_machine_without_dic<mn::N, mn::N/mn::lp/mn::size>(rank-1);
        }
    }
  }

  // Finalisation

  MPI_Finalize();
}
