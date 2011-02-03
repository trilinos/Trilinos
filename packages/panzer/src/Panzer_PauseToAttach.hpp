#ifndef PANZER_PAUSE_TO_ATTACH
#define PANZER_PAUSE_TO_ATTACH

namespace panzer {
  
  void pauseToAttach()
  {
    MPI_Comm mpicomm = MPI_COMM_WORLD;
    Teuchos::RCP<Teuchos::Comm<int> > comm = 
      Teuchos::createMpiComm<int>(Teuchos::rcp(new Teuchos::OpaqueWrapper<MPI_Comm>(mpicomm)));
    Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
    out.setShowProcRank(true);
    out.setOutputToRootOnly(-1);
    
    //out << "PID = " << getpid();
    if (comm->getRank() == 0)
      getchar();
    comm->barrier();
  }
  
}

#endif
