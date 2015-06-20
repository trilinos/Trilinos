#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Pike_MultiphysicsDistributor.hpp"

namespace pike {

  TEUCHOS_UNIT_TEST(MultiphysicsDistributor, basic)
  {
    Teuchos::RCP<Teuchos::MpiComm<int> > globalComm = 
      Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));

    // Run this on 6 processes only
    TEST_EQUALITY(globalComm->getSize(), 6);

    typedef pike::MultiphysicsDistributor::ApplicationIndex AppIndex;
    typedef pike::MultiphysicsDistributor::TransferIndex TransIndex;
    AppIndex CTF;
    AppIndex Insilico;
    AppIndex Peregrine;
    TransIndex CTF_Insilico;
    TransIndex CTF_Peregrine;
    TransIndex Insilico_Peregrine;
    
    pike::MultiphysicsDistributor dist;
    {
      CTF = dist.addApplication("CTF",0,0);
      Insilico = dist.addApplication("Insilico",1,3);
      std::vector<int> peregrineRanks;
      peregrineRanks.push_back(4);
      peregrineRanks.push_back(5);
      Peregrine = dist.addApplication("Peregrine",peregrineRanks); // vector of ranks
      
      // try to add application already registered (check for
      // duplicate name)
      TEST_THROW(dist.addApplication("CTF",1,1), std::logic_error);

      CTF_Insilico = dist.addTransfer("C_TO_I: ",CTF,Insilico);
      
      std::vector<pike::MultiphysicsDistributor::ApplicationIndex> ctfPeregrine;
      ctfPeregrine.push_back(CTF);
      ctfPeregrine.push_back(Peregrine);
      CTF_Peregrine = dist.addTransfer("C_TO_P: ",ctfPeregrine);
      
      Insilico_Peregrine = dist.addTransfer("I_TO_P: ",Insilico,Peregrine);
      
      // try to add transfer already registered (check for duplicate
      // name)
      TEST_THROW(dist.addTransfer("C_TO_I: ",CTF,Insilico), std::logic_error);
      TEST_THROW(dist.addTransfer("C_TO_P: ",ctfPeregrine), std::logic_error);

      dist.setup(globalComm,true);
    }

    if (globalComm->getRank() == 0) {
      TEST_EQUALITY(dist.appExistsOnProcess(CTF),true);
      TEST_EQUALITY(dist.appExistsOnProcess(Insilico),false);
      TEST_EQUALITY(dist.appExistsOnProcess(Peregrine),false);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Insilico),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Peregrine),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(Insilico_Peregrine),false);
    }
    if ( (globalComm->getRank() == 1) || 
	 (globalComm->getRank() == 2) ||
	 (globalComm->getRank() == 3) ) {
      TEST_EQUALITY(dist.appExistsOnProcess(CTF),false);
      TEST_EQUALITY(dist.appExistsOnProcess(Insilico),true);
      TEST_EQUALITY(dist.appExistsOnProcess(Peregrine),false);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Insilico),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Peregrine),false);
      TEST_EQUALITY(dist.transferExistsOnProcess(Insilico_Peregrine),true);
    }
    if ( (globalComm->getRank() == 4) ||
	 (globalComm->getRank() == 5) ) {
      TEST_EQUALITY(dist.appExistsOnProcess(CTF),false);
      TEST_EQUALITY(dist.appExistsOnProcess(Insilico),false);
      TEST_EQUALITY(dist.appExistsOnProcess(Peregrine),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Insilico),false);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Peregrine),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(Insilico_Peregrine),true);
    }

    TEST_ASSERT(nonnull(dist.getGlobalComm()));

    TEST_EQUALITY(dist.getApplicationIndex("CTF"),CTF);
    TEST_EQUALITY(dist.getApplicationIndex("Insilico"),Insilico);
    TEST_EQUALITY(dist.getApplicationIndex("Peregrine"),Peregrine);

    TEST_EQUALITY(dist.getApplicationName(CTF),"CTF");
    TEST_EQUALITY(dist.getApplicationName(Insilico),"Insilico");
    TEST_EQUALITY(dist.getApplicationName(Peregrine),"Peregrine");

    TEST_EQUALITY(dist.getTransferIndex("C_TO_I: "),CTF_Insilico);
    TEST_EQUALITY(dist.getTransferIndex("C_TO_P: "),CTF_Peregrine);
    TEST_EQUALITY(dist.getTransferIndex("I_TO_P: "),Insilico_Peregrine);

    TEST_EQUALITY(dist.getTransferName(CTF_Insilico),"C_TO_I: ");
    TEST_EQUALITY(dist.getTransferName(CTF_Peregrine),"C_TO_P: ");
    TEST_EQUALITY(dist.getTransferName(Insilico_Peregrine),"I_TO_P: ");

    

    Teuchos::RCP<Teuchos::FancyOStream> sout = dist.getSerialOStream();
    Teuchos::RCP<Teuchos::FancyOStream> pout = dist.getParallelOStream();

    *sout << "SERIAL OSTREAM TEST!" << std::endl;
    *dist.getApplicationOStream(CTF) << "APPLICATION OSTREAM TEST!" << std::endl;
    *dist.getApplicationOStream(Insilico) << "APPLICATION OSTREAM TEST!" << std::endl;
    *dist.getApplicationOStream(Peregrine) << "APPLICATION OSTREAM TEST!" << std::endl;
    *pout << "PARALLEL OSTREAM TEST!" << std::endl;



    // Test the translate ranks function.
    if (dist.transferExistsOnProcess(Insilico_Peregrine)) {
      int transferRankI = pike::translateMpiRank(dist.getPrintRank(Insilico),
      						*dist.getGlobalComm(),
						 *dist.getTransferComm(Insilico_Peregrine));
      TEST_EQUALITY(transferRankI, 0);
      
      int transferRankP = pike::translateMpiRank(dist.getPrintRank(Peregrine),
      						*dist.getGlobalComm(),
						 *dist.getTransferComm(Insilico_Peregrine));
      TEST_EQUALITY(transferRankP, 3);
    }

  }

  TEUCHOS_UNIT_TEST(MultiphysicsDistributor, overlapping_apps)
  {
    Teuchos::RCP<Teuchos::MpiComm<int> > globalComm = 
      Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));

    // Run this on 6 processes only
    TEST_EQUALITY(globalComm->getSize(), 6);

    typedef pike::MultiphysicsDistributor::ApplicationIndex AppIndex;
    typedef pike::MultiphysicsDistributor::TransferIndex TransIndex;
    AppIndex CTF;
    AppIndex Insilico;
    AppIndex Peregrine;
    TransIndex CTF_Insilico;
    TransIndex CTF_Peregrine;
    TransIndex Insilico_Peregrine;

    pike::MultiphysicsDistributor dist;
    {
      CTF = dist.addApplication("CTF",0,0); // rank range
      Insilico = dist.addApplication("Insilico",1,4); // rank range
      std::vector<int> peregrineRanks;
      peregrineRanks.push_back(1);
      peregrineRanks.push_back(2);
      peregrineRanks.push_back(3);
      peregrineRanks.push_back(4);
      peregrineRanks.push_back(5);
      Peregrine = dist.addApplication("Peregrine",peregrineRanks); // vector of ranks
      
      CTF_Insilico = dist.addTransfer("C_TO_I: ",CTF,Insilico);
      
      std::vector<pike::MultiphysicsDistributor::ApplicationIndex> ctfPeregrine;
      ctfPeregrine.push_back(CTF);
      ctfPeregrine.push_back(Peregrine);
      CTF_Peregrine = dist.addTransfer("C_TO_P: ",ctfPeregrine);
      
      Insilico_Peregrine = dist.addTransfer("I_TO_P: ",Insilico,Peregrine);
      
      dist.setup(globalComm,true);
    }

    if (globalComm->getRank() == 0) {
      TEST_EQUALITY(dist.appExistsOnProcess(CTF),true);
      TEST_EQUALITY(dist.appExistsOnProcess(Insilico),false);
      TEST_EQUALITY(dist.appExistsOnProcess(Peregrine),false);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Insilico),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Peregrine),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(Insilico_Peregrine),false);
    }
    if ( (globalComm->getRank() == 1) || 
	 (globalComm->getRank() == 2) ||  
	 (globalComm->getRank() == 3) ||
	 (globalComm->getRank() == 4) ) {
      TEST_EQUALITY(dist.appExistsOnProcess(CTF),false);
      TEST_EQUALITY(dist.appExistsOnProcess(Insilico),true);
      TEST_EQUALITY(dist.appExistsOnProcess(Peregrine),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Insilico),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Peregrine),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(Insilico_Peregrine),true);
    }
    if ( (globalComm->getRank() == 5) ) {
      TEST_EQUALITY(dist.appExistsOnProcess(CTF),false);
      TEST_EQUALITY(dist.appExistsOnProcess(Insilico),false);
      TEST_EQUALITY(dist.appExistsOnProcess(Peregrine),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Insilico),false);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Peregrine),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(Insilico_Peregrine),true);
    }

    TEST_ASSERT(nonnull(dist.getGlobalComm()));

    Teuchos::RCP<Teuchos::FancyOStream> sout = dist.getSerialOStream();
    Teuchos::RCP<Teuchos::FancyOStream> pout = dist.getParallelOStream();

    *sout << "SERIAL OSTREAM TEST!" << std::endl;
    *dist.getApplicationOStream(CTF) << "APPLICATION OSTREAM TEST!" << std::endl;
    *dist.getApplicationOStream(Insilico) << "APPLICATION OSTREAM TEST!" << std::endl;
    *dist.getApplicationOStream(Peregrine) << "APPLICATION OSTREAM TEST!" << std::endl;
    *pout << "PARALLEL OSTREAM TEST!" << std::endl;

  }
  
  TEUCHOS_UNIT_TEST(MultiphysicsDistributor, corner_cases)
  {
    Teuchos::RCP<Teuchos::MpiComm<int> > globalComm = 
      Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));

    // Run this on 6 processes only
    TEST_EQUALITY(globalComm->getSize(), 6);

    typedef pike::MultiphysicsDistributor::ApplicationIndex AppIndex;
    typedef pike::MultiphysicsDistributor::TransferIndex TransIndex;
    AppIndex CTF;
    AppIndex Insilico;
    AppIndex Peregrine;
    TransIndex CTF_Insilico;
    TransIndex CTF_Peregrine;
    TransIndex Insilico_Peregrine;

    pike::MultiphysicsDistributor dist;
    {
      CTF = dist.addApplication("CTF",0,0); // rank range
      Insilico = dist.addApplication("Insilico",1,4); // rank range
      std::vector<int> peregrineRanks;
      peregrineRanks.push_back(1);
      peregrineRanks.push_back(2);
      peregrineRanks.push_back(3);
      peregrineRanks.push_back(4);
      peregrineRanks.push_back(5);
      peregrineRanks.push_back(5); // test repeated ranks
      Peregrine = dist.addApplication("Peregrine",peregrineRanks);
      
      CTF_Insilico = dist.addTransfer("C_TO_I: ",CTF,Insilico);
      
      std::vector<pike::MultiphysicsDistributor::ApplicationIndex> ctfPeregrine;
      ctfPeregrine.push_back(CTF);
      ctfPeregrine.push_back(Peregrine);
      CTF_Peregrine = dist.addTransfer("C_TO_P: ",ctfPeregrine);
      
      Insilico_Peregrine = dist.addTransfer("I_TO_P: ",Insilico,Peregrine);
      
      dist.setup(globalComm,false);
    }

    if (globalComm->getRank() == 0) {
      TEST_EQUALITY(dist.appExistsOnProcess(CTF),true);
      TEST_EQUALITY(dist.appExistsOnProcess(Insilico),false);
      TEST_EQUALITY(dist.appExistsOnProcess(Peregrine),false);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Insilico),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Peregrine),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(Insilico_Peregrine),false);
    }
    if ( (globalComm->getRank() == 1) || 
	 (globalComm->getRank() == 2) ||  
	 (globalComm->getRank() == 3) ||
	 (globalComm->getRank() == 4) ) {
      TEST_EQUALITY(dist.appExistsOnProcess(CTF),false);
      TEST_EQUALITY(dist.appExistsOnProcess(Insilico),true);
      TEST_EQUALITY(dist.appExistsOnProcess(Peregrine),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Insilico),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Peregrine),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(Insilico_Peregrine),true);
    }
    if ( (globalComm->getRank() == 5) ) {
      TEST_EQUALITY(dist.appExistsOnProcess(CTF),false);
      TEST_EQUALITY(dist.appExistsOnProcess(Insilico),false);
      TEST_EQUALITY(dist.appExistsOnProcess(Peregrine),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Insilico),false);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Peregrine),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(Insilico_Peregrine),true);
    }

    TEST_ASSERT(nonnull(dist.getGlobalComm()));

    Teuchos::RCP<Teuchos::FancyOStream> sout = dist.getSerialOStream();
    Teuchos::RCP<Teuchos::FancyOStream> pout = dist.getParallelOStream();

    *sout << "SERIAL OSTREAM TEST!" << std::endl;
    *dist.getApplicationOStream(CTF) << "APPLICATION OSTREAM TEST!" << std::endl;
    *dist.getApplicationOStream(Insilico) << "APPLICATION OSTREAM TEST!" << std::endl;
    *dist.getApplicationOStream(Peregrine) << "APPLICATION OSTREAM TEST!" << std::endl;
    *pout << "PARALLEL OSTREAM TEST!" << std::endl;

  }

  TEUCHOS_UNIT_TEST(MultiphysicsDistributor, transfer_by_MPI_rank)
  {
    // This tests two parallel coupled codes where the data transfer
    // comm is not the union of the application comms, but is
    // specified by a list of mPI ranks.  This use case occurs if one
    // or more of the codes serializes the data transfer - i.e. it
    // only communicates to other codes via a single process or subset
    // of processes.

    Teuchos::RCP<Teuchos::MpiComm<int> > globalComm = 
      Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));

    // Run this on 6 processes only
    TEST_EQUALITY(globalComm->getSize(), 6);

    typedef pike::MultiphysicsDistributor::ApplicationIndex AppIndex;
    typedef pike::MultiphysicsDistributor::TransferIndex TransIndex;
    AppIndex CTF;
    AppIndex Insilico;
    TransIndex CTF_Insilico;

    pike::MultiphysicsDistributor dist;
    {
      CTF = dist.addApplication("CTF",0,2); // rank range
      Insilico = dist.addApplication("Insilico",3,5); // rank range
      
      // ctf serializes transfers - only occur from process 1.
      std::vector<int> mpiRanks;
      mpiRanks.push_back(1);
      mpiRanks.push_back(3);
      mpiRanks.push_back(4);
      mpiRanks.push_back(5);
      CTF_Insilico = dist.addTransferByRanks("C_TO_I: ",mpiRanks);
      
      dist.setup(globalComm,true);
    }

    if (globalComm->getRank() == 0) {
      TEST_EQUALITY(dist.appExistsOnProcess(CTF),true);
      TEST_EQUALITY(dist.appExistsOnProcess(Insilico),false);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Insilico),false);
    }
    if (globalComm->getRank() == 1) {
      TEST_EQUALITY(dist.appExistsOnProcess(CTF),true);
      TEST_EQUALITY(dist.appExistsOnProcess(Insilico),false);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Insilico),true);
    }
    if (globalComm->getRank() == 2) {
      TEST_EQUALITY(dist.appExistsOnProcess(CTF),true);
      TEST_EQUALITY(dist.appExistsOnProcess(Insilico),false);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Insilico),false);
    }
    if ( (globalComm->getRank() == 3) || 
	 (globalComm->getRank() == 4) ||
	 (globalComm->getRank() == 5) ) {
      TEST_EQUALITY(dist.appExistsOnProcess(CTF),false);
      TEST_EQUALITY(dist.appExistsOnProcess(Insilico),true);
      TEST_EQUALITY(dist.transferExistsOnProcess(CTF_Insilico),true);
    }

  }
  
}
