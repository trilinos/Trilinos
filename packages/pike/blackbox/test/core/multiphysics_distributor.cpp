#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "Pike_MultiphysicsDistributor.hpp"

namespace pike {

  // Applications
  const int CTF = 0;
  const int Insilico = 1;
  const int Peregrine = 2;

  //Transfers
  const int CTF_Insilico = 0;
  const int CTF_Peregrine = 1;
  const int Insilico_Peregrine = 2;

  TEUCHOS_UNIT_TEST(MultiphysicsDistributor, num_cores_version)
  {
    using pike::CTF;
    using pike::Insilico;
    using pike::Peregrine;
    using pike::CTF_Insilico;
    using pike::CTF_Peregrine;
    using pike::Insilico_Peregrine;

    Teuchos::RCP<Teuchos::MpiComm<int> > globalComm = 
      Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));

    // Run this on 6 processes only
    TEST_EQUALITY(globalComm->getSize(), 6);

    pike::MultiphysicsDistributor dist;
    {
      // Applciations are added via num cores
      dist.addApplication(CTF,"CTF",1);
      dist.addApplication(Insilico,"Insilico",3);
      dist.addApplication(Peregrine,"Peregrine",2);
      
      // try to add application already registered (check for
      // duplicate index and duplicate name)
      TEST_THROW(dist.addApplication(CTF,"CTF DUPLICATE",1), std::logic_error);
      TEST_THROW(dist.addApplication(4,"CTF",1), std::logic_error);

      dist.addTransfer(CTF_Insilico,CTF,Insilico,"C_TO_I: ");
      
      std::vector<pike::MultiphysicsDistributor::ApplicationIndex> ctfPeregrine;
      ctfPeregrine.push_back(CTF);
      ctfPeregrine.push_back(Peregrine);
      dist.addTransfer(CTF_Peregrine,ctfPeregrine,"C_TO_P: ");
      
      dist.addTransfer(Insilico_Peregrine,Insilico,Peregrine,"I_TO_P: ");
      
      // try to add transfer already registered (check for duplicate
      // index and duplicate name)
      TEST_THROW(dist.addTransfer(CTF_Insilico,CTF,Insilico,"C_TO_I: DUPLICATE"), std::logic_error);
      TEST_THROW(dist.addTransfer(4,CTF,Insilico,"C_TO_I: "), std::logic_error);
      TEST_THROW(dist.addTransfer(CTF_Peregrine,ctfPeregrine,"C_TO_P: DUPLICATE"), std::logic_error);
      TEST_THROW(dist.addTransfer(4,ctfPeregrine,"C_TO_P: "), std::logic_error);

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

  }

  TEUCHOS_UNIT_TEST(MultiphysicsDistributor, range_version)
  {
    using pike::CTF;
    using pike::Insilico;
    using pike::Peregrine;
    using pike::CTF_Insilico;
    using pike::CTF_Peregrine;
    using pike::Insilico_Peregrine;

    Teuchos::RCP<Teuchos::MpiComm<int> > globalComm = 
      Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));

    // Run this on 6 processes only
    TEST_EQUALITY(globalComm->getSize(), 6);

    pike::MultiphysicsDistributor dist;
    {
      // Applciations are added via mixed methods (rank range or vector of ranks
      dist.addApplication(CTF,"CTF",0,0); // rank range
      dist.addApplication(Insilico,"Insilico",1,3); // rank range
      std::vector<int> peregrineRanks;
      peregrineRanks.push_back(4);
      peregrineRanks.push_back(5);
      dist.addApplication(Peregrine,"Peregrine",peregrineRanks); // vector of ranks
      
      // try to add application already registered
      TEST_THROW(dist.addApplication(CTF,"CTF",1), std::logic_error);

      dist.addTransfer(CTF_Insilico,CTF,Insilico,"C_TO_I: ");
      
      std::vector<pike::MultiphysicsDistributor::ApplicationIndex> ctfPeregrine;
      ctfPeregrine.push_back(CTF);
      ctfPeregrine.push_back(Peregrine);
      dist.addTransfer(CTF_Peregrine,ctfPeregrine,"C_TO_P: ");
      
      dist.addTransfer(Insilico_Peregrine,Insilico,Peregrine,"I_TO_P: ");
      
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

    Teuchos::RCP<Teuchos::FancyOStream> sout = dist.getSerialOStream();
    Teuchos::RCP<Teuchos::FancyOStream> pout = dist.getParallelOStream();

    *sout << "SERIAL OSTREAM TEST!" << std::endl;
    *dist.getApplicationOStream(CTF) << "APPLICATION OSTREAM TEST!" << std::endl;
    *dist.getApplicationOStream(Insilico) << "APPLICATION OSTREAM TEST!" << std::endl;
    *dist.getApplicationOStream(Peregrine) << "APPLICATION OSTREAM TEST!" << std::endl;
    *pout << "PARALLEL OSTREAM TEST!" << std::endl;

  }

  TEUCHOS_UNIT_TEST(MultiphysicsDistributor, overlapping_version)
  {
    using pike::CTF;
    using pike::Insilico;
    using pike::Peregrine;
    using pike::CTF_Insilico;
    using pike::CTF_Peregrine;
    using pike::Insilico_Peregrine;

    Teuchos::RCP<Teuchos::MpiComm<int> > globalComm = 
      Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));

    // Run this on 6 processes only
    TEST_EQUALITY(globalComm->getSize(), 6);

    pike::MultiphysicsDistributor dist;
    {
      // Applciations are added via mixed methods
      dist.addApplication(CTF,"CTF",0,0); // rank range
      dist.addApplication(Insilico,"Insilico",1,4); // rank range
      std::vector<int> peregrineRanks;
      peregrineRanks.push_back(1);
      peregrineRanks.push_back(2);
      peregrineRanks.push_back(3);
      peregrineRanks.push_back(4);
      peregrineRanks.push_back(5);
      dist.addApplication(Peregrine,"Peregrine",peregrineRanks); // vector of ranks
      
      dist.addTransfer(CTF_Insilico,CTF,Insilico,"C_TO_I: ");
      
      std::vector<pike::MultiphysicsDistributor::ApplicationIndex> ctfPeregrine;
      ctfPeregrine.push_back(CTF);
      ctfPeregrine.push_back(Peregrine);
      dist.addTransfer(CTF_Peregrine,ctfPeregrine,"C_TO_P: ");
      
      dist.addTransfer(Insilico_Peregrine,Insilico,Peregrine,"I_TO_P: ");
      
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
    using pike::CTF;
    using pike::Insilico;
    using pike::Peregrine;
    using pike::CTF_Insilico;
    using pike::CTF_Peregrine;
    using pike::Insilico_Peregrine;

    Teuchos::RCP<Teuchos::MpiComm<int> > globalComm = 
      Teuchos::rcp(new Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(MPI_COMM_WORLD)));

    // Run this on 6 processes only
    TEST_EQUALITY(globalComm->getSize(), 6);

    pike::MultiphysicsDistributor dist;
    {
      // Applciations are added via num cores
      dist.addApplication(CTF,"CTF",0,0); // rank range
      dist.addApplication(Insilico,"Insilico",1,4); // rank range
      std::vector<int> peregrineRanks;
      peregrineRanks.push_back(1);
      peregrineRanks.push_back(2);
      peregrineRanks.push_back(3);
      peregrineRanks.push_back(4);
      peregrineRanks.push_back(5);
      peregrineRanks.push_back(5); // test repeated ranks
      dist.addApplication(Peregrine,"Peregrine",peregrineRanks);
      
      dist.addTransfer(CTF_Insilico,CTF,Insilico,"C_TO_I: ");
      
      std::vector<pike::MultiphysicsDistributor::ApplicationIndex> ctfPeregrine;
      ctfPeregrine.push_back(CTF);
      ctfPeregrine.push_back(Peregrine);
      dist.addTransfer(CTF_Peregrine,ctfPeregrine,"C_TO_P: ");
      
      dist.addTransfer(Insilico_Peregrine,Insilico,Peregrine,"I_TO_P: ");
      
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

}
