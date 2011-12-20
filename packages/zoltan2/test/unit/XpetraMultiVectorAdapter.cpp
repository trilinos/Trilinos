// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Basic testing of Zoltan2::XpetraMultiVectorInput 

#include <string>


#include <Zoltan2_XpetraMultiVectorInput.hpp>
#include <Zoltan2_InputTraits.hpp>
#include <UserInputForTests.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_CommHelpers.hpp>

using namespace std;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_const_cast;
using Teuchos::Comm;
using Teuchos::DefaultComm;

typedef double scalar_t;
typedef int lno_t;
typedef int gno_t;
typedef Zoltan2::default_node_t node_t;

typedef UserInputForTests<scalar_t, lno_t, gno_t> uinput_t;
typedef Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> tvector_t;
typedef Xpetra::MultiVector<scalar_t, lno_t, gno_t, node_t> xvector_t;
typedef Epetra_MultiVector evector_t;

template <typename S, typename L, typename G>
  void printMultiVector(int i, RCP<const Comm<int> > &comm, L nvtx,
    const G *vtxIds, const S *vals)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  comm->barrier();
  if (rank == 0)
    std::cout << "Vector " << i << std::endl;
  for (int p=0; p < nprocs; p++){
    if (p == rank){
      std::cout << rank << ":" << std::endl;
      for (L i=0; i < nvtx; i++){
        std::cout << " " << vtxIds[i] << ": " << vals[i] << std::endl;
      }
      std::cout.flush();
    }
    comm->barrier();
  }
  comm->barrier();
}

template <typename User>
int verifyInputAdapter(
  Zoltan2::XpetraMultiVectorInput<User> &ia, tvector_t &vector)
{
  typedef typename Zoltan2::InputTraits<User>::scalar_t S;
  typedef typename Zoltan2::InputTraits<User>::lno_t L;
  typedef typename Zoltan2::InputTraits<User>::gno_t G;

  RCP<const Comm<int> > comm = vector.getMap()->getComm();
  int fail = 0, gfail=0;

  if (!fail && ia.getLocalLength() != vector.getLocalLength())
    fail = 4;

  if (!fail && ia.getGlobalLength() != vector.getGlobalLength())
    fail = 5;

  if (!fail && ia.getNumVectors() != vector.getNumVectors())
    fail = 5;

  gfail = globalFail(comm, fail);

  const G *vtxIds=NULL;
  const S *vals=NULL;
  const S *wgts=NULL;
  size_t nvals=0;

  if (!gfail){

    for (int i=0; i < ia.getNumVectors(); i++){
      nvals = ia.getMultiVectorView(i, vtxIds, vals, wgts);
  
      if (nvals != vector.getLocalLength())
        fail = i*100 + 8;
      if (!fail && wgts != NULL)   // not implemented yet
        fail = i*100 + 10;
  
      gfail = globalFail(comm, fail);
  
      if (gfail == 0){
        printMultiVector<S, L, G>(i, comm, nvals, vtxIds, vals);
      }
    }
  }
  return fail;
}

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = DefaultComm<int>::getComm();
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  int fail = 0, gfail=0;

  // Create object that can give us test Tpetra, Xpetra
  // and Epetra multivectors for testing.

  RCP<uinput_t> uinput;

  try{
    uinput = 
      rcp(new uinput_t(std::string("../data/simple.mtx"), comm));
  }
  catch(std::exception &e){
    TEST_FAIL_AND_EXIT(*comm, 0, string("input ")+e.what(), 1);
  }

  RCP<tvector_t> tV;     // original vector (for checking)
  RCP<tvector_t> newV;   // migrated vector

  tV = uinput->getTpetraMultiVector(3);
  size_t vlen = tV->getLocalLength();
  Teuchos::ArrayView<const gno_t> rowGids = tV->getMap()->getNodeElementList();

  // To test migration in the input adapter we need a Solution
  // object.  The Solution needs an IdentifierMap.

  Zoltan2::BasicUserTypes<scalar_t, gno_t, lno_t, gno_t> UserTypes;

  typedef Zoltan2::IdentifierMap<UserTypes> idmap_t;
  typedef Zoltan2::PartitioningSolution<UserTypes> soln_t;

  RCP<const Zoltan2::Environment> env = Zoltan2::getDefaultEnvironment();

  ArrayRCP<const gno_t> gidArray = arcpFromArrayView(rowGids);
  RCP<const idmap_t> idMap = rcp(new idmap_t(env, gidArray));

  int weightDim = 1;
  scalar_t *imbal = new scalar_t [weightDim];
  imbal[0] = 1.0;
  ArrayRCP<scalar_t> metric(imbal, 0, 1, true);

  size_t *p = new size_t [nvtx];
  memset(p, 0, sizeof(size_t) * nvtx);
  ArrayRCP<size_t> solnParts(p, 0, nvtx, true);

  soln_t solution(env, idMap, weightDim);

  solution.setParts(rowGids, solnParts, metric);

  /////////////////////////////////////////////////////////////
  // User object is Tpetra::MultiVector
  if (!gfail){ 
    RCP<const tvector_t> ctV = rcp_const_cast<const tvector_t>(tV);
    RCP<Zoltan2::XpetraMultiVectorInput<tvector_t> > tVInput;
  
    try {
      tVInput = 
        rcp(new Zoltan2::XpetraMultiVectorInput<tvector_t>(ctV));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("XpetraMultiVectorInput ")+e.what(), 1);
    }
  
    if (rank==0)
      std::cout << "Input adapter for Tpetra::MultiVector" << std::endl;
    
    fail = verifyInputAdapter<tvector_t>(*tVInput, *tV);
  
    gfail = globalFail(comm, fail);
  
    if (!gfail){
      tvector_t *vMigrate = NULL;
      try{
        tVInput->applyPartitioningSolution(*tV, vMigrate, solution);
        newV = rcp(vMigrate);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(comm, fail);
  
      if (!gfail){
        RCP<const tvector_t> cnewV = rcp_const_cast<const tvector_t>(newV);
        RCP<Zoltan2::XpetraMultiVectorInput<tvector_t> > newInput;
        try{
          newInput = rcp(new Zoltan2::XpetraMultiVectorInput<tvector_t>(cnewV));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0, 
            string("XpetraMultiVectorInput 2 ")+e.what(), 1);
        }
  
        if (rank==0){
          std::cout << 
           "Input adapter for Tpetra::MultiVector migrated to proc 0" << 
           std::endl;
        }
        fail = verifyInputAdapter<tvector_t>(*newInput, *newV);
        if (fail) fail += 100;
        gfail = globalFail(comm, fail);
      }
    }
    if (gfail){
      printFailureCode(comm, fail);
    }
  }

  /////////////////////////////////////////////////////////////
  // User object is Xpetra::MultiVector
  if (!gfail){ 
    RCP<xvector_t> xV = uinput->getXpetraMultiVector(3);
    RCP<const xvector_t> cxV = rcp_const_cast<const xvector_t>(xV);
    RCP<Zoltan2::XpetraMultiVectorInput<xvector_t> > xVInput;
  
    try {
      xVInput = 
        rcp(new Zoltan2::XpetraMultiVectorInput<xvector_t>(cxV));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("XpetraMultiVectorInput 3 ")+e.what(), 1);
    }
  
    if (rank==0){
      std::cout << "Input adapter for Xpetra::MultiVector" << std::endl;
    }
    fail = verifyInputAdapter<xvector_t>(*xVInput, *tV);
  
    gfail = globalFail(comm, fail);
  
    if (!gfail){
      xvector_t *vMigrate =NULL;
      try{
        xVInput->applyPartitioningSolution(*xV, vMigrate, solution);
      }
      catch (std::exception &e){
        fail = 11;
      }
  
      gfail = globalFail(comm, fail);
  
      if (!gfail){
        RCP<const xvector_t> cnewV(vMigrate);
        RCP<Zoltan2::XpetraMultiVectorInput<xvector_t> > newInput;
        try{
          newInput = 
            rcp(new Zoltan2::XpetraMultiVectorInput<xvector_t>(cnewV));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0, 
            string("XpetraMultiVectorInput 4 ")+e.what(), 1);
        }
  
        if (rank==0){
          std::cout << 
           "Input adapter for Xpetra::MultiVector migrated to proc 0" << 
           std::endl;
        }
        fail = verifyInputAdapter<xvector_t>(*newInput, *newV);
        if (fail) fail += 100;
        gfail = globalFail(comm, fail);
      }
    }
    if (gfail){
      printFailureCode(comm, fail);
    }
  }

  /////////////////////////////////////////////////////////////
  // User object is Epetra_MultiVector
  if (!gfail){ 
    RCP<evector_t> eV = uinput->getEpetraMultiVector(3);
    RCP<const evector_t> ceV = rcp_const_cast<const evector_t>(eV);
    RCP<Zoltan2::XpetraMultiVectorInput<evector_t> > eVInput;
  
    try {
      eVInput = 
        rcp(new Zoltan2::XpetraMultiVectorInput<evector_t>(ceV));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("XpetraMultiVectorInput 5 ")+e.what(), 1);
    }
  
    if (rank==0){
      std::cout << "Input adapter for Epetra_MultiVector" << std::endl;
    }
    fail = verifyInputAdapter<evector_t>(*eVInput, *tV);
  
    gfail = globalFail(comm, fail);
  
    if (!gfail){
      evector_t *vMigrate =NULL;
      try{
        eVInput->applyPartitioningSolution(*eV, vMigrate, solution);
      }
      catch (std::exception &e){
        fail = 11;
      }
  
      gfail = globalFail(comm, fail);
  
      if (!gfail){
        RCP<const evector_t> cnewV(vMigrate, true);
        RCP<Zoltan2::XpetraMultiVectorInput<evector_t> > newInput;
        try{
          newInput = 
            rcp(new Zoltan2::XpetraMultiVectorInput<evector_t>(cnewV));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0, 
            string("XpetraMultiVectorInput 6 ")+e.what(), 1);
        }
  
        if (rank==0){
          std::cout << 
           "Input adapter for Epetra_MultiVector migrated to proc 0" << 
           std::endl;
        }
        fail = verifyInputAdapter<evector_t>(*newInput, *newV);
        if (fail) fail += 100;
        gfail = globalFail(comm, fail);
      }
    }
    if (gfail){
      printFailureCode(comm, fail);
    }
  }

  /////////////////////////////////////////////////////////////
  // DONE

  if (rank==0)
    std::cout << "PASS" << std::endl;
}
