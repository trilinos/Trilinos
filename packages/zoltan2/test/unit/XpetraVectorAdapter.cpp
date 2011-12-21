// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// ***********************************************************************
//
// Basic testing of Zoltan2::XpetraVectorInput 

#include <string>


#include <Zoltan2_XpetraVectorInput.hpp>
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
typedef Tpetra::Vector<scalar_t, lno_t, gno_t, node_t> tvector_t;
typedef Xpetra::Vector<scalar_t, lno_t, gno_t, node_t> xvector_t;
typedef Epetra_Vector evector_t;

template <typename S, typename L, typename G>
  void printVector(RCP<const Comm<int> > &comm, L vlen,
    const G *vtxIds, const S *vals)
{
  int rank = comm->getRank();
  int nprocs = comm->getSize();
  comm->barrier();
  for (int p=0; p < nprocs; p++){
    if (p == rank){
      std::cout << rank << ":" << std::endl;
      for (L i=0; i < vlen; i++){
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
  Zoltan2::XpetraVectorInput<User> &ia, tvector_t &vector)
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

  gfail = globalFail(comm, fail);

  const G *vtxIds=NULL;
  const S *vals=NULL;
  const S *wgts=NULL;
  size_t nvals=0;

  if (!gfail){

    nvals = ia.getVectorView(vtxIds, vals, wgts);

    if (nvals != vector.getLocalLength())
      fail = 8;
    if (!fail && wgts != NULL)   // not implemented yet
      fail = 10;

    gfail = globalFail(comm, fail);

    if (gfail == 0){
      printVector<S, L, G>(comm, nvals, vtxIds, vals);
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
  // and Epetra vectors for testing.

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

  tV = uinput->getTpetraVector();
  size_t vlen = tV->getLocalLength();
  Teuchos::ArrayView<const gno_t> rowGids = tV->getMap()->getNodeElementList();

  // To test migration in the input adapter we need a Solution
  // object.  The Solution needs an IdentifierMap.

  typedef Zoltan2::IdentifierMap<tvector_t> idmap_t;
  typedef Zoltan2::PartitioningSolution<tvector_t> soln_t;

  RCP<const Zoltan2::Environment> env = Zoltan2::getDefaultEnvironment();

  ArrayRCP<const gno_t> gidArray = arcpFromArrayView(rowGids);
  RCP<const idmap_t> idMap = rcp(new idmap_t(env, gidArray));

  int weightDim = 1;
  float *imbal = new float [weightDim];
  imbal[0] = 1.0;
  ArrayRCP<float> metric(imbal, 0, 1, true);

  size_t *p = new size_t [vlen];
  memset(p, 0, sizeof(size_t) * vlen);
  ArrayRCP<size_t> solnParts(p, 0, vlen, true);

  soln_t solution(env, idMap, weightDim);

  solution.setParts(rowGids, solnParts, metric);

  /////////////////////////////////////////////////////////////
  // User object is Tpetra::Vector
  if (!gfail){ 
    RCP<const tvector_t> ctV = rcp_const_cast<const tvector_t>(tV);
    RCP<Zoltan2::XpetraVectorInput<tvector_t> > tVInput;
  
    try {
      tVInput = 
        rcp(new Zoltan2::XpetraVectorInput<tvector_t>(ctV));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("XpetraVectorInput ")+e.what(), 1);
    }
  
    if (rank==0)
      std::cout << "Input adapter for Tpetra::Vector" << std::endl;
    
    fail = verifyInputAdapter<tvector_t>(*tVInput, *tV);
  
    gfail = globalFail(comm, fail);
  
    if (!gfail){
      tvector_t *vMigrate = NULL;
      try{
        tVInput->applyPartitioningSolution<tvector_t>(*tV, vMigrate, solution);
        newV = rcp(vMigrate);
      }
      catch (std::exception &e){
        fail = 11;
      }

      gfail = globalFail(comm, fail);
  
      if (!gfail){
        RCP<const tvector_t> cnewV = rcp_const_cast<const tvector_t>(newV);
        RCP<Zoltan2::XpetraVectorInput<tvector_t> > newInput;
        try{
          newInput = rcp(new Zoltan2::XpetraVectorInput<tvector_t>(cnewV));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0, 
            string("XpetraVectorInput 2 ")+e.what(), 1);
        }
  
        if (rank==0){
          std::cout << 
           "Input adapter for Tpetra::Vector migrated to proc 0" << 
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
  // User object is Xpetra::Vector
  if (!gfail){ 
    RCP<xvector_t> xV = uinput->getXpetraVector();
    RCP<const xvector_t> cxV = rcp_const_cast<const xvector_t>(xV);
    RCP<Zoltan2::XpetraVectorInput<xvector_t> > xVInput;
  
    try {
      xVInput = 
        rcp(new Zoltan2::XpetraVectorInput<xvector_t>(cxV));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("XpetraVectorInput 3 ")+e.what(), 1);
    }
  
    if (rank==0){
      std::cout << "Input adapter for Xpetra::Vector" << std::endl;
    }
    fail = verifyInputAdapter<xvector_t>(*xVInput, *tV);
  
    gfail = globalFail(comm, fail);
  
    if (!gfail){
      xvector_t *vMigrate =NULL;
      try{
        xVInput->applyPartitioningSolution<tvector_t>(*xV, vMigrate, solution);
      }
      catch (std::exception &e){
        fail = 11;
      }
  
      gfail = globalFail(comm, fail);
  
      if (!gfail){
        RCP<const xvector_t> cnewV(vMigrate);
        RCP<Zoltan2::XpetraVectorInput<xvector_t> > newInput;
        try{
          newInput = 
            rcp(new Zoltan2::XpetraVectorInput<xvector_t>(cnewV));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0, 
            string("XpetraVectorInput 4 ")+e.what(), 1);
        }
  
        if (rank==0){
          std::cout << 
           "Input adapter for Xpetra::Vector migrated to proc 0" << 
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
  // User object is Epetra_Vector
  if (!gfail){ 
    RCP<evector_t> eV = uinput->getEpetraVector();
    RCP<const evector_t> ceV = rcp_const_cast<const evector_t>(eV);
    RCP<Zoltan2::XpetraVectorInput<evector_t> > eVInput;
  
    try {
      eVInput = 
        rcp(new Zoltan2::XpetraVectorInput<evector_t>(ceV));
    }
    catch (std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("XpetraVectorInput 5 ")+e.what(), 1);
    }
  
    if (rank==0){
      std::cout << "Input adapter for Epetra_Vector" << std::endl;
    }
    fail = verifyInputAdapter<evector_t>(*eVInput, *tV);
  
    gfail = globalFail(comm, fail);
  
    if (!gfail){
      evector_t *vMigrate =NULL;
      try{
        eVInput->applyPartitioningSolution<tvector_t>(*eV, vMigrate, solution);
      }
      catch (std::exception &e){
        fail = 11;
      }
  
      gfail = globalFail(comm, fail);
  
      if (!gfail){
        RCP<const evector_t> cnewV(vMigrate, true);
        RCP<Zoltan2::XpetraVectorInput<evector_t> > newInput;
        try{
          newInput = 
            rcp(new Zoltan2::XpetraVectorInput<evector_t>(cnewV));
        }
        catch (std::exception &e){
          TEST_FAIL_AND_EXIT(*comm, 0, 
            string("XpetraVectorInput 6 ")+e.what(), 1);
        }
  
        if (rank==0){
          std::cout << 
           "Input adapter for Epetra_Vector migrated to proc 0" << 
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
