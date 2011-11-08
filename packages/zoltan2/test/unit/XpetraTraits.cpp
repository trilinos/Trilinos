// @HEADER
// ***********************************************************************
//
//         Zoltan2: Sandia Partitioning Ordering & Coloring Library
//
//                Copyright message goes here.   TODO
//
// Questions? Contact Lee Ann Riesen (lriesen@sandia.gov)
//
// ***********************************************************************
// @HEADER
//
// Basic test of the XpetraTraits definitions.
//
// TODO - a real test would figure out if the migrated objects are
// the same as the original, here we just look at them on stdout.

#include <string>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Comm.hpp>
#include "Teuchos_VerboseObject.hpp"
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <Zoltan2_XpetraTraits.hpp>
#include <UserInputForTests.hpp>

#include <Xpetra_EpetraUtils.hpp>
#ifdef HAVE_MPI
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

using namespace std;
using std::string;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::Array;
using Teuchos::rcp;
using Teuchos::Comm;

template <typename LNO, typename GNO, typename NODE>
  ArrayRCP<GNO> roundRobinMap(
    const RCP<const Tpetra::Map<LNO, GNO, NODE> > &m)
{
  const RCP<const Comm<int> > &comm = m->getComm();
  int proc = comm->getRank();
  int nprocs = comm->getSize();
  GNO base = m->getMinAllGlobalIndex();
  GNO max = m->getMaxAllGlobalIndex();
  size_t globalrows = m->getGlobalNumElements();
  if (globalrows != max - base + 1){
    TEST_FAIL_AND_EXIT(*comm, 0, 
      string("Map is invalid for test - fix test"), 1);
  }
  RCP<Array<GNO> > mygids = rcp(new Array<GNO>);
  GNO firstGNO = proc; 
  if (firstGNO < base){
    GNO n = base % proc;
    if (n>0)
      firstGNO = base - n + proc;
    else
      firstGNO = base;
  }
  for (GNO gid=firstGNO; gid <= max; gid+=nprocs){
    (*mygids).append(gid);
  }

  ArrayRCP<GNO> newIdArcp = Teuchos::arcp(mygids);

  return newIdArcp;
}

typedef int lno_t;
typedef long gno_t;
typedef float scalar_t;
typedef Zoltan2::default_node_t node_t;

typedef int epetra_lno_t;
typedef int epetra_gno_t;
typedef double epetra_scalar_t;

int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv);
  RCP<const Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
  int rank = comm->getRank();

  Teuchos::RCP<Teuchos::FancyOStream> out = 
    Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::EVerbosityLevel v=Teuchos::VERB_EXTREME;

  typedef UserInputForTests<scalar_t, lno_t, gno_t> uinput_t;
  typedef Tpetra::CrsMatrix<scalar_t,lno_t,gno_t,node_t> tmatrix_t;
  typedef Tpetra::Vector<scalar_t,lno_t,gno_t,node_t> tvector_t;
  typedef Tpetra::MultiVector<scalar_t,lno_t,gno_t,node_t> tmvector_t;

  // Create an object that can give us test input.

  RCP<uinput_t> uinput;

  try{
    uinput = rcp(new uinput_t(4, 4, 4, comm));
  }
  catch(std::exception &e){
    TEST_FAIL_AND_EXIT(*comm, 0, string("input ")+e.what(), 1);
  }

    /////////////////////////////////////////////////////////////////
    // XpetraTraits<Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> > 
  {
    RCP<tmatrix_t> M;
  
    try{
      M = uinput->getTpetraCrsMatrix();
    }
    catch(std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("getTpetraCrsMatrix ")+e.what(), 1);
    }
  
    if (rank== 0)
      std::cout << "Original matrix" << std::endl;
  
    Teuchos::RCP<Teuchos::FancyOStream>
      out = Teuchos::VerboseObjectBase::getDefaultOStream();
    Teuchos::EVerbosityLevel v=Teuchos::VERB_EXTREME;
    M->describe(*out,v);
  
    ArrayRCP<gno_t> newRowIds = 
      roundRobinMap<lno_t,gno_t,node_t>(M->getRowMap());
  
    gno_t localNumRows = newRowIds.size();
    gno_t base = M->getIndexBase();
  
    RCP<tmatrix_t> newM = Zoltan2::XpetraTraits<tmatrix_t>::doImport(
      rcp_const_cast<const tmatrix_t>(M),
      localNumRows, newRowIds.getRawPtr(), base);
  
    if (rank== 0)
      std::cout << "Migrated matrix" << std::endl;
  
    newM->describe(*out,v);
  }

    /////////////////////////////////////////////////////////////////
    // XpetraTraits<Tpetra::CrsGraph<scalar_t, lno_t, gno_t, node_t> > 
  {
    RCP<tgraph_t> G;
  
    try{
      G = uinput->getTpetraCrsGraph();
    }
    catch(std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("getTpetraCrsGraph ")+e.what(), 1);
    }
  
    if (rank== 0)
      std::cout << "Original graph" << std::endl;
  
    Teuchos::RCP<Teuchos::FancyOStream>
      out = Teuchos::VerboseObjectBase::getDefaultOStream();
    Teuchos::EVerbosityLevel v=Teuchos::VERB_EXTREME;
    G->describe(*out,v);
  
    ArrayRCP<gno_t> newRowIds = 
      roundRobinMap<lno_t,gno_t,node_t>(G->getRowMap());
  
    gno_t localNumRows = newRowIds.size();
    gno_t base = G->getIndexBase();
  
    RCP<tgraph_t> newG = Zoltan2::XpetraTraits<tgraph_t>::doImport(
      rcp_const_cast<const tgraph_t>(G),
      localNumRows, newRowIds.getRawPtr(), base);
  
    if (rank== 0)
      std::cout << "Migrated graph" << std::endl;
  
    newG->describe(*out,v);
  }

    /////////////////////////////////////////////////////////////////
    // XpetraTraits<Tpetra::Vector<scalar_t, lno_t, gno_t, node_t>> 
  {
    RCP<tvector_t> V;
  
    try{
      V = uinput->getTpetraVector();
    }
    catch(std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("getTpetraVector")+e.what(), 1);
    }
  
    if (rank== 0)
      std::cout << "Original vector" << std::endl;
  
    V->describe(*out,v);
  
    ArrayRCP<gno_t> newRowIds = 
      roundRobinMap<lno_t,gno_t,node_t>(V->getMap());
  
    gno_t localNumRows = newRowIds.size();
    gno_t base = V->getMap()->getIndexBase();
  
    RCP<tvector_t> newV = Zoltan2::XpetraTraits<tvector_t>::doImport(
      rcp_const_cast<const tvector_t>(V),
      localNumRows, newRowIds.getRawPtr(), base);
  
    if (rank== 0)
      std::cout << "Migrated vector" << std::endl;
  
    newV->describe(*out,v);
  }

    /////////////////////////////////////////////////////////////////
    // XpetraTraits<Tpetra::MultiVector<scalar_t, lno_t, gno_t, node_t>> 
  {
    RCP<tmvector_t> MV;
  
    try{
      MV = uinput->getTpetraMultiVector(3);
    }
    catch(std::exception &e){
      TEST_FAIL_AND_EXIT(*comm, 0, 
        string("getTpetraMultiVector")+e.what(), 1);
    }
  
    if (rank== 0)
      std::cout << "Original vector" << std::endl;
  
    MV->describe(*out,v);
  
    ArrayRCP<gno_t> newRowIds = 
      roundRobinMap<lno_t,gno_t,node_t>(MV->getMap());
  
    gno_t localNumRows = newRowIds.size();
    gno_t base = MV->getMap()->getIndexBase();
  
    RCP<tmvector_t> newMV = Zoltan2::XpetraTraits<tmvector_t>::doImport(
      rcp_const_cast<const tmvector_t>(MV),
      localNumRows, newRowIds.getRawPtr(), base);
  
    if (rank== 0)
      std::cout << "Migrated multi vector" << std::endl;
  
    newMV->describe(*out,v);
  }

  if (rank==0)
    std::cout << "PASS" << std::endl;
}

