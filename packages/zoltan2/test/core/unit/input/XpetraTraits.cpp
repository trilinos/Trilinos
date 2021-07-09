// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
//
// Basic test of the XpetraTraits definitions.
//
// TODO - a real test would figure out if the migrated objects are
// the same as the original, here we just look at them on stdout.
// TODO look at number of diagonals and max number of entries in
//   Tpetra and Xpetra migrated graphs.  They're garbage.

#include <Zoltan2_XpetraTraits.hpp>
#include <Zoltan2_TestHelpers.hpp>

#include <string>
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

#ifdef HAVE_ZOLTAN_EPETRA
#include <Xpetra_EpetraUtils.hpp>
#endif

using std::string;
using Teuchos::RCP;
using Teuchos::ArrayRCP;
using Teuchos::ArrayView;
using Teuchos::Array;
using Teuchos::rcp;
using Teuchos::Comm;

ArrayRCP<zgno_t> roundRobinMapShared(
  int proc,
  int nprocs,
  zgno_t basegid,
  zgno_t maxgid,
  size_t nglobalrows
)
{
  if (nglobalrows != size_t(maxgid - basegid + 1)){
    std::cout << "Error:  Map is invalid for test - fix test" << std::endl;
    std::cerr << "Error:  Map is invalid for test - fix test" << std::endl;
    std::cout << "FAIL" << std::endl;
    exit(1);
  }
  RCP<Array<zgno_t> > mygids = rcp(new Array<zgno_t>);
  zgno_t firstzgno_t = proc;
  if (firstzgno_t < basegid){
    zgno_t n = basegid % proc;
    if (n>0)
      firstzgno_t = basegid - n + proc;
    else
      firstzgno_t = basegid;
  }
  for (zgno_t gid=firstzgno_t; gid <= maxgid; gid+=nprocs){
    (*mygids).append(gid);
  }

  ArrayRCP<zgno_t> newIdArcp = Teuchos::arcp(mygids);

  return newIdArcp;
}

#ifdef HAVE_EPETRA_DATA_TYPES
ArrayRCP<zgno_t> roundRobinMap(const Epetra_BlockMap &emap)
{
  const Epetra_Comm &comm = emap.Comm();
  int proc = comm.MyPID();
  int nprocs = comm.NumProc();
  zgno_t basegid = emap.MinAllGID();
  zgno_t maxgid = emap.MaxAllGID();
  size_t nglobalrows = emap.NumGlobalElements();

  return roundRobinMapShared(proc, nprocs, basegid, maxgid, nglobalrows);
}
#endif

ArrayRCP<zgno_t> roundRobinMap(const Tpetra::Map<zlno_t, zgno_t, znode_t> &tmap)
{
  const RCP<const Comm<int> > &comm = tmap.getComm();
  int proc = comm->getRank();
  int nprocs = comm->getSize();
  zgno_t basegid = tmap.getMinAllGlobalIndex();
  zgno_t maxgid = tmap.getMaxAllGlobalIndex();
  size_t nglobalrows = tmap.getGlobalNumElements();

  return roundRobinMapShared(proc, nprocs, basegid, maxgid, nglobalrows);
}

ArrayRCP<zgno_t> roundRobinMap(const Xpetra::Map<zlno_t, zgno_t, znode_t> &xmap)
{
  const RCP<const Comm<int> > &comm = xmap.getComm();
  int proc = comm->getRank();
  int nprocs = comm->getSize();
  zgno_t basegid = xmap.getMinAllGlobalIndex();
  zgno_t maxgid = xmap.getMaxAllGlobalIndex();
  size_t nglobalrows = xmap.getGlobalNumElements();

  return roundRobinMapShared(proc, nprocs, basegid, maxgid, nglobalrows);
}

int main(int narg, char *arg[])
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  int rank = comm->getRank();
  bool aok = true;

  Teuchos::RCP<Teuchos::FancyOStream> outStream =
    Teuchos::VerboseObjectBase::getDefaultOStream();
  Teuchos::EVerbosityLevel v=Teuchos::VERB_EXTREME;

  typedef Tpetra::CrsMatrix<zscalar_t,zlno_t,zgno_t,znode_t> tmatrix_t;
  typedef Tpetra::CrsGraph<zlno_t,zgno_t,znode_t> tgraph_t;
  typedef Tpetra::Vector<zscalar_t,zlno_t,zgno_t,znode_t> tvector_t;
  typedef Tpetra::MultiVector<zscalar_t,zlno_t,zgno_t,znode_t> tmvector_t;
  typedef Xpetra::CrsMatrix<zscalar_t,zlno_t,zgno_t,znode_t> xmatrix_t;
  typedef Xpetra::CrsGraph<zlno_t,zgno_t,znode_t> xgraph_t;
  typedef Xpetra::Vector<zscalar_t,zlno_t,zgno_t,znode_t> xvector_t;
  typedef Xpetra::MultiVector<zscalar_t,zlno_t,zgno_t,znode_t> xmvector_t;

  // Create object that can give us test Tpetra and Xpetra input.

  RCP<UserInputForTests> uinput;
  Teuchos::ParameterList params;
  params.set("input file", "simple");
  params.set("file type", "Chaco");

  try{
    uinput = rcp(new UserInputForTests(params, comm));
  }
  catch(std::exception &e){
    aok = false;
    std::cout << e.what() << std::endl;
  }
  TEST_FAIL_AND_EXIT(*comm, aok, "input ", 1);

  /////////////////////////////////////////////////////////////////
  //   Tpetra::CrsMatrix
  //   Tpetra::CrsGraph
  //   Tpetra::Vector
  //   Tpetra::MultiVector
  /////////////////////////////////////////////////////////////////


  // XpetraTraits<Tpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t, znode_t> >
  {
    RCP<tmatrix_t> M;

    try{
      M = uinput->getUITpetraCrsMatrix();
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "getTpetraCrsMatrix ", 1);

    if (rank== 0)
      std::cout << "Original Tpetra matrix " << M->getGlobalNumRows()
        << " x " << M->getGlobalNumCols() << std::endl;

    M->describe(*outStream,v);

    ArrayRCP<zgno_t> newRowIds = roundRobinMap(*(M->getRowMap()));

    zgno_t localNumRows = newRowIds.size();

    RCP<const tmatrix_t> newM;
    try{
      newM = Zoltan2::XpetraTraits<tmatrix_t>::doMigration(*M,
        localNumRows, newRowIds.getRawPtr());
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok,
        " Zoltan2::XpetraTraits<tmatrix_t>::doMigration ", 1);

    if (rank== 0)
      std::cout << "Migrated Tpetra matrix" << std::endl;

    newM->describe(*outStream,v);
  }

  // XpetraTraits<Tpetra::CrsGraph<zscalar_t, zlno_t, zgno_t, znode_t> >
  {
    RCP<tgraph_t> G;

    try{
      G = uinput->getUITpetraCrsGraph();
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "getTpetraCrsGraph ", 1);

    if (rank== 0)
      std::cout << "Original Tpetra graph" << std::endl;

    G->describe(*outStream,v);

    ArrayRCP<zgno_t> newRowIds = roundRobinMap(*(G->getRowMap()));

    zgno_t localNumRows = newRowIds.size();

    RCP<const tgraph_t> newG;
    try{
      newG = Zoltan2::XpetraTraits<tgraph_t>::doMigration(*G,
        localNumRows, newRowIds.getRawPtr());
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok,
        " Zoltan2::XpetraTraits<tgraph_t>::doMigration ", 1);

    if (rank== 0)
      std::cout << "Migrated Tpetra graph" << std::endl;

    newG->describe(*outStream,v);
  }

  // XpetraTraits<Tpetra::Vector<zscalar_t, zlno_t, zgno_t, znode_t>>
  {
    RCP<tvector_t> V;

    try{
      V = rcp(new tvector_t(uinput->getUITpetraCrsGraph()->getRowMap(),  1));
      V->randomize();
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "getTpetraVector", 1);

    if (rank== 0)
      std::cout << "Original Tpetra vector" << std::endl;

    V->describe(*outStream,v);

    ArrayRCP<zgno_t> newRowIds = roundRobinMap(*(V->getMap()));

    zgno_t localNumRows = newRowIds.size();

    RCP<const tvector_t> newV;
    try{
      newV = Zoltan2::XpetraTraits<tvector_t>::doMigration(*V,
        localNumRows, newRowIds.getRawPtr());
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok,
        " Zoltan2::XpetraTraits<tvector_t>::doMigration ", 1);

    if (rank== 0)
      std::cout << "Migrated Tpetra vector" << std::endl;

    newV->describe(*outStream,v);
  }

  // XpetraTraits<Tpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t>>
  {
    RCP<tmvector_t> MV;

    try{
      MV = rcp(new tmvector_t(uinput->getUITpetraCrsGraph()->getRowMap(), 3));
      MV->randomize();
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "getTpetraMultiVector", 1);

    if (rank== 0)
      std::cout << "Original Tpetra multivector" << std::endl;

    MV->describe(*outStream,v);

    ArrayRCP<zgno_t> newRowIds = roundRobinMap(*(MV->getMap()));

    zgno_t localNumRows = newRowIds.size();

    RCP<const tmvector_t> newMV;
    try{
      newMV = Zoltan2::XpetraTraits<tmvector_t>::doMigration(*MV,
        localNumRows, newRowIds.getRawPtr());
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok,
        " Zoltan2::XpetraTraits<tmvector_t>::doMigration ", 1);

    if (rank== 0)
      std::cout << "Migrated Tpetra multivector" << std::endl;

    newMV->describe(*outStream,v);
  }

  /////////////////////////////////////////////////////////////////
  //   Xpetra::CrsMatrix
  //   Xpetra::CrsGraph
  //   Xpetra::Vector
  //   Xpetra::MultiVector
  /////////////////////////////////////////////////////////////////

  // XpetraTraits<Xpetra::CrsMatrix<zscalar_t, zlno_t, zgno_t, znode_t> >
  {
    RCP<xmatrix_t> M;

    try{
      M = uinput->getUIXpetraCrsMatrix();
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "getXpetraCrsMatrix ", 1);

    if (rank== 0)
      std::cout << "Original Xpetra matrix" << std::endl;

    M->describe(*outStream,v);

    ArrayRCP<zgno_t> newRowIds = roundRobinMap(*(M->getRowMap()));

    zgno_t localNumRows = newRowIds.size();

    RCP<const xmatrix_t> newM;
    try{
      newM = Zoltan2::XpetraTraits<xmatrix_t>::doMigration(*M,
        localNumRows, newRowIds.getRawPtr());
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok,
        " Zoltan2::XpetraTraits<xmatrix_t>::doMigration ", 1);

    if (rank== 0)
      std::cout << "Migrated Xpetra matrix" << std::endl;

    newM->describe(*outStream,v);
  }

  // XpetraTraits<Xpetra::CrsGraph<zscalar_t, zlno_t, zgno_t, znode_t> >
  {
    RCP<xgraph_t> G;

    try{
      G = uinput->getUIXpetraCrsGraph();
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "getXpetraCrsGraph ", 1);

    if (rank== 0)
      std::cout << "Original Xpetra graph" << std::endl;

    G->describe(*outStream,v);

    ArrayRCP<zgno_t> newRowIds = roundRobinMap(*(G->getRowMap()));

    zgno_t localNumRows = newRowIds.size();

    RCP<const xgraph_t> newG;
    try{
      newG = Zoltan2::XpetraTraits<xgraph_t>::doMigration(*G,
        localNumRows, newRowIds.getRawPtr());
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok,
        " Zoltan2::XpetraTraits<xgraph_t>::doMigration ", 1);

    if (rank== 0)
      std::cout << "Migrated Xpetra graph" << std::endl;

    newG->describe(*outStream,v);
  }

  // XpetraTraits<Xpetra::Vector<zscalar_t, zlno_t, zgno_t, znode_t>>
  {
    RCP<xvector_t> V;

    try{
      RCP<tvector_t> tV =
          rcp(new tvector_t(uinput->getUITpetraCrsGraph()->getRowMap(),  1));
      tV->randomize();
      V = Zoltan2::XpetraTraits<tvector_t>::convertToXpetra(tV);
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "getXpetraVector", 1);

    if (rank== 0)
      std::cout << "Original Xpetra vector" << std::endl;

    V->describe(*outStream,v);

    ArrayRCP<zgno_t> newRowIds = roundRobinMap(*(V->getMap()));

    zgno_t localNumRows = newRowIds.size();

    RCP<const xvector_t> newV;
    try{
      newV = Zoltan2::XpetraTraits<xvector_t>::doMigration(*V,
        localNumRows, newRowIds.getRawPtr());
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok,
        " Zoltan2::XpetraTraits<xvector_t>::doMigration ", 1);

    if (rank== 0)
      std::cout << "Migrated Xpetra vector" << std::endl;

    newV->describe(*outStream,v);
  }

  // XpetraTraits<Xpetra::MultiVector<zscalar_t, zlno_t, zgno_t, znode_t>>
  {
    RCP<xmvector_t> MV;

    try{
      RCP<tmvector_t> tMV =
          rcp(new tmvector_t(uinput->getUITpetraCrsGraph()->getRowMap(), 3));
      tMV->randomize();
      MV = Zoltan2::XpetraTraits<tmvector_t>::convertToXpetra(tMV);
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "getXpetraMultiVector", 1);

    if (rank== 0)
      std::cout << "Original Xpetra multivector" << std::endl;

    MV->describe(*outStream,v);

    ArrayRCP<zgno_t> newRowIds = roundRobinMap(*(MV->getMap()));

    zgno_t localNumRows = newRowIds.size();

    RCP<const xmvector_t> newMV;
    try{
      newMV = Zoltan2::XpetraTraits<xmvector_t>::doMigration(*MV,
        localNumRows, newRowIds.getRawPtr());
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok,
        " Zoltan2::XpetraTraits<xmvector_t>::doMigration ", 1);

    if (rank== 0)
      std::cout << "Migrated Xpetra multivector" << std::endl;

    newMV->describe(*outStream,v);
  }

#ifdef HAVE_EPETRA_DATA_TYPES
  /////////////////////////////////////////////////////////////////
  //   Epetra_CrsMatrix
  //   Epetra_CrsGraph
  //   Epetra_Vector
  //   Epetra_MultiVector
  /////////////////////////////////////////////////////////////////

  typedef Epetra_CrsMatrix ematrix_t;
  typedef Epetra_CrsGraph egraph_t;
  typedef Epetra_Vector evector_t;
  typedef Epetra_MultiVector emvector_t;
  typedef Epetra_BlockMap emap_t;

  // XpetraTraits<Epetra_CrsMatrix>
  {
    RCP<ematrix_t> M;

    try{
      M = uinput->getUIEpetraCrsMatrix();
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "getEpetraCrsMatrix ", 1);

    if (rank== 0)
      std::cout << "Original Epetra matrix" << std::endl;

    M->Print(std::cout);

    RCP<const emap_t> emap = Teuchos::rcpFromRef(M->RowMap());
    ArrayRCP<zgno_t> newRowIds = roundRobinMap(*emap);

    zgno_t localNumRows = newRowIds.size();

    RCP<const ematrix_t> newM;
    try{
      newM = Zoltan2::XpetraTraits<ematrix_t>::doMigration(*M,
        localNumRows, newRowIds.getRawPtr());
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok,
        " Zoltan2::XpetraTraits<ematrix_t>::doMigration ", 1);

    if (rank== 0)
      std::cout << "Migrated Epetra matrix" << std::endl;

    newM->Print(std::cout);
  }

  // XpetraTraits<Epetra_CrsGraph>
  {
    RCP<egraph_t> G;

    try{
      G = uinput->getUIEpetraCrsGraph();
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "getEpetraCrsGraph ", 1);

    if (rank== 0)
      std::cout << "Original Epetra graph" << std::endl;

    G->Print(std::cout);

    RCP<const emap_t> emap = Teuchos::rcpFromRef(G->RowMap());
    ArrayRCP<zgno_t> newRowIds = roundRobinMap(*emap);

    zgno_t localNumRows = newRowIds.size();

    RCP<const egraph_t> newG;
    try{
      newG = Zoltan2::XpetraTraits<egraph_t>::doMigration(*G,
        localNumRows, newRowIds.getRawPtr());
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok,
        " Zoltan2::XpetraTraits<egraph_t>::doMigration ", 1);

    if (rank== 0)
      std::cout << "Migrated Epetra graph" << std::endl;

    newG->Print(std::cout);
  }

  // XpetraTraits<Epetra_Vector>
  {
    RCP<evector_t> V;

    try{
      V = rcp(new Epetra_Vector(uinput->getUIEpetraCrsGraph()->RowMap()));
      V->Random();
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "getEpetraVector", 1);

    if (rank== 0)
      std::cout << "Original Epetra vector" << std::endl;

    V->Print(std::cout);

    RCP<const emap_t> emap = Teuchos::rcpFromRef(V->Map());
    ArrayRCP<zgno_t> newRowIds = roundRobinMap(*emap);

    zgno_t localNumRows = newRowIds.size();

    RCP<const evector_t> newV;
    try{
      newV = Zoltan2::XpetraTraits<evector_t>::doMigration(*V,
        localNumRows, newRowIds.getRawPtr());
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok,
        " Zoltan2::XpetraTraits<evector_t>::doMigration ", 1);

    if (rank== 0)
      std::cout << "Migrated Epetra vector" << std::endl;

    newV->Print(std::cout);
  }

  // XpetraTraits<Epetra_MultiVector>
  {
    RCP<emvector_t> MV;

    try{
      MV =
        rcp(new Epetra_MultiVector(uinput->getUIEpetraCrsGraph()->RowMap(),3));
      MV->Random();
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok, "getEpetraMultiVector", 1);

    if (rank== 0)
      std::cout << "Original Epetra multivector" << std::endl;

    MV->Print(std::cout);

    RCP<const emap_t> emap = Teuchos::rcpFromRef(MV->Map());
    ArrayRCP<zgno_t> newRowIds = roundRobinMap(*emap);

    zgno_t localNumRows = newRowIds.size();

    RCP<const emvector_t> newMV;
    try{
      newMV = Zoltan2::XpetraTraits<emvector_t>::doMigration(*MV,
        localNumRows, newRowIds.getRawPtr());
    }
    catch(std::exception &e){
      aok = false;
      std::cout << e.what() << std::endl;
    }
    TEST_FAIL_AND_EXIT(*comm, aok,
        " Zoltan2::XpetraTraits<emvector_t>::doMigration ", 1);

    if (rank== 0)
      std::cout << "Migrated Epetra multivector" << std::endl;

    newMV->Print(std::cout);
  }
#endif   // have epetra data types (int, int, double)

  /////////////////////////////////////////////////////////////////
  // DONE
  /////////////////////////////////////////////////////////////////

  if (rank==0)
    std::cout << "PASS" << std::endl;
}

