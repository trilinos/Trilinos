/*
// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER
*/

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_as.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_Comm.hpp"
#include "Teuchos_UnitTestHarness.hpp"


namespace {


//
// Helper code and declarations
//

using Teuchos::as;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Ptr;
using Teuchos::outArg;
using Teuchos::Array;
using Teuchos::Comm;
typedef Teuchos_Ordinal Ordinal;


int g_localDim = 4;
bool g_dumpAll = false;
bool g_show_all_tests = false;


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "local-dim", &g_localDim, "Local dimension of each vector." );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "show-all-tests", "no-show-all-tests", &g_show_all_tests,
    "Set if all tests are shown or not." );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "dump-all", "no-dump-all", &g_dumpAll,
    "Dump lots of data" );
}


RCP<const Epetra_Comm> getEpetraComm()
{
#ifdef HAVE_MPI
  return rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  return rcp(new Epetra_SerialComm());
#endif
}


void createEpetraVsAndMap(const Thyra::Ordinal localDim_in,
  const Ptr<RCP<const Thyra::VectorSpaceBase<double> > > &vs,
  const Ptr<RCP<const Epetra_Map> > &epetra_map,
  const int emptyProcRootRank = -1
  )
{
  const RCP<const Epetra_Comm> epetra_comm = getEpetraComm();
  const int procRank = epetra_comm->MyPID();
  const Thyra::Ordinal localDim = (procRank == emptyProcRootRank ? 0 : localDim_in);
  *epetra_map = rcp(new Epetra_Map(-1, as<int>(localDim), 0, *epetra_comm));
  *vs =  Thyra::create_VectorSpace(*epetra_map);
}


RCP<Epetra_CrsMatrix> getEpetraMatrix(int numRows, int numCols, double shift=0.0) 
{

  const RCP<const Epetra_Comm> comm = getEpetraComm();

  const Epetra_Map rowMap(numRows, 0, *comm);
  const Epetra_Map domainMap(numCols, numCols, 0, *comm);
 
  const RCP<Epetra_CrsMatrix> epetraCrsM =
    rcp(new Epetra_CrsMatrix(Copy, rowMap, numCols));

  Array<double> rowEntries(numCols);
  Array<int> columnIndices(numCols);
  for (int j = 0; j < numCols; ++j) {
    columnIndices[j] = j;
  }

  const int numLocalRows = rowMap.NumMyElements();

  for (int i = 0; i < numLocalRows; ++i) {
    
    for (int j = 0; j < numCols; ++j) {
      rowEntries[j] = as<double>(i+1) + as<double>(j+1) / 10 + shift;
    }

    epetraCrsM->InsertMyValues( i, numCols, &rowEntries[0], &columnIndices[0] );

  }

  epetraCrsM->FillComplete(domainMap, rowMap);

  return epetraCrsM;
}


} // namespace
