/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER
*/

#include <Tpetra_ConfigDefs.hpp>
#include <Teuchos_Array.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Version.hpp"

typedef double scalar_type;
typedef int local_ordinal_type;
//typedef long global_ordinal_type;  //<<<<<<<<   valgrind is clean
typedef int global_ordinal_type;     //<<<<<<<<   valgrind complains
typedef Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;
typedef Tpetra::Map<local_ordinal_type, global_ordinal_type> map_type;
typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type> mat_type;

void
GetNeighboursCartesian2d (const global_ordinal_type i,
                          const global_ordinal_type nx,
                          const global_ordinal_type ny,
                          global_ordinal_type& left,
                          global_ordinal_type& right,
                          global_ordinal_type& lower,
                          global_ordinal_type& upper)
{
  global_ordinal_type ix, iy;
  ix = i % nx;
  iy = (i - ix) / nx;

  if (ix == 0)      left  = -1;
  else              left  = i - 1;
  if (ix == nx - 1) right = -1;
  else              right = i + 1;
  if (iy == 0)      lower = -1;
  else              lower = i - nx;
  if (iy == ny - 1) upper = -1;
  else              upper = i + nx;
}

int
main (int argc, char *argv[])
{
  // global_size_t: Tpetra defines this unsigned integer type big
  // enough to hold any global dimension or amount of data.
  typedef Tpetra::global_size_t GST;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackhole);
  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  const size_t myRank = comm->getRank ();

  // global_ordinal_type nx=10, ny=10;
  // scalar_type a=4,b=-1,c=-1,d=-1,e=-1;

  Array<global_ordinal_type> myElts;

  switch (myRank) {
  case 0:
    myElts.push_back(0);
    myElts.push_back(1);
    break;
  case 1:
    myElts.push_back(2);
    myElts.push_back(3);
    break;
  }

  const GST GST_INV = Teuchos::OrdinalTraits<GST>::invalid ();
  const global_ordinal_type indexBase = 0;
  RCP<const map_type> map =
    rcp (new map_type (GST_INV, myElts (), indexBase, comm));

  RCP<Teuchos::FancyOStream> fos =
    Teuchos::fancyOStream (Teuchos::rcpFromRef (std::cout));
  fos->setOutputToRootOnly (-1);
  map->describe (*fos, Teuchos::VERB_EXTREME);

  const size_t numMyElements = map->getNodeNumElements ();
  switch (myRank) {
  case 0:
    assert(numMyElements==2);
    break;
  case 1:
    assert(numMyElements==2);
    break;
  }

  // FIXME (mfh 19 Mar 2014) Once you set this
  // ("setOutputToRootOnly"), you can't unset it.
  fos->setOutputToRootOnly (0);
  *fos << std::endl << "Creating the sparse matrix" << std::endl;

  const local_ordinal_type nnz = 4;
  RCP<mat_type> A = rcp (new mat_type (map, nnz));

  ArrayView<const global_ordinal_type> myGlobalElements =
    map->getNodeElementList ();

  // global_ordinal_type center, left, right, lower, upper;

  Teuchos::Array<scalar_type> vals (nnz);
  Teuchos::Array<global_ordinal_type> inds (nnz);
  for (int i = 0; i < nnz; ++i) {
    inds[i] = i;
    vals[i] = 1.0;
  }

  comm->barrier ();
  for (int k = 0; k< comm->getSize (); ++k) {
    if (comm->getRank () == k) {
      std::cout << "pid " << k << " inserting global rows" << std::endl;
      for (size_t i = 0; i < numMyElements; ++i)  {
        //size_t n = 0;

        std::cout << "   grow " << myGlobalElements[i] << " : ";
        for (int jj = 0; jj < 4; ++jj) {
          std::cout << inds[jj] << " ";
        }
        std::cout << std::endl;
        A->insertGlobalValues (myGlobalElements[i], inds (), vals ());
      }
    }
    // mfh 01 Apr 2014: sleep() is a POSIX function, not a C++
    // standard function, and thus not suitable for Tpetra tests.
    //
    //sleep (1);
    comm->barrier ();
  } //k

  A->fillComplete ();

  // Make sure that all processes finished fillComplete, before
  // reporting success.
  comm->barrier ();
  if (comm->getRank () == 0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
    return EXIT_SUCCESS;
  }
}


