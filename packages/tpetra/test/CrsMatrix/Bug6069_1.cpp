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

// Some Macro Magic to ensure that if CUDA and KokkosCompat is enabled
// only the .cu version of this file is actually compiled
#include <Tpetra_config.h>

#include <Teuchos_Array.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

typedef double scalar_type;
typedef int local_ordinal_type;
//typedef long global_ordinal_type;  //<<<<<<<<   valgrind is clean
typedef int global_ordinal_type;     //<<<<<<<<   valgrind complains

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
  using Tpetra::global_size_t;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::ArrayRCP;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Put these typedefs here, to avoid global shadowing warnings.
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type> map_type;
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type, global_ordinal_type> mat_type;

  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackhole);
  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();
  const size_t myRank = comm->getRank ();

  global_ordinal_type nx=10, ny=10;
  scalar_type a=4,b=-1,c=-1,d=-1,e=-1;

  //global_size_t numMyElts;
  Array<global_ordinal_type> myElts;

  switch (myRank) {

  case 0:
    //numMyElts=40;
    for (int i=0; i<10; ++i) {
      for (int j=0; j<4; ++j) {
        myElts.push_back(i+10*j);
      }
    }
    break;

  case 1:
    //numMyElts=30;
    for (int i=40; i<50; ++i) {
      for (int j=0; j<3; ++j) {
        myElts.push_back(i+10*j);
      }
    }
    break;

  case 2:
    //numMyElts=30;
    for (int i=70; i<80; ++i) {
      for (int j=0; j<3; ++j) {
        myElts.push_back(i+10*j);
      }
    }
    break;
  }

/*
  My global indices: [0, 10, 20, 30, 1, 11, 21, 31, 2, 12, 22, 32, 3, 13, 23, 33, 4, 14, 24, 34, 5, 15, 25, 35, 6, 16, 26, 36, 7, 17, 27, 37, 8, 18, 28, 38, 9, 19, 29, 39]
 Process 1:
  My number of entries: 30
  My minimum global index: 40
  My maximum global index: 69
  My global indices: [40, 50, 60, 41, 51, 61, 42, 52, 62, 43, 53, 63, 44, 54, 64, 45, 55, 65, 46, 56, 66, 47, 57, 67, 48, 58, 68, 49, 59, 69]
 Process 2:
  My number of entries: 30
  My minimum global index: 70
  My maximum global index: 99
  My global indices: [70, 80, 90, 71, 81, 91, 72, 82, 92, 73, 83, 93, 74, 84, 94, 75, 85, 95, 76, 86, 96, 77, 87, 97, 78, 88, 98, 79, 89, 99]
*/

  const global_ordinal_type indexBase = 0;
  RCP<const map_type> map = rcp (new map_type (100, myElts (), indexBase, comm));

  RCP<Teuchos::FancyOStream> fos =
    Teuchos::fancyOStream (Teuchos::rcpFromRef (std::cout));
  //fos->setOutputToRootOnly(-1);
  //map->describe(*fos,Teuchos::VERB_EXTREME);

  const size_t numMyElements = map->getNodeNumElements ();
  switch (myRank) {
  case 0:
    assert(numMyElements==40);
    break;
  case 1:
    assert(numMyElements==30);
    break;
  case 2:
    assert(numMyElements==30);
    break;
  }

  // FIXME (mfh 19 Mar 2014) Once you set this
  // ("setOutputToRootOnly"), you can't unset it.
  fos->setOutputToRootOnly(0);
  *fos << std::endl << "Creating the sparse matrix" << std::endl;

  local_ordinal_type nnz=5;
  RCP<mat_type> A = rcp (new mat_type (map, nnz));

  ArrayView<const global_ordinal_type> myGlobalElements = map->getNodeElementList();

  global_ordinal_type center, left, right, lower, upper;
  std::vector<scalar_type>        vals(nnz);
  std::vector<global_ordinal_type> inds(nnz);

  //    e
  //  b a c
  //    d
  for (size_t i = 0; i < numMyElements; ++i)  {
    size_t n = 0;

    center = myGlobalElements[i] - indexBase;
    GetNeighboursCartesian2d(center, nx, ny, left, right, lower, upper);

    if (left  != -1) { inds[n] = left;  vals[n++] = b; }
    if (right != -1) { inds[n] = right; vals[n++] = c; }
    if (lower != -1) { inds[n] = lower; vals[n++] = d; }
    if (upper != -1) { inds[n] = upper; vals[n++] = e; }

    // diagonal
    scalar_type z = a;
    inds[n]   = center;
    vals[n++] = z;

    for (size_t j = 0; j < n; j++)
      inds[j] += indexBase;

    ArrayView<global_ordinal_type> iv(&inds[0], n);
    ArrayView<scalar_type>        av(&vals[0], n);
    A->insertGlobalValues(myGlobalElements[i], iv, av);
  }

  A->fillComplete ();

  // Make sure that all processes finished fillComplete, before
  // reporting success.
  comm->barrier ();
  if (comm->getRank () == 0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
    return EXIT_SUCCESS;
  }
}


