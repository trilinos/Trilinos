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

#include <Teuchos_DefaultComm.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <MatrixMarket_Tpetra.hpp>

#include <string>
#include <iostream>
#include <exception>

using std::string;
using std::cerr;
using std::endl;
using Teuchos::RCP;

typedef Kokkos::DefaultNode::DefaultNodeType node_t;
typedef float scalar_t;
typedef int lno_t;
typedef int gno_t;

typedef Tpetra::CrsMatrix<scalar_t, lno_t, gno_t, node_t> tcrsMatrix_t;


int main(int argc, char *argv[])
{
  Teuchos::GlobalMPISession session(&argc, &argv, NULL);
  RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

  // Error: it's a pattern matrix.  (A graph, there are non non-zeros.)
  // packages/tpetra/inout/MatrixMarket_Tpetra.hpp line 1003
  //
  // string fname("commanche_dual.mtx"); 

  // Crash: 
  string fname("USAir97.mtx");

  //string fname("simple.mtx"); // This file is read without error

  RCP<tcrsMatrix_t> M;
  RCP<Kokkos::DefaultNode::DefaultNodeType> dnode
    = Kokkos::DefaultNode::getDefaultNode();

  try{
    M = Tpetra::MatrixMarket::Reader<tcrsMatrix_t>::readSparseFile(
      fname, comm, dnode);
  }
  catch (std::invalid_argument &e){
    cerr << e.what() << endl;
  }
  catch (...){
    cerr << "error" << endl;
  }
}
