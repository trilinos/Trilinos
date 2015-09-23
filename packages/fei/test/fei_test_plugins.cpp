/*
// @HEADER
// ************************************************************************
//             FEI: Finite Element Interface to Linear Solvers
//                  Copyright (2005) Sandia Corporation.
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation, the
// U.S. Government retains certain rights in this software.
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
// Questions? Contact Alan Williams (william@sandia.gov) 
//
// ************************************************************************
// @HEADER
*/


#include <fei_iostream.hpp>

#include <fei_utils.hpp>
#include <test_utils/fei_test_utils.hpp>

#include <test_utils/test_Factory.hpp>
#include <test_utils/test_Vector.hpp>
#include <test_utils/test_Matrix.hpp>

#include <fei_Factory_Trilinos.hpp>

#undef fei_file
#define fei_file "fei_test_plugins.cpp"
#include <fei_ErrMacros.hpp>


//--------------------------------------------------
// main
//--------------------------------------------------
int main(int argc, char** argv) {

  int numProcs, localProc;
  CHK_ERR( fei_test_utils::initialize_mpi(argc, argv, localProc, numProcs) );

  if (localProc == 0) {
    FEI_COUT << "FEI version: " << fei::utils::version() << FEI_ENDL;
  }

  fei::SharedPtr<fei::Factory> factory(new Factory_Trilinos(MPI_COMM_WORLD));

  fei::ParameterSet pset;
  pset.add(fei::Param("LPM_EpetraBasic", true));
  factory->parameters(pset);

  test_Factory ftester(MPI_COMM_WORLD);

  ftester.factory_test1(factory);

  test_Vector vtester(MPI_COMM_WORLD);

  fei::SharedPtr<fei::Vector> fei_vec = vtester.create_vector(factory);

  vtester.vector_test1(fei_vec);

  test_Matrix mtester(MPI_COMM_WORLD);

  fei::SharedPtr<fei::Matrix> fei_mat = mtester.create_matrix(factory);

  mtester.matrix_test1(fei_mat);

#ifndef FEI_SER
  if (MPI_Finalize() != MPI_SUCCESS) ERReturn(-1);
#endif

  return(0);
}

