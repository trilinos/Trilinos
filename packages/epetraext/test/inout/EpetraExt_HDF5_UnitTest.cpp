//@HEADER
// ***********************************************************************
//
//     EpetraExt: Epetra Extended - Linear Algebra Services Package
//                 Copyright (2011) Sandia Corporation
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
// ***********************************************************************
//@HEADER

#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_GlobalMPISession.hpp"

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif

#include "EpetraExt_HDF5.h"

namespace EpetraExt {

TEUCHOS_UNIT_TEST(EpetraExt_HDF5, WriteReadInt)
{
#ifdef EPETRA_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  HDF5 file1(comm);
  file1.Create("HDF5_test.h5");

  const int value1 = 5;
  file1.Write("data", "int", value1);
  file1.Close();

  HDF5 file2(comm);
  file2.Open("HDF5_test.h5");

  int value2 = -1;
  file2.Read("data", "int", value2);
  file2.Close();

  TEST_EQUALITY(value1, value2);
}

TEUCHOS_UNIT_TEST(EpetraExt_HDF5, WriteReadSerialData)
{
#ifdef EPETRA_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  HDF5 file1(comm);
  file1.Create("HDF5_test.h5");

  const int data1[2] = {1, 2};
  file1.Write("data", "values", H5T_NATIVE_INT, 2, data1);
  file1.Close();

  HDF5 file2(comm);
  file2.Open("HDF5_test.h5");

  int data2[2] = {-1, -1};
  file2.Read("data", "values", H5T_NATIVE_INT, 2, data2);
  file2.Close();

  TEST_EQUALITY(data1[0], data2[0]);
  TEST_EQUALITY(data1[1], data2[1]);
}

TEUCHOS_UNIT_TEST(EpetraExt_HDF5, NestedGroups)
{
#ifdef EPETRA_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  HDF5 file1(comm);
  file1.Create("HDF5_test.h5");

  const int value1 = 5;
  file1.CreateGroup("group 1");
  file1.CreateGroup("group 1/group 2");
  file1.Write("group 1/group 2/data", "int", value1);
  file1.Close();

  HDF5 file2(comm);
  file2.Open("HDF5_test.h5");

  int value2 = -1;
  file2.Read("group 1/group 2/data", "int", value2);
  file2.Close();

  TEST_EQUALITY(value1, value2);
}

} // namespace EpetraExt

int main(int argc, char* argv[])
{
  Teuchos::GlobalMPISession mpiSession(&argc, &argv);
  return Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
}
