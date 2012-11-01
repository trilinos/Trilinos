// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_DefaultComm.hpp"
#include "Teuchos_CommHelpers.hpp"
#include "Teuchos_GlobalMPISession.hpp"

int main(int argc, char* argv[]) {

  using Teuchos::RCP;

  Teuchos::GlobalMPISession mpiSession(&argc, &argv);

  int ierr = 0;

  RCP<const Teuchos::Comm<int> >
    comm = Teuchos::DefaultComm<int>::getComm();

  const int size = Teuchos::size(*comm);
  const int rank = Teuchos::rank(*comm);

  if (rank == 0) {
    std::cout
      << "Rank: " << rank 
      << "\t\tSize: " << size << std::endl;
  }

  const int one  = Teuchos::OrdinalTraits<int>::one();

  std::vector<int> sendbuf(size,one);
  std::vector<int> recvcounts(size,one);

  // Try straight MPI version
#ifdef HAVE_MPI
  {
    std::vector<int> sbuf(sendbuf);
    std::vector<int> rcnt(recvcounts);
    int rbuf;

    MPI_Comm mpi_comm = MPI_COMM_WORLD;
    MPI_Reduce_scatter(&sbuf[0],&rbuf,&rcnt[0],MPI_INT,MPI_SUM,mpi_comm);

    // should be size
    if (rank == 0) {
      std::cout << "Received from direct MPI call: " << rbuf << std::endl;
    }
    ierr += (rbuf == size ? 0 : 1);
  }
#endif

  // Try straight MPI version with user-defined type
#ifdef HAVE_MPI
  {
    std::vector<int> sbuf(sendbuf);
    std::vector<int> rcnt(recvcounts);
    int rbuf;

    MPI_Comm mpi_comm = MPI_COMM_WORLD;
    MPI_Datatype _chars_type;
    MPI_Type_contiguous(sizeof(int), MPI_CHAR, &_chars_type);
    MPI_Reduce_scatter(&sbuf[0], &rbuf, &rcnt[0], _chars_type, MPI_SUM, mpi_comm);

    // should be size
    if (rank == 0) {
      std::cout << "Received from direct MPI call with user type: " << rbuf << std::endl;
    }
    ierr += (rbuf == size ? 0 : 1);
  }
#endif

  // Try Teuchos version
  { 
    std::vector<int> sbuf(sendbuf);
    std::vector<int> rcnt(recvcounts);
    int rbuf;

    Teuchos::reduceAllAndScatter<int>(*comm,Teuchos::REDUCE_SUM,(int)sbuf.size(),&sbuf[0],&rcnt[0],&rbuf);

    // should be size
    if (rank == 0) {
      std::cout << "Received from MPI-via-Teuchos call: " << rbuf << std::endl;
    }
    ierr += (rbuf == size ? 0 : 1);
  }

  if (rank == 0) {
    if (ierr) {
      std::cerr << "Test FAILED." << std::endl;
    }
    else {
      std::cerr << "Test passed." << std::endl;
    }
  }
  return ierr;
}
