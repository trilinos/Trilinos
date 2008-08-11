// @HEADER
// ***********************************************************************
// 
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
