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

#ifndef TPETRA_MPIDATA_HPP
#define TPETRA_MPIDATA_HPP

#include <mpi.h>
#include "Tpetra_Object.hpp"

namespace Tpetra {
	
  //! MpiData: Inner data class for MpiPlatform, MpiComm, and MpiDistributor.
  /*! Basically a wrapper for an MPI Communicator object which can be shared
      among multiple Tpetra objects.
  */

  class MpiData : public Object {
  public:
    // default constructor
    MpiData(MPI_Comm Comm)
      : Object("Tpetra::MpiData")
      , MpiComm_(Comm)
      , curTag_(minTag_)
  	{
      // we would prefer to do this in a member initialization, but that's not possible
      MPI_Comm_size(Comm, &size_);
      MPI_Comm_rank(Comm, &rank_);
    };
    
    // destructor
    ~MpiData() {};

    int getMyImageID() const {return(rank_);};
    int getNumImages() const {return(size_);};
    MPI_Comm getMpiComm() const {return(MpiComm_);};
    int getMpiTag() {
      if(curTag_ > maxTag_)
        curTag_ = minTag_;
      return(curTag_++);
    };
    
  private:
    MPI_Comm MpiComm_;
    int rank_;
    int size_;
    
    static int const minTag_ = 26000; // these are somewhat
    static int const maxTag_ = 26099; // magic numbers
    int curTag_;
    
    //! Copy constructor (declared but not defined, do not use)
    MpiData(MpiData const& rhs);
    //! Assignment operator (declared but not defined, do not use)
    MpiData& operator = (MpiData const& rhs);
    
  }; // class MpiData
  
} // namespace Tpetra

#endif // TPETRA_MPIDATA_HPP
