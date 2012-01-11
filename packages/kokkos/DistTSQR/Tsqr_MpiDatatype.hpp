/*
//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2008) Sandia Corporation
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
//@HEADER
*/

#ifndef __TSQR_MpiDatatype_hpp
#define __TSQR_MpiDatatype_hpp

#include <Teuchos_ConfigDefs.hpp> // HAVE_MPI

#ifdef HAVE_MPI
#  include <mpi.h>
#  include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace MPI {

    /// \brief Clone an MPI_Datatype object.
    ///
    MPI_Datatype 
    cloneRawDatatype (MPI_Datatype in, const bool needsFree);

    /// Implements a map from C++ datatype to MPI_Datatype.
    ///
    /// MpiDatatype manages type creation and freeing (via
    /// MPI_Type_free()) automatically, for types for which it's
    /// necessary.  For other types, MpiDatatype::get() just returns
    /// the predefined constant (e.g., MPI_DOUBLE or MPI_INT).  Note
    /// also that MpiDatatype< Datum > has only been defined for
    /// certain types of Datum (see MpiDatatype.cpp for which ones).
    template< class Datum >
    class MpiDatatype {
    public:
      /// Constructor, specialized by hand for typical scalar data
      /// types.
      MpiDatatype ();

      /// Copy constructor
      ///
      MpiDatatype (const MpiDatatype& rhs) 
      {
	clone (this->type_, this->needsFree_, rhs);
      }

      /// Assignment operator
      ///
      MpiDatatype& operator= (const MpiDatatype& rhs)
      {
	if (this != &rhs)
	  {
	    if (needsFree_)
	      {
		// Return value doesn't matter...
		(void) MPI_Type_free (&type_);
		needsFree_ = false;
	      }
	    clone (this->type_, this->needsFree_, rhs);
	  }
	return *this;
      }

      /// Destructor
      ///
      ~MpiDatatype () 
      { 
	if (needsFree_)
	  {
	    // Return value doesn't matter...
	    (void) MPI_Type_free (&type_);
	    needsFree_ = false;
	  }
      }

      MPI_Datatype get() const { return type_; }

    private:
      static void 
      clone (MPI_Datatype& newType, 
	     bool& needsFree,
	     const MpiDatatype& rhs)
      {
	newType = cloneRawDatatype (rhs.get(), rhs.needsFree_);
	needsFree = rhs.needsFree_;
      }

      /// The actual MPI_Datatype object corresponding to Datum.
      ///
      MPI_Datatype type_;

      /// Whether or not MPI_Type_free() needs to be called on type_, in
      /// ~MpiDatatype().
      bool needsFree_;
    };

  } // namespace MPI
} // namespace TSQR


#endif // HAVE_MPI
#endif // __TSQR_MpiDatatype_hpp
