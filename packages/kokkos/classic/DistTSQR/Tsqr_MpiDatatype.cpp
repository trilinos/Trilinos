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

#include <Tsqr_ConfigDefs.hpp>
#include <Tsqr_MpiDatatype.hpp>

#ifdef HAVE_MPI // #defined (or not) via Teuchos_ConfigDefs.hpp
#  ifdef HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
#    include <complex>
#  endif // HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
#  include <utility> // std::pair


namespace TSQR {
  namespace MPI {

    MPI_Datatype 
    cloneRawDatatype (MPI_Datatype in, const bool needsFree)
    {
      MPI_Datatype out;
      if (needsFree)
	{
	  // Non-simple MPI_Datatype objects (i.e., those on which
	  // we need to call MPI_Type_free() after we're done using
	  // them) cannot be copied directly.  Instead, we clone
	  // them using MPI_Type_contiguous() with a count argument
	  // of 1.
	  const int err = MPI_Type_contiguous (1, in, &out);
	  if (err != MPI_SUCCESS)
	    throw std::runtime_error("Failed to clone MPI_Datatype object");
	}
      else
	{
	  // Simple MPI_Datatype objects can be copied directly.
	  // Here, "simple" means that we don't have to call
	  // MPI_Type_free() on them when we are done using them.
	  out = in;
	}
      return out;
    }

    /// Return a new MPI_Datatype, corresponding to a pair of
    /// contiguously placed objects, each of which is represented by
    /// the MPI_Datatype "in".  
    ///
    /// \warning MPI_Type_free() needs to be called on the return
    /// value, once you're done using it.
    ///
    /// \note We need this function because the C MPI binding doesn't
    /// define MPI_2FLOAT and MPI_2DOUBLE, which would be nice things
    /// to have defined... Now, the Fortran MPI binding defines
    /// MPI_2REAL, MPI_2DOUBLE_PRECISION, and MPI_2COMPLEX, but the C
    /// binding doesn't (well, not in the standard, though they may be
    /// defined in particular implementations, I suppose -- but I
    /// haven't found them yet in a particular implementation).
    static MPI_Datatype
    mpi_pair_of (MPI_Datatype in)
    {
      // This is just a handle, remember?  It's ok to return it; it
      // won't fall out of scope.  The MPI implementation is responsible
      // for keeping track of these things.
      MPI_Datatype new_type;
      int err = MPI_Type_contiguous (2, in, &new_type);
      if (err != MPI_SUCCESS)
	throw std::logic_error ("Failed to create MPI_Datatype");
      return new_type;
    }

    template<>
    MpiDatatype< double >::MpiDatatype () : 
      type_ (MPI_DOUBLE),
      needsFree_ (false)
    {}

    template<>
    MpiDatatype< float >::MpiDatatype () :
      type_ (MPI_FLOAT),
      needsFree_ (false)
    {}

    template<>
    MpiDatatype< std::pair< double, double > >::MpiDatatype () : 
      type_ (mpi_pair_of (MPI_DOUBLE)),
      needsFree_ (true)
    {}

    template<>
    MpiDatatype< std::pair< float, float > >::MpiDatatype () : 
      type_ (mpi_pair_of (MPI_FLOAT)),
      needsFree_ (true)
    {}

#ifdef HAVE_KOKKOSCLASSIC_TSQR_COMPLEX
    template<>
    MpiDatatype< std::complex<double> >::MpiDatatype () :
      type_ (mpi_pair_of (MPI_DOUBLE)),
      needsFree_ (true)
    {}

    template<>
    MpiDatatype< std::complex<float> >::MpiDatatype () :
      type_ (mpi_pair_of (MPI_FLOAT)),
      needsFree_ (true)
    {}
#endif // HAVE_KOKKOSCLASSIC_TSQR_COMPLEX

    template<>
    MpiDatatype< int >::MpiDatatype () : 
      type_ (MPI_INT),
      needsFree_ (false)
    {}

    template<>
    MpiDatatype< std::pair<int, int> >::MpiDatatype () : 
      type_ (MPI_2INT),
      needsFree_ (false)
    {}

    template<>
    MpiDatatype< unsigned long >::MpiDatatype () : 
      type_ (MPI_UNSIGNED_LONG),
      needsFree_ (false)
    {}

    template<>
    MpiDatatype< std::pair< unsigned long, unsigned long > >::MpiDatatype () : 
      type_ (mpi_pair_of (MPI_UNSIGNED_LONG)),
      needsFree_ (true)
    {}

  } // namespace MPI
} // namespace TSQR


#endif // HAVE_MPI
