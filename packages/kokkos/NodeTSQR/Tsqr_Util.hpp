//@HEADER
// ************************************************************************
// 
//          Kokkos: Node API and Parallel Node Kernels
//              Copyright (2009) Sandia Corporation
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
// ************************************************************************
//@HEADER

#ifndef __TSQR_Tsqr_Util_hpp
#define __TSQR_Tsqr_Util_hpp

#include "Tsqr_ScalarTraits.hpp"

#include <algorithm>
#include <complex>
#include <ostream>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  /// \brief Print a Scalar value to the given output stream
  /// 
  /// C++ annoyingly doesn't let me do partial template specialization
  /// of functions.  Because of that, I can't use a template function;
  /// instead, I have to reify the function into a class ("function
  /// object"), in the worst Java style (where everything is a noun
  /// with a "run()" method).
  ///
  template< class Scalar, bool isComplex >
  class ScalarPrinter {
  public:
    ///
    /// Print elt to out
    void operator() (std::ostream& out, const Scalar& elt) const;
  };

  // Partial specialization for real Scalar
  template< class Scalar >
  class ScalarPrinter< Scalar, false > {
  public:
    void operator() (std::ostream& out, const Scalar& elt) const {
      out << elt;
    }
  };

  // Partial specialization for complex Scalar
  template< class Scalar >
  class ScalarPrinter< Scalar, true > {
  public:
    void operator() (std::ostream& out, const Scalar& elt) const {
      typedef typename ScalarTraits< Scalar >::magnitude_type magnitude_type;

      const magnitude_type ZERO (0);
      const magnitude_type& realPart = std::real (elt);
      const magnitude_type& imagPart = std::imag (elt);

      out << realPart;
      if (imagPart < ZERO)
	out << "-" << ScalarTraits< Scalar >::abs(imagPart) << "*i";
      else if (imagPart > ZERO)
	out << "+" << imagPart << "*i";
    }
  };

  template< class LocalOrdinal, class Scalar >
  void
  print_local_matrix (std::ostream& out,
		      const LocalOrdinal nrows_local,
		      const LocalOrdinal ncols,
		      const Scalar A[],
		      const LocalOrdinal lda)
  {
    ScalarPrinter< Scalar, ScalarTraits< Scalar >::is_complex > printer;
    for (LocalOrdinal i = 0; i < nrows_local; i++)
      {
	for (LocalOrdinal j = 0; j < ncols; j++)
	  {
	    const Scalar& curElt = A[i + j*lda];
	    printer (out, curElt);
	    if (j < ncols - 1)
	      out << ", ";
	  }
	out << ";" << std::endl;
      }
  }

  template< class Ordinal, class Scalar >
  void
  copy_matrix (const Ordinal nrows,
	       const Ordinal ncols,
	       Scalar* const A,
	       const Ordinal lda,
	       const Scalar* const B,
	       const Ordinal ldb)
  {
    for (Ordinal j = 0; j < ncols; j++)
      {
	Scalar* const A_j = &A[j*lda];
	const Scalar* const B_j = &B[j*ldb];
	std::copy (B_j, B_j + nrows, A_j);
      }
  }

  template< class Ordinal, class Scalar >
  void
  fill_matrix (const Ordinal nrows,
	       const Ordinal ncols,
	       Scalar* const A,
	       const Ordinal lda,
	       const Scalar& default_val)
  {
    for (Ordinal j = 0; j < ncols; j++)
      {
	Scalar* const A_j = &A[j*lda];
	std::fill (A_j, A_j + nrows, default_val);
      }
  }


  template< class Ordinal, class Scalar, class Generator >
  void
  generate_matrix (const Ordinal nrows,
		   const Ordinal ncols,
		   Scalar* const A,
		   const Ordinal lda,
		   Generator gen)
  {
    for (Ordinal j = 0; j < ncols; j++)
      {
	Scalar* const A_j = &A[j*lda];
	std::generate (A_j, A_j + nrows, gen);
      }
  }

  template< class Ordinal, class Scalar >
  void
  copy_upper_triangle (const Ordinal nrows,
		       const Ordinal ncols,
		       Scalar* const R_out,
		       const Ordinal ldr_out,
		       const Scalar* const R_in,
		       const Ordinal ldr_in)
  {
    if (nrows >= ncols)
      {
	for (Ordinal j = 0; j < ncols; j++)
	  {
	    Scalar* const A_j = &R_out[j*ldr_out];
	    const Scalar* const B_j = &R_in[j*ldr_in];
	    for (Ordinal i = 0; i <= j; i++)
	      A_j[i] = B_j[i];
	  }
      }
    else
      {
	copy_upper_triangle (nrows, nrows, R_out, ldr_out, R_in, ldr_in);
	for (Ordinal j = nrows; j < ncols; j++)
	  {
	    Scalar* const A_j = &R_out[j*ldr_out];
	    const Scalar* const B_j = &R_in[j*ldr_in];
	    for (Ordinal i = 0; i < nrows; i++)
	      A_j[i] = B_j[i];
	  }
      }
  }


  template< class Scalar >
  class SumSquare {
  public:
    Scalar operator() (const Scalar& result, const Scalar& x) const { 
      return result + x*x; 
    }
  };

  // Specialization for complex numbers
  template< class Scalar >
  class SumSquare< std::complex< Scalar > >  {
  public:
    Scalar operator() (const std::complex<Scalar>& result, 
		       const std::complex<Scalar>& x) const { 
      const Scalar absval = std::norm (x);
      return result + absval * absval; 
    }
  };

  template< class Ordinal, class Scalar >
  void
  pack_R_factor (const Ordinal nrows, 
		 const Ordinal ncols, 
		 const Scalar R_in[], 
		 const Ordinal ldr_in,
		 Scalar buffer[])
  {
    Ordinal count = 0; // current position in output buffer
    if (nrows >= ncols)
      for (Ordinal j = 0; j < ncols; j++)
	for (Ordinal i = 0; i <= j; i++)
	  buffer[count++] = R_in[i + j*ldr_in];
    else
      for (Ordinal j = 0; j < nrows; j++)
	for (Ordinal i = 0; i <= j; i++)
	  buffer[count++] = R_in[i + j*ldr_in];
  }

  template< class Ordinal, class Scalar >
  void
  unpack_R_factor (const Ordinal nrows, 
		   const Ordinal ncols, 
		   Scalar R_out[], 
		   const Ordinal ldr_out,
		   const Scalar buffer[])
  {
    Ordinal count = 0; // current position in input buffer
    if (nrows >= ncols)
      for (Ordinal j = 0; j < ncols; j++)
	for (Ordinal i = 0; i <= j; i++)
	  R_out[i + j*ldr_out] = buffer[count++];
    else
      for (Ordinal j = 0; j < nrows; j++)
	for (Ordinal i = 0; i <= j; i++)
	  R_out[i + j*ldr_out] = buffer[count++];
  }

} // namespace TSQR

#endif // __TSQR_Tsqr_Util_hpp
