#ifndef __TSQR_Tsqr_Util_hpp
#define __TSQR_Tsqr_Util_hpp

#include <algorithm>
#include <complex>
#include <ostream>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {

  template< class LocalOrdinal, class Scalar >
  void
  print_local_matrix (std::ostream& out,
		      const LocalOrdinal nrows_local,
		      const LocalOrdinal ncols,
		      const Scalar A[],
		      const LocalOrdinal lda)
  {
    for (LocalOrdinal i = 0; i < nrows_local; i++)
      {
	for (LocalOrdinal j = 0; j < ncols; j++)
	  out << A[i + j*lda] << " ";
	out << std::endl;
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
