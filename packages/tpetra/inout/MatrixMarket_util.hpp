//@HEADER
// ************************************************************************
// 
//               Tpetra: Linear Algebra Services Package 
//                 Copyright 2011 Sandia Corporation
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

#ifndef __MatrixMarket_util_hpp
#define __MatrixMarket_util_hpp

#include <string>

namespace Tpetra {
  namespace MatrixMarket {

    namespace details {

      /// \class SetScientific
      /// \brief Politely make an output stream use scientific notation.
      /// \author Mark Hoemmen
      ///
      /// On construction, apply the necessary flags to the given
      /// output stream so that floating-point numbers are written in
      /// scientific notation with precision appropriate for the
      /// Scalar type.  On destruction, restore the original
      /// (pre-construction) flags to the output stream.  This makes
      /// SetScientific good for scope-protected alteration of the
      /// output stream's flags; no matter how the scope exits
      /// (normally or by a thrown exception), the original flags will
      /// be restored.  Hence, "polite": we restore the original flags
      /// on scope exit.
      ///
      /// \tparam Scalar A floating-point type, either real or
      ///   complex, for which Teuchos::ScalarTraits<Scalar> has a
      ///   specialization.  Currently we also require that
      ///   std::log10() take arguments of type Scalar, but this may
      ///   be relaxed in the future if Teuchos::ScalarTraits gets its
      ///   own log10() class method.
      template<class Scalar>
      class SetScientific {
      public:
	typedef Scalar scalar_type;

	/// \brief Constructor.
	///
	/// \param out [in/out] Output stream to which to apply the
	///   scientific notation flags.
	SetScientific (std::ostream& out) : 
	  out_ (out), originalFlags_ (out.flags()) 
	{
	  typedef Teuchos::ScalarTraits<scalar_type> STS;
	  typedef typename STS::magnitudeType magnitude_type;
	  typedef Teuchos::ScalarTraits<magnitude_type> STM;

	  // Print floating-point values in scientific notation.
	  out << std::scientific;

	  // We're writing decimal digits, so compute the number of
	  // digits we need to get reasonable accuracy when reading
	  // values back in.
	  //
	  // There is actually an algorithm, due to Guy Steele (yes,
	  // _that_ Guy Steele) et al., for idempotent printing of
	  // finite-length floating-point values.  We should actually
	  // implement that algorithm, but I don't have time for that
	  // now.  Currently, I just print no more than (one decimal
	  // digit more than (the number of decimal digits justified
	  // by the precision of magnitude_type)).

	  // FIXME (mfh 06 Apr 2011) This will only work if
	  // std::log10() is defined for magnitude_type inputs.
	  // Teuchos::ScalarTraits does not currently have an log10()
	  // class method.
	  const magnitude_type numDecDigits = STM::t() * std::log10 (STM::base());

	  // Round and add one. Hopefully this doesn't overflow...
	  const magnitude_type one = STM::one();
	  const magnitude_type two = one + one;
	  const int prec = 1 + static_cast<int> ((two*numDecDigits + one) / two);
	  
	  // Set the number of (decimal) digits after the decimal
	  // point to print.
	  out.precision (prec);
	}	

	/// \brief Destructor.
	///
	/// The destructor sets the output stream's flags back to
	/// their original state, that is, the state before the
	/// constructor of this object was called.
	~SetScientific () {
	  out_.flags (originalFlags_);
	}

      private:
	//! The output stream to which to apply flags.
	std::ostream& out_;

	//! The output stream's original flags.
	std::ios_base::fmtflags originalFlags_;
      };

      // We define a class because template functions can't (in the
      // current C++ standard) have default template parameters.
      template<class Scalar, bool isComplex=Teuchos::ScalarTraits<Scalar>::isComplex>
      class ScalarAssigner {
      public:
	static void
	assign (Scalar& val,
		const typename Teuchos::ScalarTraits<Scalar>::magnitudeType& real,
		const typename Teuchos::ScalarTraits<Scalar>::magnitudeType& imag);
      };

      template<class RealType>
      class ScalarAssigner<RealType, false> {
      public:
	static void
	assign (RealType& val,
		const typename Teuchos::ScalarTraits<RealType>::magnitudeType& real,
		const typename Teuchos::ScalarTraits<RealType>::magnitudeType& imag)
	{
	  // imag had better be zero.  We're ignoring it regardless.
	  (void) imag;
	  val = real; 
	}
      };

#ifdef HAVE_TEUCHOS_COMPLEX
      template<class MagType>
      class ScalarAssigner<std::complex<MagType>, true> {
      public:
	static void
	assign (std::complex<MagType>& val,
		const typename Teuchos::ScalarTraits<std::complex<MagType> >::magnitudeType& real,
		const typename Teuchos::ScalarTraits<std::complex<MagType> >::magnitudeType& imag)
	{
	  val = std::complex<MagType> (real, imag);
	}
      };
#endif // HAVE_TEUCHOS_COMPLEX

      // \fn assignScalar
      // \brief val = S(real, imag).
      //
      // We have to template it because we don't know that S is a
      // complex type; if we write S(real,imag), the compiler will
      // complain if S is a real type.
      template<class Scalar>
      void 
      assignScalar (Scalar& val, 
		    const typename Teuchos::ScalarTraits<Scalar>::magnitudeType& real,
		    const typename Teuchos::ScalarTraits<Scalar>::magnitudeType& imag)
      {
	ScalarAssigner<Scalar>::assign (val, real, imag);
      }

    } // namespace details


    static bool isSkew (const std::string& symmType) {
      return symmType.size() >= 4 && symmType.substr(0,4) == "skew";
    }
    static bool isConj (const std::string& symmType) {
      return std::string::npos != symmType.find ("hermitian");
    }
    static bool needsSymmetrization (const std::string& symmType) {
      return symmType != "general";
    }

    /// \class SymmetrizingAdder
    /// \author Mark Hoemmen
    /// \brief Adds entries with optional symmetry to a sparse matrix
    ///
    /// This class wraps any existing class (AdderType) that defines the
    /// index_type and value_type typedefs, and a "void operator()
    /// (const index_type, const index_type, const value_type&)" (that
    /// conceptually adds an entry to a sparse matrix).  Given the
    /// Matrix Market symmetry type, this class' corresponding
    /// operator() may invoke AdderType's operator() twice, in order to
    /// add entry (j,i) if entry (i,j) is to be added.
    template<class AdderType>
    class SymmetrizingAdder {
    public:
      typedef typename AdderType::index_type index_type;
      typedef typename AdderType::value_type value_type;
    
      SymmetrizingAdder (const Teuchos::RCP<AdderType>& adder, 
			 const std::string& symmType) :
	adder_ (adder),
	symmetrize_ (needsSymmetrization (symmType)),
	conjugate_ (isConj (symmType)),
	skew_ (isSkew (symmType))
      {}


    
      void operator() (const index_type i, const index_type j, const value_type& Aij) {
	AdderType& theAdder = *adder_;

	theAdder (i, j, Aij);
	if (symmetrize_ && i != j) {
	  typedef Teuchos::ScalarTraits<value_type> STS;
	  const value_type Aji = skew_ ? 
	    -(conjugate_ ? STS::conjugate(Aij) : Aij) : 
	    +(conjugate_ ? STS::conjugate(Aij) : Aij);
	  theAdder (j, i, Aji);
	}
      }

      /// \brief Persisting non-const view of the underlying adder object.
      ///
      /// This violates encapsulation, so please be careful with this.
      Teuchos::RCP<AdderType> getAdder() const {
	return adder_;
      }

    private:
      Teuchos::RCP<AdderType> adder_;
      bool symmetrize_, conjugate_, skew_;
    };

  } // namespace MatrixMarket
} // namespace Tpetra

#endif // __MatrixMarket_util_hpp
