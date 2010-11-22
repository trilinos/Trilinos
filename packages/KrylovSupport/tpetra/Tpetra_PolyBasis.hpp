#ifndef __Tpetra_PolyBasis_hpp
#define __Tpetra_PolyBasis_hpp

#include <Krylov_PolyBasis.hpp>
#include <Tpetra_MultiVector_decl.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_TypeTraits.hpp>

#include <algorithm>
#ifdef HAVE_TEUCHOS_COMPLEX
#  include <complex>
#endif // HAVE_TEUCHOS_COMPLEX
#include <cctype> // tolower
#include <stdexcept>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace Tpetra {

  /// \class PolyBasis
  /// \brief Interface to the matrix powers kernel
  /// 
  /// This class defines the interface to the matrix powers kernel,
  /// for Tpetra::MultiVector inputs and outputs.
  ///
  /// Template parameters:
  /// - DomainScalar: Scalar type of input multivector
  /// - RangeScalar:  Scalar type of output multivector
  /// - LocalOrdinal: Local ordinal type of input and output multivectors
  /// - GlobalOrdinal: Global ordinal type of input and output multivectors
  /// - Node: Kokkos Node type of input and output multivectors
  template< class DomainScalar, 
	    class RangeScalar,
	    class LocalOrdinal,
	    class GlobalOrdinal,
	    class Node >
  class PolyBasis : public Krylov::PolyBasis< RangeScalar, Tpetra::MultiVector< RangeScalar, LocalOrdinal, GlobalOrdinal, Node >, DomainScalar, Tpetra::MultiVector< DomainScalar, LocalOrdinal, GlobalOrdinal, Node > >
  {
  public:
    typedef DomainScalar domain_scalar_type;
    typedef typename Teuchos::ScalarTraits< domain_scalar_type >::magnitudeType domain_magnitude_type;
    typedef RangeScalar range_scalar_type;
    typedef typename Teuchos::ScalarTraits< range_scalar_type >::magnitudeType range_magnitude_type;
    typedef LocalOrdinal local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node node_type;

    /// \brief Multivector input type.
    ///
    /// PolyBasis defines an interface for a mapping between
    /// MultiVectors.  The input and output MultiVector types may be
    /// the same, or may differ only in their Scalar types.
    typedef MultiVector< domain_scalar_type, local_ordinal_type, global_ordinal_type, node_type > domain_mv_type;

    /// \brief Multivector output type.  
    /// 
    /// Should be the same, except perhaps for a different Scalar
    /// type, as domain_mv_type (see above).  The best reason to have
    /// RangeScalar != DomainScalar is to improve precision of the
    /// matrix powers kernel by accumulating the output vectors in
    /// higher precision than the input.
    typedef MultiVector< range_scalar_type, local_ordinal_type, global_ordinal_type, node_type > range_mv_type;

    virtual ~PolyBasis() = 0;

    /// \brief Compute Y = [A*X, A^2 * X, ..., A^power * X]
    ///
    /// \param X [in]   Input multivector
    /// \param Y [out]  Output multivector:  
    ///   \fn$Y = [A*X, A^2 * X, ..., A^power * X]\fn$.
    ///   In Y, columns of \fn$A^k X\fn$ are adjacent.
    /// \param power [in] Number of times to "apply" the operator A.
    virtual void 
    applyMonomial (const domain_mv_type& X,
		   range_mv_type& Y,
		   const int power) const = 0;

    /// We store possibly complex shifts for the Newton basis as pairs
    /// of magnitudes (real,imag).  This means no std::complex
    /// arithmetic support is needed (even though we wouldn't actually
    /// use complex arithmetic).
    typedef std::pair< range_magnitude_type, range_magnitude_type > shift_type;

    /// Compute 'power' steps of the "modified Newton" basis.
    ///
    /// \param X [in]
    /// \param Y [out]
    /// \param power [in] Number of times to "apply" the operator A
    ///
    virtual void 
    applyModifiedNewton (const domain_mv_type& X,
			 range_mv_type& Y,
			 const std::vector< shift_type >& shifts,
			 const int power) const = 0;

  protected:

    /// Check whether the arguments are valid arguments to methods
    /// like applyMonomial().  If so, do nothing, else throw
    /// std::invalid_argument.
    ///
    /// \note We supply a reasonable default implementation, which
    /// considers power > maxPower an error, and does not attempt to
    /// resize Y.  Other implementations might prefer to accept any
    /// positive integer for the 'power' argument, if necessary
    /// adjusting their internal data structures, and might also
    /// resize Y if it needs more columns.  This is why we make this
    /// method virtual and protected.
    virtual void 
    applyCheck (const domain_mv_type& X,
		range_mv_type& Y,
		const int power,
		const int maxPower) const
    {
      using Teuchos::Range1D;
      using Teuchos::RCP;
      using Teuchos::ScalarTraits;
      
      if (maxPower < 0)
	{
	  std::ostringstream os;
	  os << "Matrix powers operator only supports nonnegative integer powers, but you specified maxPower = " << maxPower;
	}
      else if (power < 0)
	{
	  std::ostringstream os;
	  os << "Matrix powers operator only supports nonnegative integer powers, but you specified power = " << power;
	  throw std::invalid_argument(os.str());
	}
      else if (power > maxPower)
	{
	  std::ostringstream os;
	  os << "Matrix powers operator only supports up to A^{power} with power=" << maxPower << ", but you specified power = " << power;
	  throw std::invalid_argument(os.str());
	}

      // X and Y must have been constructed using the same Map.  If
      // compatibility would suffice, isCompatible would be the
      // corresponding test.
      if (! X.getMap()->isSameAs (Y.getMap()))
	{
	  std::ostringstream os;
	  os << "X and Y must use the same Map in order to work with the matrix powers kernel";
	  throw std::invalid_argument(os.str());
	}

      // Check that Y has enough columns to hold 'power' copies of X.
      // Leave any extra columns of Y alone.
      if (Y.getNumVectors() < power * X.getNumVectors())
	{
	  std::ostringstream os;
	  os << "Y needs at least " << (power * X.getNumVectors()) << " columns, but only has " << Y.getNumVectors() << " columns";
	  throw std::invalid_argument("Y doesn't have enough columns");
	}
    }
  };


  /// \class PolyBasisFactory
  /// \brief Factory for PolyBasis specializations
  ///
  /// Interface to a factory for creating a specialization of
  /// PolyBasis, parametrized by the sparse matrix / operator type
  /// (MatrixType) and the domain and range Scalar types.
  ///
  /// \note We would have used a nonmember function for this, rather
  ///   than a factory class, but C++ (until the latest revision of
  ///   the standard) doesn't allow nonmember template functions with
  ///   default template parameter values.  I dislike the
  ///   "nounification" approach to object-oriented programming, but
  ///   the deficiencies of C++ compilers require it.
  ///
  template< class MatrixType,
	    class RangeScalar = typename MatrixType::scalar_type,
	    class DomainScalar = typename MatrixType::scalar_type >
  class PolyBasisFactory {
  public:
    typedef MatrixType matrix_type;
    typedef DomainScalar domain_scalar_type;
    typedef RangeScalar range_scalar_type;
    typedef typename matrix_type::local_ordinal_type local_ordinal_type;
    typedef typename matrix_type::global_ordinal_type global_ordinal_type;
    typedef typename matrix_type::node_type node_type;
    typedef PolyBasis< domain_scalar_type, range_scalar_type, local_ordinal_type, global_ordinal_type, node_type > poly_basis_type;

    static Teuchos::RCP< poly_basis_type >
    makePolyBasis (const Teuchos::RCP< const matrix_type >& A,
		   const int maxPower,
		   const Teuchos::ParameterList& tuningHints);

    // Default tuning hints to pass to the nonmember constructor.
    //static Teuchos::ParameterList defaultHints();

    /// Like makePolyBasis above, but using default tuning hints.
    ///
    static Teuchos::RCP< poly_basis_type >
    makePolyBasis (const Teuchos::RCP< const matrix_type >& A,
		   const int maxPower);
  };


} // namespace Tpetra


#endif // __Tpetra_PolyBasis_hpp
