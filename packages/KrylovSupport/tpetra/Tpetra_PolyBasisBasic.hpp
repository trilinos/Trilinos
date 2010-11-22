#include "Tpetra_PolyBasis.hpp"
#include "Tpetra_Operator.hpp"

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace Tpetra {

  // Forward declaration
  template< class MatrixScalar, 
	    class LocalOrdinal = int, 
	    class GlobalOrdinal = LocalOrdinal, 
	    class Node = Kokkos::DefaultNode::DefaultNodeType, 
	    class RangeScalar = MatrixScalar,
	    class DomainScalar = MatrixScalar >
  class PolyBasisBasic;

  // PolyBasisFactory for Tpetra::Operator, that provides a sensible
  // default implementation.
  template< class MatrixScalar, 
	    class LocalOrdinal = int, 
	    class GlobalOrdinal = LocalOrdinal, 
	    class Node = Kokkos::DefaultNode::DefaultNodeType, 
	    class RangeScalar = MatrixScalar, 
	    class DomainScalar = MatrixScalar >
  class PolyBasisFactory< Tpetra::Operator< MatrixScalar, LocalOrdinal, GlobalOrdinal, Node >, RangeScalar, DomainScalar >
  {
  public:
    typedef Tpetra::Operator< MatrixScalar, LocalOrdinal, GlobalOrdinal, Node > matrix_type;
    typedef DomainScalar domain_scalar_type;
    typedef RangeScalar range_scalar_type;
    typedef LocalOrdinal local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node node_type;
    typedef PolyBasis< domain_scalar_type, range_scalar_type, local_ordinal_type, global_ordinal_type, node_type > poly_basis_type;

    static Teuchos::RCP< poly_basis_type >
    makePolyBasis (const Teuchos::RCP< const matrix_type >& A,
		   const int maxPower,
		   const Teuchos::ParameterList& tuningHints)
    {
      typedef PolyBasisBasic< MatrixScalar, LocalOrdinal, GlobalOrdinal, Node, RangeScalar, DomainScalar > impl_type;
      return Teuchos::rcp_implicit_cast< poly_basis_type > (Teuchos::rcp (new impl_type (A, maxPower)));
    }

    static Teuchos::RCP< poly_basis_type >
    makePolyBasis (const Teuchos::RCP< const matrix_type >& A,
		   const int maxPower)
    {
      Teuchos::ParameterList emptyHints;
      return makePolyBasis (A, maxPower, emptyHints);
    }
  };


  /// \class PolyBasisBasic
  /// \brief Default implementation of PolyBasis, for Tpetra::Operator.
  ///
  template< class MatrixScalar, 
	    class LocalOrdinal = int, 
	    class GlobalOrdinal = LocalOrdinal, 
	    class Node = Kokkos::DefaultNode::DefaultNodeType, 
	    class RangeScalar = MatrixScalar,
	    class DomainScalar = MatrixScalar >
  class PolyBasisBasic : 
    public PolyBasis< DomainScalar, RangeScalar, LocalOrdinal, GlobalOrdinal, Node >
  {
  public:
    typedef PolyBasis< DomainScalar, RangeScalar, LocalOrdinal, GlobalOrdinal, Node > base_type;
    typedef Tpetra::Operator< MatrixScalar, LocalOrdinal, GlobalOrdinal, Node > matrix_type;
    typedef typename base_type::domain_mv_type domain_mv_type;
    typedef typename base_type::range_mv_type range_mv_type;
    typedef typename base_type::range_magnitude_type range_magnitude_type;
    typedef typename base_type::shift_type shift_type;

    PolyBasisBasic (const Teuchos::RCP< const matrix_type >& A,
		    const int maxPower) :
      A_ (A), maxPower_ (maxPower) {}

    PolyBasisBasic (const Teuchos::RCP< const matrix_type >& A,
		    const int maxPower,
		    const Teuchos::ParameterList& tuningHints) :
      A_ (A), maxPower_ (maxPower) {}

    virtual ~PolyBasisBasic() {}

    virtual void 
    applyMonomial (const domain_mv_type& X,
		   range_mv_type& Y,
		   const int power) const;

    virtual void 
    applyModifiedNewton (const domain_mv_type& X,
			 range_mv_type& Y,
			 const std::vector< shift_type >& shifts,
			 const int power) const;
    
  private:
    Teuchos::RCP< const matrix_type > A_;
    int maxPower_;
  };

  template< class MatrixScalar, 
	    class LocalOrdinal = int, 
	    class GlobalOrdinal = LocalOrdinal, 
	    class Node = Kokkos::DefaultNode::DefaultNodeType,
	    class RangeScalar = MatrixScalar,
	    class DomainScalar = MatrixScalar >
  void 
  PolyBasisBasic< MatrixScalar, LocalOrdinal, GlobalOrdinal, Node, RangeScalar, DomainScalar >::
  applyMonomial (const domain_mv_type& X,
		 range_mv_type& Y,
		 const int power) const
  {
    using Teuchos::Range1D;
    using Teuchos::RCP;
    using Teuchos::ScalarTraits;
    
    applyCheck (X, Y, power, maxPower_);
    if (power == 0)
      return; // Do nothing
    
    Range1D Y_prvRange (0, X.getNumVectors() - 1);
    // Will be changed later.  Range1D only has op+= and not op+,
    // otherwise we would make Y_curRange local to the for loop.
    Range1D Y_curRange = Y_prvRange + power;
    
    for (int k = 0; k < power; ++k)
      {
	if (k == 0)
	  {
	    A_->apply (X, Y.subView (0, X.getNumVectors() - 1));
	  }
	else
	  {
	    RCP< range_mv_type > Y_prv = Y.subView (Y_prvRange);
	    RCP< range_mv_type > Y_cur = Y.subView (Y_curRange);
	    A_->apply (Y_prv, Y_cur);
	  }
	Y_prvRange += power;
	Y_curRange += power;
      }
  }


  template< class MatrixScalar, 
	    class LocalOrdinal = int, 
	    class GlobalOrdinal = LocalOrdinal, 
	    class Node = Kokkos::DefaultNode::DefaultNodeType,
	    class RangeScalar = MatrixScalar,
	    class DomainScalar = MatrixScalar >
  void 
  PolyBasisBasic< MatrixScalar, LocalOrdinal, GlobalOrdinal, Node, RangeScalar, DomainScalar >::
  applyRealNewton (const domain_mv_type& X,
		   range_mv_type& Y,
		   const std::vector< shift_type >& shifts,
		   const int power) const
  {
    using Teuchos::Range1D;
    using Teuchos::RCP;
    using Teuchos::ScalarTraits;

    applyCheck (X, Y, power, maxPower_);
    if (power == 0)
      return; // Do nothing

    Range1D Y_prvRange (0, X.getNumVectors() - 1);
    // Will be changed later.  Range1D only has op+= and not op+,
    // otherwise we would make Y_curRange local to the for loop.
    Range1D Y_curRange = Y_prvRange;
    // Will also be changed later.
    Range1D Y_nxtRange = Y_prvRange; 

    const range_magnitude_type ZERO = ScalarTraits< range_scalar_type >::zero();
    const range_magnitude_type ONE = ScalarTraits< range_scalar_type >::one();

    Range1D Y_prvRange (0, X.getNumVectors() - 1);
    Range1D Y_curRange (0, X.getNumVectors() - 1);
    int k = 0;
    while (k < shifts.size())
      {
	RCP< range_mv_type > Y_prv;
	if (k == 0)
	  Y_prv = X.subView (0, X.getNumVectors() - 1);
	else
	  {
	    Y_prvRange = Y_curRange -= power;
	    Y_prv = Y.subView (Y_prvRange);
	  }

	if (shifts[k].second == ZERO)
	  {
	    const range_magnitude_type Re_shift_k = shifts[k].first;
	    RCP< range_mv_type > Y_cur = Y.subView (Y_curRange);
	    A_->apply (Y_prv, Y_cur, Teuchos::NO_TRANS, Re_shift_k);
	    Y_curRange += power;
	    k = k + 1;
	  }
	else
	  {
	    // We should require shifts[k].second ==
	    // -shifts[k+1].second, but depending on how the shifts
	    // are computed, this may only be true approximately.
	    // We _assume_ that it is true, but don't test for it.
	    if (k+1 >= shifts.size() || shifts[k].second < 0 || shifts[k+1].second > 0)
	      {
		std::ostringstream os;
		os << "Complex shifts must occur in conjugate pairs, and the"
		  " member of each pair with positive imaginary part must "
		  "occur first";
		throw std::invalid_argument (os.str());
	      }

	    const range_magnitude_type Re_shift_k = shifts[k].first;
	    const range_magnitude_type Im_shift_k = shifts[k].second;
	    RCP< range_mv_type > Y_cur = Y.subView (Y_curRange);

	    Range1D Y_nxtRange = Y_curRange; // no op+, alas
	    Y_nxtRange += power;
	    RCP< range_mv_type > Y_nxt = Y.subView (Y_nxtRange);

	    // $Y_{\text{cur}} = (A - \Re{\theta_k} I) Y_{\text{prv}}$
	    A_->apply (Y_prv, Y_cur, Teuchos::NO_TRANS, Re_shift_k);
	    // $Y_{\text{nxt}} = (A - \Re{\theta_k} I) Y_{\text{cur}} - \Im{\theta_k}^2 Y_{\text{prv}}$
	    A_->apply (Y_cur, Y_nxt, Teuchos::NO_TRANS, Re_shift_k);
	    Y_nxt.update (Im_shift_k*Im_shift_k, Y_prv, ONE);
	    Y_curRange += power;
	    k = k + 2; // skip over $\theta_{k+1} = \overline{\theta_k}$
	  }
      }
  }

} // namespace Tpetra
