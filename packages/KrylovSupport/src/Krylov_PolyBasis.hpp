#ifndef __Krylov_PolyBasis_hpp
#define __Krylov_PolyBasis_hpp

#include <Teuchos_ScalarTraits.hpp>
#include <utility> // std::pair

namespace Krylov {

  /// \class PolyBasis
  /// \brief Interface to the matrix powers kernel
  /// 
  /// This class defines the interface to the matrix powers kernel.
  /// The simplest form (monomial basis) of the matrix powers kernel
  /// takes a MultiVector X, and computes \fn$[A X, A^2 X, \dots, A^s
  /// X]\fn$, for some positive integer power \fn$s\fn$: thus, a basis
  /// of a (block) Krylov subspace.  (Other bases besides this
  /// monomial basis are supported; see references below for details
  /// on what that means.)  The point of the matrix powers kernel is
  /// to compute this basis with minimal communication.
  /// "Communication" includes all kinds of data movement, both
  /// between processors in parallel (either via message passing or a
  /// shared address space), and between levels of the memory
  /// hierarchy.
  ///
  /// This interface is meant to be used as an abstract base class for
  /// concrete implementations.  Unlike Krylov::OperatorTraits,
  /// implementations of subclasses of Krylov::PolyBasis are allowed
  /// to have state, and the methods are instance methods rather than
  /// class ("static" in C++ terms) methods.  That's why we don't call
  /// PolyBasis a "traits" class; traits classes have no state.
  ///
  /// Template parameters:
  /// - RangeScalar: Type of the elements of the output multivector
  /// - RangeMultiVector: Type of the output multivector
  /// - DomainMultiVector: Type of the input multivector
  ///
  /// References:
  ///
  /// For mathematical details on the matrix powers kernel, see:
  ///
  /// Mark Hoemmen, "Communication-avoiding Krylov subspace methods,"
  /// PhD thesis, University of California Berkeley EECS, 2010.
  ///
  /// For details on the efficient implementation of the matrix powers
  /// kernel in a shared-memory parallel environment with nonuniform
  /// memory access (NUMA), see
  ///
  /// Marghoob Mohiyuddin, Mark Hoemmen, James Demmel, and Kathy
  /// Yelick, "Minimizing Communication in Sparse Matrix Solvers," in
  /// Proceeedings of Supercomputing 2009.
  ///
  /// For details on the distributed-memory implementation of the
  /// matrix powers kernel, see
  ///
  /// James Demmel, Mark Hoemmen, Marghoob Mohiyuddin, and Katherine
  /// Yelick, "Avoiding Communication in Sparse Matrix Computations,"
  /// IEEE International Parallel and Distributed Processing Symposium
  /// (IPDPS), April 2008.
  /// 
  /// as well as the following University of California Berkeley EECS
  /// technical report:
  ///
  /// James Demmel, Mark Hoemmen, Marghoob Mohiyuddin, and Katherine
  /// Yelick, "Avoiding Communication in Computing Krylov Subspaces,"
  /// UCB/EECS-2007-123, October 2007.
  ///
  template< class RangeScalar,
	    class RangeMultiVector,
	    class DomainScalar,
	    class DomainMultiVector >
  class PolyBasis {
  public:
    typedef DomainScalar domain_scalar_type;
    typedef RangeScalar range_scalar_type;

    /// Type of the magnitude (absolute value) of a DomainScalar.
    /// If RangeScalar is complex, range_magnitude_type is real.
    typedef typename Teuchos::ScalarTraits< domain_scalar_type >::magnitudeType domain_magnitude_type;

    /// Type of the magnitude (absolute value) of a RangeScalar.
    /// If RangeScalar is complex, range_magnitude_type is real.
    typedef typename Teuchos::ScalarTraits< range_scalar_type >::magnitudeType range_magnitude_type;

    /// \brief Type of shifts for the (modified) Newton basis.
    ///
    /// We store possibly complex shifts for the Newton basis as pairs
    /// of (real, imaginary), where real and imaginary are each a
    /// range_magnitude_type.  This is because the shifts come from
    /// the eigenvalues of a nonsymmetric matrix, and so they might be
    /// complex-valued, even if the matrix itself is real.  Storing
    /// the shifts as pairs means that you can use the matrix powers
    /// kernel, even if you are compiling Trilinos without complex
    /// arithmetic support.  (Our implementations take care to do no
    /// complex arithmetic when all the Scalar types being used for
    /// the matrix and vectors are real.)
    ///
    /// \note Shifts are stored using the range's magnitude type,
    ///   because the typical case is that the range's scalar type has
    ///   at least as much (if not more) precision as the domain's
    ///   scalar type.  If the shifts have the domain's scalar type,
    ///   precision would be lost in the multiplication.
    typedef std::pair< range_magnitude_type, range_magnitude_type > shift_type;

    /// \brief Multivector input type.
    ///
    /// PolyBasis defines an interface for a mapping between
    /// MultiVectors.  The input and output MultiVector types may be
    /// the same, or may differ only in their Scalar types.
    typedef DomainMultiVector domain_mv_type;

    /// \brief Multivector output type.  
    /// 
    /// Should be the same as domain_mv_type (see above), except
    /// perhaps for a different Scalar type. The best reason to have
    /// different Scalar types for input and output is to improve
    /// precision of the matrix powers kernel by accumulating the
    /// output vectors in higher precision than the input.
    typedef RangeMultiVector range_mv_type;

    virtual ~PolyBasis () = 0;

    /// \brief Compute Y = [A*X, A^2 * X, ..., A^power * X]
    ///
    /// \param X [in]   Input multivector
    /// \param Y [out]  Output multivector:  
    ///   \fn$Y = [A*X, A^2 * X, ..., A^{\text{power}} * X]\fn$.
    ///   In Y, te columns of \fn$A^k X\fn$ are stored adjacently 
    ///   for all k.
    /// \param power [in] Number of times to "apply" the operator A.
    ///   If A is stored in a particular way, efficient
    ///   implementations will rearrange this computation in order to
    ///   avoid communication.
    virtual void 
    applyMonomial (const domain_mv_type& X,
		   range_mv_type& Y,
		   const int power) const = 0;

    /// Compute 'power' steps of the "modified Newton" basis.  If all
    /// Scalar types are real, this basis avoids complex arithmetic,
    /// even if the shifts for the Newton basis are complex.  If there
    /// are shifts with nonzero imaginary part, they must occur in
    /// consecutive complex conjugate pairs, and the first element in
    /// each pair must have positive imaginary part.
    ///
    /// \param X [in]
    /// \param Y [out]
    /// \param power [in] Number of times to "apply" the operator A
    ///
    /// \note See Mark Hoemmen's PhD dissertation, Chapter 7, Section
    /// 3.2 for a definition of the "modified Newton" basis.
    ///
    /// \note This is the same thing as the regular Newton basis, if
    ///   both Scalar types and all the shifts are real.
    virtual void 
    applyModifiedNewton (const domain_mv_type& X,
			 range_mv_type& Y,
			 const std::vector< shift_type >& shifts,
			 const int power) const = 0;
  };

} // namespace Krylov


#endif // __Krylov_PolyBasis_hpp
