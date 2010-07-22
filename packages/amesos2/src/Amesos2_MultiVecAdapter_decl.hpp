/**
  \file   Amesos2_MultiVecAdapter.hpp
  \author Eric T Bavier <etbavier@sandia.gov>
  \date   Wed Feb 10 14:53:29 2010

  \brief  A templated adapter/wrapper class for Trilinos Multivector
          type classes.  Provides the functions necessary for Amesos2
          to function.  Other wrapper methods may of course be added
          for special cases.
*/

#ifndef AMESOS2_MULTIVEC_ADAPTER_DECL_HPP
#define AMESOS2_MULTIVEC_ADAPTER_DECL_HPP


namespace Amesos {


/**
 * \brief A templated MultiVector class adapter for Amesos2
 *
 * Specializations of this templated class provide a unified interface
 * to MultiVector types for Amesos2.  Any specializations are expected to
 * implement the following methods:
 *
 * <br>
 * <b>Implementation Requirements:</b>
 * <ul>
 * <li>Default constructor.
 *
 * \code
 * MultiVecAdapter<MultiVecType>();
 * \endcode
 * </li>
 *
 * <li>Wrapper constructor
 * \code MultiVecAdapter<MultiVecType>(const Teuchos::RCP<MultiVecType>& mat); \endcode
 * </li>
 *
 * <li>Copy constructor.
 *
 * \code
 * MultiVecAdapter<MultiVecType>(const MultiVecAdapter<MultiVecType>& const);
 * \endcode
 * </li>
 *
 * <li> MultiVec scaling operation.  Scales each element by \c alpha .
 *
 * \code
 * MultiVecAdapter<multivec_type>& scale(const scalar_type alpha);
 * \endcode
 * </li>
 *
 * <li>MultiVec update operation.  Replace each element in \c this with
 * \f$alpha*this + beta*B\f$.
 *
 * \code
 * MultiVecAdapter<multivec_type>&
 * add(const scalar_type beta, MultiVecAdapter<multivec_type>& B, const scalar_type alpha);
 * \endcode
 * </li>
 *
 * <li> Method to get locality of multivec, either globally or locally indexed.
 *
 * \code
 * bool isLocal() const;
 * \endcode
 * </li>
 *
 * <li> Method to get multi-vector communicator.
 *
 * \code
 * const Teuchos::RCP<const Teuchos::Comm<int> >& getComm() const;
 * \endcode
 * </li>
 *
 * <li> Methods to get the local and global length of vectors and number of
 * vectors.
 *
 * \code
 * size_t getLocalLength();
 *
 * size_t  getLocalNumVectors();
 *
 * global_size_type getGlobalLength();
 *
 * global_size_type getGlobalNumVectors();
 * \endcode
 * </li>
 *
 * <li> Get access to multi-vector stride
 *
 * \code
 * size_t getStride() const;
 * \endcode
 * </li>
 *
 * <li> Ask whether this multi-vector has constant stride between vectors on
 * this node.
 *
 * \code
 * bool isConstantStride() const;
 * \endcode
 * </li>
 *
 * <li> Vector access methods.
 *
 * \code
 * Teuchos::RCP<const Tpetra::Vector<Scalar,LO,GO,Node> >
 * getVector( size_t j ) const;
 *
 * Teuchos::RCP<Tpetra::Vector<Scalar,LO,GO,Node> >
 * getVectorNonConst( size_t j );
 * \endcode
 * </li>
 *
 * <li> Data access methods
 *
 * \code
 * void get1dCopy( const Teuchos::ArrayView<scalar_type>& A, size_t lda) const;
 *
 * Teuchos::ArrayRCP<scalar_type> get1dViewNonConst( bool local = false ) const;
 *
 * void get2dCopy( Teuchos::ArrayView<const Teuchos::ArrayView<scalar_type> > A ) const;
 *
 * Teuchos::ArrayRCP<Teuchos::ArrayRCP<scalar_type> > get2dViewNonConst() const;
 * \endcode
 * </li>
 *
 * <li> Method to export any local changes to this multi-vector back to the
 * global space.
 *
 * \code
 * void globalize();
 * \endcode
 * </li>
 *
 * <li> Method to export an array of new values into the global multi-vector.
 *
 * \code
 * template<typename Value_t>
 * void globalize( const Teuchos::ArrayView<Value_t>& newVals )
 * \endcode
 * </li>
 *
 * <li> Method to revert MultiVec back to distributed.  Should reverse the
 * changes made to any method which localizes the MultiVec.
 *
 * <li> Get a description of this adapter.
 *
 * \code
 * std::string description() const;
 * \endcode
 * </li>
 *
 * <li> Print the multivec to the \c os output stream.
 *
 * \code
 * void describe(
 *   Teuchos::FancyOStream& os,
 *   const Teuchos::EVerbosityLevel verblevel) const;
 * \endcode
 * </li>
 *
 * <li> A \c private method to localize the multi-vector.  Solver interfaces
 * themselves should not have to concern themselves with the details of how
 * this localization process is performed.  The method should create and/or
 * initialize private local variables as necessary.
 *
 * \code
 * void localize();
 * \endcode
 * </li>
 *
 */
template <class MultiVecType>
struct MultiVecAdapter {

  // /**
  //  * \brief Default constructor throws compiler error
  //  *
  //  * Without specialization, the Amesos2::MultiVecAdapter simply throws a
  //  * compiler error upon instantiation.
  //  */
  // MultiVecAdapter()
  //   {
  //     MultiVecType::multivec_adapter_missing_for_this_type();
  //   }

  // MultiVecAdapter(const Teuchos::RCP<MultiVecType>& vector)
  //   {
  //     MultiVecType::multivec_adapter_missing_for_this_type();
  //   }
};


} // end namespace Amesos

#endif  // AMESOS2_MULTIVEC_ADAPTER_DECL_HPP
