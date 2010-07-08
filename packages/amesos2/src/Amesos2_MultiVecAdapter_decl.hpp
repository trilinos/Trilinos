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

#include <Teuchos_RCP.hpp>

namespace Amesos {

// TODO: Update method requirement list

/** \brief A templated MultiVec class adapter for Amesos2
 *
 * Specializations of this templated class provide a unified interface
 * to MultiVec types for Amesos2.  Any specializations are expected to
 * implement the following methods:
 *
 * <br><b>Implementation Requirements:</B>
 * <ul>
 * <li>Default constructor
 * \code MultiVecAdapter<MultiVecType>(); \endcode
 * </li>
 *
 * <li>Wrapper constructor
 * \code MultiVecAdapter<MultiVecType>(MultiVecType* const mat); \endcode
 * </li>
 *
 *<li>Copy constructor
 * \code MultiVecAdapter<MultiVecType>(MultiVecAdapter<MultiVecType>& const); \endcode
 * </li>
 *
 * <li>MultiVec scaling operation.  Scales each element by \c alpha
 * \code void scale(const scalar_type alpha); \endcode
 * </li>
 *
 * <li>MultiVec add operation.  Replace each element in \c this with \f$alpha*this + beta*B\f$
 * \code void add(const scalar_type beta, MultiVecAdapter<MultiVecType> B, const scalar_type alpha); \endcode
 * </li>
 *
 * <li>Method to get locality of multivec, either globally or locally indexed
 * \code bool isLocal() \endcode
 * </li>
 *
 * <li>Method to get multivec communicator
 * \code virtual const Teuchos::RCP<const Teuchos::Comm<int> >& getComm() const; \endcode
 * </li>
 *
 * <li>Methods to get the local and global length of vectors and number
 * of vectors
 * \code local_ordinal_type  getLocalLength();
 *       local_ordinal_type  getLocalNumVectors();
 *       global_ordinal_type getGlobalLength();
 *	global_ordinal_type getGlobalNumVectors(); \endcode
 * </li>
 *
 * <li>Method to access number of nonzero entries
 * \code local_ordinal_type  getLocalNNZ();
 *       global_ordinal_type getGlobalNNZ(); \endcode
 * </li>
 *
 * <li>Method to get the maximum number of nonzeros in all rows.
 * \code global_ordinal_type getMaxNNZ(); \endcode
 * </li>
 *
 * <li>Map methods.  We should make this method more general, in that we
 * should not be constraining ourselves to a particular Map object, and I
 * don't believe they all inherit from the same base class, and they all have
 * different interfaces.  Is there some way that we can get the functionality
 * from them, but make the interface more abstract?
 * \code Teuchos::RCP<Map> getRowMap();
 *       Teuchos::RCP<Map> getRangeMap();
 *	Teuchos::RCP<Map> getDomainMap(); \endcode
 * </li>
 *
 * <li>Row access methods
 * \code void getLocalRowCopy(local_ordinal_type localRow,
 *                            const Teuchos::ArrayView<local_ordinal_type> &indices,
 *                            const Teuchos::ArrayView<scalar_type> &values,
 *                            size_t &numEntries) const;
 * void getGlobalRowCopy(global_ordinal_type globalRow,
 *                       const Teuchos::ArrayView<global_ordinal_type> &indices,
 *                       const Teuchos::ArrayView<scalar_type> &values,
 *                       size_t &numEntries) const;
 * \endcode
 * </li>
 *
 * <li>Data access methods
 * \code void get1dCopy(scalar_type* A, int lda) const;
 *       void get1dView(scalar_type** A, int* lda) const;
 *       void get2dCopy(scalar_type** A, int lda) const;
 *       void get2dView(scalar_type*** A, int* lda) const;
 * \endcode
 * </li>
 *
 * <li>Method to convert MultiVec to serial.  Amesos2 itself should not have
 * to deal with the Maps (Epetra_Map, Tpetra::Map, etc) for this to work, the
 * method should simply return a reference to a serial *copy* of the original
 * multivector.
 * 
 * \code RCP<multivec_type>& serial(); \endcode
 * </li>
 *
 * <li>Method to revert MultiVec back to distributed.  Should reverse the
 * changes made to any method which localizes the MultiVec.
 *
 * <li>Get a description of this adapter
 * \code std::string description(); \endcode
 * </li>
 * 
 * <li>Print the multivec to the \c os output stream
 * \code void describe(Teuchos::FancyOStream& os) const; \endcode
 * </li>
 */
template <class MultiVecType>
struct MultiVecAdapter {

  /** \brief Default constuctor throws compiler error
   *
   * Without specialization, the Amesos2::MultiVecAdapter simply throws a
   * compiler error upon instantiation.
   */
  MultiVecAdapter()
    {
      MultiVecType::multivec_adapter_missing_for_this_type();
    }

  MultiVecAdapter(const Teuchos::RCP<MultiVecType>& vector)
    {
      MultiVecType::multivec_adapter_missing_for_this_type();
    }
};


} // end namespace Amesos

#endif  // AMESOS2_MULTIVEC_ADAPTER_DECL_HPP
