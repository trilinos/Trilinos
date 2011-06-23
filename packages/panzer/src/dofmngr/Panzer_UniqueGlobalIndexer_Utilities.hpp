#ifndef __Panzer_UniqueGlobalIndexer_Utilities_hpp__
#define __Panzer_UniqueGlobalIndexer_Utilities_hpp__

#include <map>
#include <string>

#include <boost/unordered_map.hpp>

#include <Kokkos_DefaultNode.hpp>
#include <Tpetra_Vector.hpp>

#include "Panzer_UniqueGlobalIndexer.hpp"

namespace panzer {

/** Construct a vector that contains a reduced set of field numbers.
  * The ordering is based on the ordering from <code>ugi.getOwnedAndSharedIndices()</code>.
  * The term "reduced" means that this processor must be able to fully determine the
  * field number for each global number. There are some cases where processors are split
  * across element blocks with differing fields that this is nontrivial. 
  *
  * \param[in] ugi Global indexer to use.
  * 
  * \returns Reduced vector containing the field numbers.
  *
  * \note The description and use of this function are equally confusing...
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
Teuchos::RCP<const Tpetra::Vector<int,std::size_t,GlobalOrdinalT,Node> >
buildGhostedFieldReducedVector(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi);

/** This function builds a vector that defines fields for each global unknown.
  * Notice that requires global communication and uses (underneath) the <code>Tpetra</code>
  * vector hence the required node type template parameter.
  *
  * This function returns a vector that serves as a map between Global indices ordered from
  * <code>ugi.getOwnedAndSharedIndices()</code> and the corresponding field number.
  *
  * \param[in] ugi Unique global indexer object that defines the ordering, global ids and field numbers.
  * \param[in] reducedVec Reduced field vector to use.  If none is passed it is compute by <code>buildGhostedFieldReducedVector</code>
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
Teuchos::RCP<const Tpetra::Vector<int,std::size_t,GlobalOrdinalT,Node> >
buildGhostedFieldVector(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi,
                        const Teuchos::RCP<const Tpetra::Vector<int,std::size_t,GlobalOrdinalT,Node> > & reducedVec=Teuchos::null);

/** This function builds a vector that defines fields for each global unknown.
  * Notice that requires global communication and uses (underneath) the <code>Tpetra</code>
  * vector hence the required node type template parameter.
  *
  * This function returns a vector that serves as a map between Global indices ordered from
  * <code>ugi.getOwnedAndSharedIndices()</code> and the corresponding field number.
  *
  * \param[in] ugi Unique global indexer object that defines the ordering, global ids and field numbers.
  * \param[out] fieldNumbers Field numbers ordered to match up with a call to <code>ugi.getOwnedAndSharedIndices()</code>.
  *                          Meaning that <code>fieldNumbers.size()</code> is the same length as the vector
  *                          build by a call to <code>getOwnedAndSharedIndices()</code>.
  * \param[in] reducedVec Reduced field vector to use.  If none is passed it is compute by <code>buildGhostedFieldReducedVector</code>
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void buildGhostedFieldVector(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi,
                             std::vector<int> & fieldNumbers,
                             const Teuchos::RCP<const Tpetra::Vector<int,std::size_t,GlobalOrdinalT,Node> > & reducedVec=Teuchos::null);

/** Convenience function default to the basic Kokkos node type.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
void buildGhostedFieldVector(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi,
                             std::vector<int> & fieldNumbers,
                             const Teuchos::RCP<const Tpetra::Vector<int,std::size_t,GlobalOrdinalT,Kokkos::DefaultNode::DefaultNodeType> > & reducedVec=Teuchos::null)
{ buildGhostedFieldVector<LocalOrdinalT,GlobalOrdinalT,Kokkos::DefaultNode::DefaultNodeType>(ugi,fieldNumbers,reducedVec); }

} // end namspace panzer

#include "Panzer_UniqueGlobalIndexer_UtilitiesT.hpp"

#endif
