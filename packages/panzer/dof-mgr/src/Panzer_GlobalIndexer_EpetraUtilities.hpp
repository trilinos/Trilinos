// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_GlobalIndexer_EpetraUtilities_decl_hpp__
#define __Panzer_GlobalIndexer_EpetraUtilities_decl_hpp__

#include <map>
#include <string>

#include <unordered_map>

#include <Panzer_NodeType.hpp>
#include <Epetra_Vector.h>
#include <Epetra_IntVector.h>
#include <Epetra_MultiVector.h>

#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_FieldPattern.hpp"

namespace panzer {


/** Construct a vector that contains a reduced set of field numbers.
  * The ordering is based on the ordering from <code>ugi.getOwnedAndGhostedIndices()</code>.
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
Teuchos::RCP<Epetra_IntVector>
buildGhostedFieldReducedVectorEpetra(const GlobalIndexer & ugi);

/** This function builds a vector that defines fields for each global unknown.
  * Notice that requires global communication and uses (underneath) the <code>Tpetra</code>
  * vector hence the required node type template parameter.
  *
  * This function returns a vector that serves as a map between Global indices ordered from
  * <code>ugi.getOwnedAndGhostedIndices()</code> and the corresponding field number.
  *
  * \param[in] ugi Unique global indexer object that defines the ordering, global ids and field numbers.
  * \param[in] reducedVec Reduced field vector to use.  If none is passed it is compute by <code>buildGhostedFieldReducedVector</code>
  */
Teuchos::RCP<const Epetra_IntVector>
buildGhostedFieldVectorEpetra(const GlobalIndexer & ugi,
                              const Teuchos::RCP<const Epetra_IntVector> & reducedVec=Teuchos::null);

/** This function builds a vector that defines fields for each global unknown.
  * Notice that requires global communication and uses (underneath) the <code>Tpetra</code>
  * vector hence the required node type template parameter.
  *
  * This function returns a vector that serves as a map between Global indices ordered from
  * <code>ugi.getOwnedAndGhostedIndices()</code> and the corresponding field number.
  *
  * \param[in] ugi Unique global indexer object that defines the ordering, global ids and field numbers.
  * \param[out] fieldNumbers Field numbers ordered to match up with a call to <code>ugi.getOwnedAndGhostedIndices()</code>.
  *                          Meaning that <code>fieldNumbers.size()</code> is the same length as the vector
  *                          build by a call to <code>getOwnedAndGhostedIndices()</code>.
  * \param[in] reducedVec Reduced field vector to use.  If none is passed it is compute by <code>buildGhostedFieldReducedVector</code>
  */
void buildGhostedFieldVectorEpetra(const GlobalIndexer & ugi,
                                   std::vector<int> & fieldNumbers,
                                   const Teuchos::RCP<const Epetra_IntVector> & reducedVec=Teuchos::null);

/** Build a reduced data vector using the reduced field vector. Here reduced is meant in the
  * exact same context as for the field vectors.
  *
  * \param[in] fieldName Name of field data should be ordered as
  * \param[in] ugi Unique global indexer object that defines the ordering, global ids and field numbers.
  * \param[in] reducedFieldVec A reduced vector containing the correctly ordered field numbers.
  *                            Likely computed by <code>buildGhostedFieldReducedVector</code>.
  * \param[in] data Data to put in vector
  * \param[out] A Tpetra vector containing the data in the reduced format.  This is
  *             now available for an import to construct the true ghosted data vector. This
  *             map must match the reducedFieldVec map.
  */
template <typename ScalarT,typename ArrayT>
void updateGhostedDataReducedVectorEpetra(const std::string & fieldName,const std::string blockId,
                                          const GlobalIndexer & ugi,
                                          const ArrayT & data,
                                          Epetra_MultiVector & dataVector);

/** Construct a map that only uses a certain field.
  *
  * \param[in] fieldNum Field ID to use to build the map
  * \param[in] fieldVector An listing of the fields. For instance [0,1,2,0,1,2,0,1,2...]
  *                        for a three dimensional interlaced vector field.
  *
  * \returns A vector that contains the indices of the field requested. Again for the
  *          three dimensional vector field if fieldNum==1, it would be [1,4,7,10,13,...].
  */
Teuchos::RCP<const Epetra_BlockMap>
getFieldMapEpetra(int fieldNum,const Epetra_IntVector & fieldVector);


/** This class assists in mapping arrays of field data to field vectors.
  */
class ArrayToFieldVectorEpetra {
public:
   /** Construct information for the unique global indexer. Notice that this
     * requires global communication.
     */
   ArrayToFieldVectorEpetra(const Teuchos::RCP<const GlobalIndexer> & ugi);

   /** Get a Tpetra vector containing the data ordered according to 
     * the ordering from <code>UGI::getOwnedAndGhostedIndices</code>.
     *
     * \param[in] fieldName Name of field this data is from
     * \param[in] data Array of data
     *
     * \returns Returns a vector populated with the data. This vector
     *          is related to the <code>UGI::getOwnedAndGhostedIndices</code>.
     */
   template <typename ScalarT,typename ArrayT>
   Teuchos::RCP<Epetra_MultiVector>
   getGhostedDataVector(const std::string & fieldName,const std::map<std::string,ArrayT> & data) const;

   /** Get a Tpetra vector containing the data ordered according to 
     * the ordering from <code>UGI::getOwnedIndices</code>.
     *
     * \param[in] fieldName Name of field this data is from
     * \param[in] data Array of data
     *
     * \returns Returns a vector populated with the data. This vector
     *          is related to the <code>UGI::getOwnedIndices</code>.
     */
   template <typename ScalarT,typename ArrayT>
   Teuchos::RCP<Epetra_MultiVector>
   getDataVector(const std::string & fieldName,const std::map<std::string,ArrayT> & data) const;

   /** Build a map that contains only global IDs associated with a particular field.
     * This serves to go from a unique vector of all fields, to a vector
     * containing the uniquely owned global ids for a single field.
     */
   Teuchos::RCP<const Epetra_Map>
   getFieldMap(const std::string & fieldName) const;

   /** Build a map that contains only global IDs associated with a particular field.
     * This serves to go from a unique vector of all fields, to a vector
     * containing the uniquely owned global ids for a single field.
     */
   Teuchos::RCP<const Epetra_Map>
   getFieldMap(int fieldNum) const;

protected:
   typedef Epetra_IntVector IntVector;
   typedef Epetra_MultiVector MultiVector;
   typedef Epetra_BlockMap Map;

   //! build unghosted field vector from ghosted field vector
   void buildFieldVector(const Epetra_IntVector & source) const;

   //! DOF mapping
   Teuchos::RCP<const GlobalIndexer> ugi_;

   Teuchos::RCP<const IntVector> gh_reducedFieldVector_; //! ghosted reduced field vector
   Teuchos::RCP<const IntVector> gh_fieldVector_;        //! ghosted field vector

   mutable std::map<int,Teuchos::RCP<const Map>> gh_reducedFieldMaps_; //! Maps for each field (as needed)
   mutable std::map<int,Teuchos::RCP<const Map>> gh_fieldMaps_;        //! Maps for each field (as needed)

   mutable Teuchos::RCP<const IntVector> fieldVector_;               //! (unghosted) field vector (as needed)
   mutable std::map<int,Teuchos::RCP<const Map>> fieldMaps_;        //! Maps for each field (as needed)

private:
   // hide some constructors
   ArrayToFieldVectorEpetra();
   ArrayToFieldVectorEpetra(const ArrayToFieldVectorEpetra &);
};


} // end namspace panzer

#include "Panzer_GlobalIndexer_EpetraUtilities_impl.hpp"

#endif
