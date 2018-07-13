// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_UniqueGlobalIndexer_Utilities_decl_hpp__
#define __Panzer_UniqueGlobalIndexer_Utilities_decl_hpp__

#include <map>
#include <string>

#include <unordered_map>

#include <Panzer_NodeType.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>

#include "Panzer_UniqueGlobalIndexer.hpp"
#include "Panzer_FieldPattern.hpp"

namespace panzer {

/** Convert a nonconst to a constat vector. This works around an RCP issue where the compiler
  * gets confused by const.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > > 
nc2c_vector(const std::vector<Teuchos::RCP<UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > > & ugis);

/** Get the block associated with a particular field. This is an exhaustive (e.g. expensive)
  * search. This returns the first found index into the <code>ugis</code> that contains the 
  * required field.
  *
  * \param[in] fieldName The field to look for.
  * \param[in] ugis The global indexers to look in.
  *
  * \returns The index that this field is in. If this returns -1, no field was found.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
int getFieldBlock(const std::string & fieldName,
                  const std::vector<Teuchos::RCP<UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > > & ugis);

/** Get the block associated with a particular field. This is an exhaustive (e.g. expensive)
  * search. This returns the first found index into the <code>ugis</code> that contains the 
  * required field.
  *
  * \param[in] fieldName The field to look for.
  * \param[in] ugis The global indexers to look in.
  *
  * \returns The index that this field is in. If this returns -1, no field was found.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
int getFieldBlock(const std::string & fieldName,
                  const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > > & ugis);

/** Compute the offsets for the global indexer for a particular block id. This
  * is useful for unknown numbering into a block format. The local element matrix
  * will be logically ordered like the monolithic operator. These offsets are the
  * local start and end of a particular block. For instance if you are intersted in
  * block <code>i</code> then the start index would be at <code>blockOffsets[i]</code>
  * and the end would be at <code>blockOffsets[i+1]</code> Note that this requires a
  * sentinnel in the block offsets size. (Note that the use of element block and block
  * are confusing but different in this context)
  *
  * \param[in] blockId Element block name that to do this for.
  * \param[in] ugis Unique global indexers containing fields to be blocked
  * \param[out] blockOffsets Result vector with length <code>ugis.size()+1</code>.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
void computeBlockOffsets(const std::string & blockId,
                         const std::vector<Teuchos::RCP<UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > > & ugis,
                         std::vector<int> & blockOffsets);

/** Compute the offsets for the global indexer for a particular block id. This
  * is useful for unknown numbering into a block format. The local element matrix
  * will be logically ordered like the monolithic operator. These offsets are the
  * local start and end of a particular block. For instance if you are intersted in
  * block <code>i</code> then the start index would be at <code>blockOffsets[i]</code>
  * and the end would be at <code>blockOffsets[i+1]</code> Note that this requires a
  * sentinnel in the block offsets size. (Note that the use of element block and block
  * are confusing but different in this context)
  *
  * \param[in] blockId Element block name that to do this for.
  * \param[in] ugis Unique global indexers containing fields to be blocked
  * \param[out] blockOffsets Result vector with length <code>ugis.size()+1</code>.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
void computeBlockOffsets(const std::string & blockId,
                         const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > > & ugis,
                         std::vector<int> & blockOffsets);

/** Print out unique global indexer load balancing information. This includes
  * the minimum unknown, maximum unknown count, mean unknown and standard deviation of
  * the unknowns for both owned and owned and ghosted.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
std::string printUGILoadBalancingInformation(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi);

/** Print all GIDs and their associated elements to the screen. Note if a FancyOStream is used
  * the correct prefixes this method will label the processors as well. This can print out an 
  * extreme amount of information so it is only useful for debugging.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
void printMeshTopology(std::ostream & os,const panzer::UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi);

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
template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
Teuchos::RCP<Tpetra::Vector<int,int,GlobalOrdinalT,Node> >
buildGhostedFieldReducedVector(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi);

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
template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
Teuchos::RCP<const Tpetra::Vector<int,int,GlobalOrdinalT,Node> >
buildGhostedFieldVector(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi,
                        const Teuchos::RCP<const Tpetra::Vector<int,int,GlobalOrdinalT,Node> > & reducedVec=Teuchos::null);

/** Convenience function default to the basic Kokkos node type.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
Teuchos::RCP<const Tpetra::Vector<int,int,GlobalOrdinalT,panzer::TpetraNodeType> >
buildGhostedFieldVector(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi,
                        const Teuchos::RCP<const Tpetra::Vector<int,int,GlobalOrdinalT,panzer::TpetraNodeType> > & reducedVec=Teuchos::null)
{ return buildGhostedFieldVector<LocalOrdinalT,GlobalOrdinalT,panzer::TpetraNodeType>(ugi,reducedVec); }

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
template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void buildGhostedFieldVector(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi,
                             std::vector<int> & fieldNumbers,
                             const Teuchos::RCP<const Tpetra::Vector<int,int,GlobalOrdinalT,Node> > & reducedVec=Teuchos::null);

/** Convenience function default to the basic Kokkos node type.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT>
void buildGhostedFieldVector(const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi,
                             std::vector<int> & fieldNumbers,
                             const Teuchos::RCP<const Tpetra::Vector<int,int,GlobalOrdinalT,panzer::TpetraNodeType> > & reducedVec=Teuchos::null)
{ buildGhostedFieldVector<LocalOrdinalT,GlobalOrdinalT,panzer::TpetraNodeType>(ugi,fieldNumbers,reducedVec); }

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
template <typename ScalarT,typename ArrayT,typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
void updateGhostedDataReducedVector(const std::string & fieldName,const std::string blockId,
                                    const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> & ugi,
                                    const ArrayT & data,
                                    Tpetra::MultiVector<ScalarT,int,GlobalOrdinalT,Node> & dataVector);

/** Construct a map that only uses a certain field.
  *
  * \param[in] fieldNum Field ID to use to build the map
  * \param[in] fieldVector An listing of the fields. For instance [0,1,2,0,1,2,0,1,2...]
  *                        for a three dimensional interlaced vector field.
  *                        
  * \returns A vector that contains the indices of the field requested. Again for the
  *          three dimensional vector field if fieldNum==1, it would be [1,4,7,10,13,...].
  */
template <typename GlobalOrdinalT,typename Node>
Teuchos::RCP<const Tpetra::Map<int,GlobalOrdinalT,Node> >
getFieldMap(int fieldNum,const Tpetra::Vector<int,int,GlobalOrdinalT,Node> & fieldVector);


namespace orientation_helpers {

/** For a given field pattern compute the offsets that give 
  * the beginning and end dimension-0-subcell index for each edge.
  *
  * \param[in] pattern Pattern specifying the layout of IDs. Note that
  *                    this pattern must have dimension-0-subcell indices.
  * \param[in,out] edgeIndices Empty vector that on exit will have a pair
  *                            containing start and end indices for each
  *                            edge in a cell.
  */
void computePatternEdgeIndices(const FieldPattern & pattern,std::vector<std::pair<int,int> > & edgeIndices);

/** This function computes the edge orientation for a cell given its global topology.
  * It is most often called in conjunction with <code>computePatternEdgeIndices</code>. The
  * general model is to call <code>computePatternEdgeIndices</code> once for a given topology
  * field pattern. These edge indices are used with the topology vector of GIDs (which is laid out 
  * using the topology field pattern) for particular to define the orientation for that element.
  * The layout of the orientation vector is defined to satisfy yet another field pattern whose
  * structure defines the global unknowns for that element. This function can then be called
  * repeatedly for each element that satisfies the topology used in the 
  * <code>computePatternEdgeIndices</code> call.
  *
  * \param[in] topEdgeIndices Computed by <code>computePatternEdgeIndices</code>
  * \param[in] topology Global topology of this element satisfying field pattern used to construct
  *                     the <code>edgeIndices</code> vector.
  * \param[in] fieldPattern Field pattern used to define the orientation layout
  * \param[in,out] orientation Orientation vector satisfying the field pattern layout
  */
template <typename GlobalOrdinalT>
void computeCellEdgeOrientations(const std::vector<std::pair<int,int> > & topEdgeIndices,
                                 const std::vector<GlobalOrdinalT> & topology,
                                 const FieldPattern & fieldPattern, 
                                 std::vector<signed char> & orientation);

/** For a given field pattern compute the offsets that give 
  * the dimension-0-subcell indices for each face. Note that the assumption is
  * made that the node ordering returned by Shards for each face is counter-clockwise.
  * This is how the determination of inward or outward facing normals is made.
  *
  * \param[in] pattern Pattern specifying the layout of IDs. Note that
  *                    this pattern must have dimension-0-subcell indices.
  * \param[in,out] faceIndices Empty vector that on exit will have a vector
  *                            containing start and end indices for each
  *                            face in a cell.
  *
  * \note In 2D a "Face" is defined to be the element (not the edge).
  */
void computePatternFaceIndices(const FieldPattern & pattern,std::vector<std::vector<int> > & faceIndices);

/** This function computes the face orientation for a cell given its global topology.
  * It is most often called in conjunction with <code>computePatternFaceIndices</code>. The
  * general model is to call <code>computePatternFaceIndices</code> once for a given topology
  * field pattern. These face indices are used with the topology vector of GIDs (which is laid out 
  * using the topology field pattern) to define the orientation for that element.
  * The layout of the orientation vector is defined to satisfy yet another field pattern whose
  * structure defines the global unknowns for that element. This function can then be called
  * repeatedly for each element that satisfies the topology used in the 
  * <code>computePatternFaceIndices</code> call.
  *
  * \param[in] topFaceIndices Computed by <code>computePatternFaceIndices</code>
  * \param[in] topology Global topology of this element satisfying field pattern used to construct
  *                     the <code>edgeIndices</code> vector.
  * \param[in] fieldPattern Field pattern used to define the orientation layout
  * \param[in,out] orientation Orientation vector satisfying the field pattern layout
  */
template <typename GlobalOrdinalT>
void computeCellFaceOrientations(const std::vector<std::pair<int,int> > & topEdgeIndices,
                                 const std::vector<GlobalOrdinalT> & topology,
                                 const FieldPattern & fieldPattern, 
                                 std::vector<signed char> & orientation);

} // end orientations_helpers


/** This class assists in mapping arrays of field data to field vectors.
  */
template <typename LocalOrdinalT,typename GlobalOrdinalT,typename Node>
class ArrayToFieldVector {
public:
   /** Construct information for the unique global indexer. Notice that this
     * requires global communication.
     */
   ArrayToFieldVector(const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > & ugi);

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
   Teuchos::RCP<Tpetra::MultiVector<ScalarT,int,GlobalOrdinalT,Node> >
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
   Teuchos::RCP<Tpetra::MultiVector<ScalarT,int,GlobalOrdinalT,Node> >
   getDataVector(const std::string & fieldName,const std::map<std::string,ArrayT> & data) const;

   /** Build a map that contains only global IDs associated with a particular field.
     * This serves to go from a unique vector of all fields, to a vector
     * containing the uniquely owned global ids for a single field.
     */
   Teuchos::RCP<const Tpetra::Map<int,GlobalOrdinalT,Node> >
   getFieldMap(const std::string & fieldName) const;

   /** Build a map that contains only global IDs associated with a particular field.
     * This serves to go from a unique vector of all fields, to a vector
     * containing the uniquely owned global ids for a single field.
     */
   Teuchos::RCP<const Tpetra::Map<int,GlobalOrdinalT,Node> >
   getFieldMap(int fieldNum) const;

protected:
   typedef Tpetra::Vector<int,int,GlobalOrdinalT,Node> IntVector;
   typedef Tpetra::Map<int,GlobalOrdinalT,Node> Map;

   //! build unghosted field vector from ghosted field vector
   void buildFieldVector(const Tpetra::Vector<int,int,GlobalOrdinalT,Node> & source) const;

   //! DOF mapping
   Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,GlobalOrdinalT> > ugi_;

   Teuchos::RCP<const IntVector> gh_reducedFieldVector_; //! ghosted reduced field vector
   Teuchos::RCP<const IntVector> gh_fieldVector_;        //! ghosted field vector

   mutable std::map<int,Teuchos::RCP<const Map> > gh_reducedFieldMaps_; //! Maps for each field (as needed)
   mutable std::map<int,Teuchos::RCP<const Map> > gh_fieldMaps_;        //! Maps for each field (as needed)

   mutable Teuchos::RCP<const IntVector> fieldVector_;               //! (unghosted) field vector (as needed)
   mutable std::map<int,Teuchos::RCP<const Map> > fieldMaps_;        //! Maps for each field (as needed)

private:
   // hide some constructors
   ArrayToFieldVector();
   ArrayToFieldVector(const ArrayToFieldVector &);
};


} // end namspace panzer

#include "Panzer_UniqueGlobalIndexer_Utilities_impl.hpp"

#endif
