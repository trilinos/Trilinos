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

#ifndef __Panzer_TpetraVector_ReadOnly_GlobalEvaluationData_hpp__
#define __Panzer_TpetraVector_ReadOnly_GlobalEvaluationData_hpp__

#include "Tpetra_Import.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_Map.hpp"

#include "Teuchos_RCP.hpp"

#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"

#include "Panzer_NodeType.hpp"
#include "Panzer_ReadOnlyVector_GlobalEvaluationData.hpp"

namespace panzer {

/** This class provides a boundary exchange communication mechanism for vectors.
  * Not this provides a "read only" (RO) interface for parameters (so vectors are write protected).
  */
template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,
          typename NodeT=panzer::TpetraNodeType>
class TpetraVector_ReadOnly_GlobalEvaluationData : public ReadOnlyVector_GlobalEvaluationData {
public:
   typedef Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> VectorType;
   typedef Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> MapType;
   typedef Tpetra::Import<LocalOrdinalT,GlobalOrdinalT,NodeT> ImportType;

   //! Default constructor
   TpetraVector_ReadOnly_GlobalEvaluationData()
      : isInitialized_(false) { }

   TpetraVector_ReadOnly_GlobalEvaluationData(const TpetraVector_ReadOnly_GlobalEvaluationData & src)
      : isInitialized_(false) { initialize(src.importer_, src.ghostedMap_, src.ownedMap_); }

   /** Initialize this object with some Tpetra communication objects. This method
     * must be called before an object of this type can be used.
     *
     * \param[in] importer Importer for doing communication from the owned 
     *                     to the ghosted vector.
     * \param[in] ghostedMap Map describing the ghosted vector.
     * \param[in] ownedMap Map describing the ghosted vector.
     */
   TpetraVector_ReadOnly_GlobalEvaluationData(const Teuchos::RCP<const ImportType>& importer,
                                              const Teuchos::RCP<const MapType>&    ghostedMap,
                                              const Teuchos::RCP<const MapType>&    ownedMap)
      : isInitialized_(false) { initialize(importer, ghostedMap, ownedMap); }

   /** Choose a few GIDs and instead of zeroing them out in the ghosted vector set
     * them to a specified value. Note that this is only useful for GIDs in the
     * ghosted map that are not in the owned map.
     *
     * This must be called before initialize. Also note that no attempt to synchronize
     * these values a crossed processor is made. So its up to the user to be consistent.
     */
   void useConstantValues(const std::vector<GlobalOrdinalT> & indices,double value);

   /** Initialize this object with some Tpetra communication objects. This method
     * must be called before an object of this type can be used.
     *
     * \param[in] importer Importer for doing communication from the owned 
     *                     to the ghosted vector.
     * \param[in] ghostedMap Map describing the ghosted vector.
     * \param[in] ownedMap Map describing the ghosted vector.
     */
   void initialize(const Teuchos::RCP<const ImportType>& importer,
                   const Teuchos::RCP<const MapType>&    ghostedMap,
                   const Teuchos::RCP<const MapType>&    ownedMap);

   /** For this class, this method does the halo exchange for the 
     * vector.
     */
   virtual void globalToGhost(int mem);

   //! Clear out the ghosted vector 
   virtual void initializeData(); 
  
   //! For this class this method does nothing.
   virtual void ghostToGlobal(int /* mem */) {} 

   //! Nothing to do (its read only)
   virtual bool requiresDirichletAdjustment() const { return false; }

   //! Set the owned vector (Tpetra version)
   void setOwnedVector_Tpetra(const Teuchos::RCP<const VectorType>& ownedVector);

   //! Get the owned vector (Tpetra version)
   Teuchos::RCP<const VectorType> getOwnedVector_Tpetra() const;

   //! Get the ghosted vector (Tpetra version)
   Teuchos::RCP<VectorType> getGhostedVector_Tpetra() const;

   //! Set the owned vector (Thyra version)
   void setOwnedVector(const Teuchos::RCP<const Thyra::VectorBase<double> >& ownedVector);

   //! Get the owned vector (Thyra version)
   Teuchos::RCP<const Thyra::VectorBase<double> > getOwnedVector() const;

   //! Get the ghosted vector (Thyra version)
   Teuchos::RCP<Thyra::VectorBase<double> > getGhostedVector() const;

   //! Is this object initialized
   bool isInitialized() const { return isInitialized_; }

   //! Diagnostic function
   void print(std::ostream & os) const;

private:
   bool isInitialized_;

   Teuchos::RCP<const MapType> ghostedMap_;
   Teuchos::RCP<const MapType> ownedMap_;

   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ghostedSpace_;
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ownedSpace_;

   Teuchos::RCP<const ImportType> importer_;
   Teuchos::RCP<VectorType>       ghostedVector_;
   Teuchos::RCP<const VectorType> ownedVector_;

   typedef std::pair<std::vector<GlobalOrdinalT>,double> FilteredGlobalPair;
   typedef std::pair<std::vector<LocalOrdinalT>,double> FilteredLocalPair;
   std::vector<FilteredGlobalPair> globalFilteredPairs_;
   std::vector<FilteredLocalPair> filteredPairs_;
};

}

#endif
