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

#ifndef __Panzer_EpetraVector_Write_GlobalEvaluationData_hpp__
#define __Panzer_EpetraVector_Write_GlobalEvaluationData_hpp__

#include "Epetra_Export.h"
#include "Epetra_Vector.h"
#include "Epetra_Map.h"

#include "Teuchos_RCP.hpp"

#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"

#include "Panzer_WriteVector_GlobalEvaluationData.hpp"

namespace panzer {

/** This class provides a boundary exchange communication mechanism for vectors.
  * Note this provides a "write" interface for parameters 
  */
class EpetraVector_Write_GlobalEvaluationData : public WriteVector_GlobalEvaluationData {
public:
   //! Default constructor
   EpetraVector_Write_GlobalEvaluationData()
      : isInitialized_(false) { }

   EpetraVector_Write_GlobalEvaluationData(const EpetraVector_Write_GlobalEvaluationData & src)
      : isInitialized_(false) { initialize(src.exporter_, src.ghostedMap_, src.ownedMap_); }

   /** Initialize this object with some Epetra communication objects. This method
     * must be called before an object of this type can be used.
     *
     * \param[in] exporter Exporter for doing communication from the owned 
     *                     to the ghosted vector.
     * \param[in] ghostedMap Map describing the ghosted vector.
     * \param[in] ownedMap Map describing the ghosted vector.
     */
   EpetraVector_Write_GlobalEvaluationData(const Teuchos::RCP<const Epetra_Export>& exporter,
                                           const Teuchos::RCP<const Epetra_Map>&    ghostedMap,
                                           const Teuchos::RCP<const Epetra_Map>&    ownedMap)
      : isInitialized_(false) { initialize(exporter, ghostedMap, ownedMap); }

   /** Initialize this object with some Epetra communication objects. This method
     * must be called before an object of this type can be used.
     *
     * \param[in] exporter Exporter for doing communication from the owned 
     *                     to the ghosted vector.
     * \param[in] ghostedMap Map describing the ghosted vector.
     * \param[in] ownedMap Map describing the ghosted vector.
     */
   void initialize(const Teuchos::RCP<const Epetra_Export>& exporter,
                   const Teuchos::RCP<const Epetra_Map>&    ghostedMap,
                   const Teuchos::RCP<const Epetra_Map>&    ownedMap);

   /** For this class, this method does the halo exchange for the 
     * vector.
     */
   virtual void ghostToGlobal(int mem);

   //! Clear out the ghosted vector 
   virtual void initializeData(); 

   //! Nothing to do (its read only)
   virtual bool requiresDirichletAdjustment() const { return false; }

   //! Set the owned vector (Epetra version)
   void setOwnedVector_Epetra(const Teuchos::RCP<Epetra_Vector>&
      ownedVector);

   //! Get the ghosted vector (Epetra version)
   Teuchos::RCP<Epetra_Vector> getGhostedVector_Epetra() const;

   //! Set the owned vector (Thyra version)
   void setOwnedVector(const Teuchos::RCP<Thyra::VectorBase<double> >&
      ownedVector);

   //! Get the owned vector (Thyra version)
   Teuchos::RCP<Thyra::VectorBase<double> > getOwnedVector() const;

   //! Get the ghosted vector (Thyra version)
   Teuchos::RCP<Thyra::VectorBase<double> > getGhostedVector() const;

   //! Is this object initialized
   virtual bool isInitialized() const { return isInitialized_; }

   //! Diagnostic function
   void print(std::ostream & os) const;

private:
   bool isInitialized_;

   Teuchos::RCP<const Epetra_Map> ghostedMap_;
   Teuchos::RCP<const Epetra_Map> ownedMap_;

   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ghostedSpace_;
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ownedSpace_;

   Teuchos::RCP<const Epetra_Export> exporter_;
   Teuchos::RCP<Epetra_Vector>       ghostedVector_;

   Teuchos::RCP<Thyra::VectorBase<double> > ownedVector_;
};

}

#endif
