/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER

*/

#ifndef __Teko_TpetraReorderedMappingStrategy_hpp__
#define __Teko_TpetraReorderedMappingStrategy_hpp__

// stl includes
#include <vector>

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"

// Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_CrsMatrix.hpp"

// Teko includes
#include "Teko_BlockedReordering.hpp"
#include "Teko_TpetraOperatorWrapper.hpp"
#include "Teko_ConfigDefs.hpp"

namespace Teko {
namespace TpetraHelpers {

class TpetraReorderedMappingStrategy : public MappingStrategy {
public:
   //! \name Constructors
   //@{

   /** Creates a reordered mapping strategy. This class is useful
     * for using a reordering manager to build a new mapping strategy.
     * This is new strategy is based on some prior mapping.
     *
     * \param[in]      brm   How the vectors are to be reordered
     * \param[in]      map   Original mapping strategy to base the reordering on
     */
   TpetraReorderedMappingStrategy(const BlockReorderManager & brm,const Teuchos::RCP<const MappingStrategy> & map);
   //@}

   //!\name Member functions inherited from Teko::Tpetra::MappingStrategy
   //@{

   /** Virtual function defined in MappingStrategy.  This copies
     * an Tpetra::MultiVector<ST,LO,GO,NT> into a Thyra::MultiVectorBase with
     * blocking handled by the strides defined in the constructor.
     *
     * \param[in]     tpetra_X  source Tpetra::MultiVector<ST,LO,GO,NT>
     * \param[in,out]     thyra_X   destination Thyra::MultiVectorBase
     */
   virtual void copyTpetraIntoThyra(const Tpetra::MultiVector<ST,LO,GO,NT>& tpetra_X, 
                                    const Teuchos::Ptr<Thyra::MultiVectorBase<ST> > & thyra_X) const;

   /** Virtual function defined in MappingStrategy.  This copies
     * an Tpetra::MultiVector<ST,LO,GO,NT> into a Thyra::MultiVectorBase with
     * blocking handled by the strides defined in the constructor.
     *
     * \param[in]     thyra_Y  source Thyra::MultiVectorBase
     * \param[in,out]     tpetra_Y destination Tpetra::MultiVector<ST,LO,GO,NT>
     */
   virtual void copyThyraIntoTpetra(const Teuchos::RCP<const Thyra::MultiVectorBase<ST> > & thyra_Y, 
                                    Tpetra::MultiVector<ST,LO,GO,NT>& tpetra_Y) const;

   /** Returns the domain and range maps used by this class.
     * This faciliates building an Tpetra_Operator around this
     * class with its core functionality being a Thyra::LinearOpBase
     * operator
     *
     * \returns Range map corresponding to this class
     */
   virtual const Teuchos::RCP<const Tpetra::Map<LO,GO,NT> > domainMap() const
   { return domainMap_; }

   /** Returns the domain and range maps used by this class.
     * This faciliates building an Tpetra_Operator around this
     * class with its core functionality being a Thyra::LinearOpBase
     * operator
     *
     * \returns Range map corresponding to this class
     */
   virtual const Teuchos::RCP<const Tpetra::Map<LO,GO,NT> > rangeMap() const
   { return rangeMap_; }

   /** A function for my sanity
     *
     * \returns String with description of this class
     */
   virtual std::string toString() const
   { return std::string("TpetraReorderedMappingStrategy"); }
 
   //@}

   //! Get the underlying maping strategy used by this object
   virtual RCP<const MappingStrategy> GetPrimaryMapStrategy() const
   { return mapStrategy_; }

protected:
   // member variables

   //! \name storage for sanity
   //@{
   Teuchos::RCP<const Tpetra::Map<LO,GO,NT> > domainMap_; 
   Teuchos::RCP<const Tpetra::Map<LO,GO,NT> > rangeMap_;
   //@}

   //! \name Reordering information
   //@{
   const BlockReorderManager reorderManager_;
   RCP<const MappingStrategy> mapStrategy_;
   //@}
};

} // end namespace TpetraHelpers
} // end namespace Teko

#endif
