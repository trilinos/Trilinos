// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_ReorderedMappingStrategy_hpp__
#define __Teko_ReorderedMappingStrategy_hpp__

// stl includes
#include <vector>

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"

// Epetra includes
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"

// Teko includes
#include "Teko_BlockedReordering.hpp"
#include "Teko_EpetraOperatorWrapper.hpp"

namespace Teko {
namespace Epetra {

class ReorderedMappingStrategy : public MappingStrategy {
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
  ReorderedMappingStrategy(const BlockReorderManager& brm,
                           const Teuchos::RCP<const MappingStrategy>& map);
  //@}

  //!\name Member functions inherited from Teko::Epetra::MappingStrategy
  //@{

  /** Virtual function defined in MappingStrategy.  This copies
   * an Epetra_MultiVector into a Thyra::MultiVectorBase with
   * blocking handled by the strides defined in the constructor.
   *
   * \param[in]     epetra_X  source Epetra_MultiVector
   * \param[in,out]     thyra_X   destination Thyra::MultiVectorBase
   */
  virtual void copyEpetraIntoThyra(
      const Epetra_MultiVector& epetra_X,
      const Teuchos::Ptr<Thyra::MultiVectorBase<double> >& thyra_X) const;

  /** Virtual function defined in MappingStrategy.  This copies
   * an Epetra_MultiVector into a Thyra::MultiVectorBase with
   * blocking handled by the strides defined in the constructor.
   *
   * \param[in]     thyra_Y  source Thyra::MultiVectorBase
   * \param[in,out]     epetra_Y destination Epetra_MultiVector
   */
  virtual void copyThyraIntoEpetra(
      const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& thyra_Y,
      Epetra_MultiVector& epetra_Y) const;

  /** Returns the domain and range maps used by this class.
   * This faciliates building an Epetra_Operator around this
   * class with its core functionality being a Thyra::LinearOpBase
   * operator
   *
   * \returns Range map corresponding to this class
   */
  virtual const Teuchos::RCP<const Epetra_Map> domainMap() const { return domainMap_; }

  /** Returns the domain and range maps used by this class.
   * This faciliates building an Epetra_Operator around this
   * class with its core functionality being a Thyra::LinearOpBase
   * operator
   *
   * \returns Range map corresponding to this class
   */
  virtual const Teuchos::RCP<const Epetra_Map> rangeMap() const { return rangeMap_; }

  /** A function for my sanity
   *
   * \returns String with description of this class
   */
  virtual std::string toString() const { return std::string("ReorderedMappingStrategy"); }

  //@}

  //! Get the underlying maping strategy used by this object
  virtual RCP<const MappingStrategy> GetPrimaryMapStrategy() const { return mapStrategy_; }

 protected:
  // member variables

  //! \name storage for sanity
  //@{
  Teuchos::RCP<const Epetra_Map> domainMap_;
  Teuchos::RCP<const Epetra_Map> rangeMap_;
  //@}

  //! \name Reordering information
  //@{
  const BlockReorderManager reorderManager_;
  RCP<const MappingStrategy> mapStrategy_;
  //@}
};

}  // end namespace Epetra
}  // end namespace Teko

#endif
