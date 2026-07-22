// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_BasicMappingStrategy_hpp__
#define __Teko_BasicMappingStrategy_hpp__

// stl includes
#include <vector>

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_MultiVectorBase.hpp"

// Epetra includes
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"

// Teko includes
#include "Teko_EpetraOperatorWrapper.hpp"

namespace Teko {
namespace Epetra {

class BasicMappingStrategy : public MappingStrategy {
 public:
  //! \name Constructors
  //@{

  /** Creates a strided mapping strategy. This class is useful
   * for breaking up nodally ordered matrices (i.e. the unknowns
   * in a FEM problem are ordered [u0,v0,p0,u1,v1,p1,...]). Current
   * implimentation only supports a fixed number of variables
   *
   * \param[in]      rMap   Range map
   * \param[in]      dMap   Domain map
   * \param[in]      comm   Epetra_Comm object related to the map
   */
  BasicMappingStrategy(const Teuchos::RCP<const Epetra_Map>& rMap,
                       const Teuchos::RCP<const Epetra_Map>& dMap, const Epetra_Comm& comm);
  //@}

  //!\name Member functions inherited from Teko::Epetra::MappingStrategy
  //@{

  /** Virtual function defined in MappingStrategy.  This copies
   * an Epetra_MultiVector into a Thyra::MultiVectorBase with
   * blocking handled by the strides defined in the constructor.
   *
   * \param[in]     epetra_X  source Epetra_MultiVector
   * \param[in,out]     thyra_X   destination Thyra::MultiVectorBase
   * \param[in]     eow       Operator that defines the transition
   */
  virtual void copyEpetraIntoThyra(
      const Epetra_MultiVector& epetra_X,
      const Teuchos::Ptr<Thyra::MultiVectorBase<double> >& thyra_X) const;
  // const Teko::Epetra::EpetraOperatorWrapper & eow) const;

  /** Virtual function defined in MappingStrategy.  This copies
   * an Epetra_MultiVector into a Thyra::MultiVectorBase with
   * blocking handled by the strides defined in the constructor.
   *
   * \param[in]     thyra_Y  source Thyra::MultiVectorBase
   * \param[in,out]     epetra_Y destination Epetra_MultiVector
   * \param[in]     eow      Operator that defines the transition
   */
  virtual void copyThyraIntoEpetra(
      const Teuchos::RCP<const Thyra::MultiVectorBase<double> >& thyra_Y,
      Epetra_MultiVector& epetra_Y) const;
  // const Teko::Epetra::EpetraOperatorWrapper & eow) const;

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
  virtual std::string toString() const { return std::string("BasicMappingStrategy"); }

  //@}

 protected:
  // member variables

  //! \name storage for sanity
  //@{
  Teuchos::RCP<const Epetra_Map> domainMap_;
  Teuchos::RCP<const Epetra_Map> rangeMap_;
  //@}
};

}  // end namespace Epetra
}  // end namespace Teko

#endif
