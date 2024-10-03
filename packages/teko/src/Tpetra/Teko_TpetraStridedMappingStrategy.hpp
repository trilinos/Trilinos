// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_TpetraStridedMappingStrategy_hpp__
#define __Teko_TpetraStridedMappingStrategy_hpp__

// stl includes
#include <vector>

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"

// Tpetra includes
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_Export.hpp"
#include "Tpetra_CrsMatrix.hpp"

// Teko includes
#include "Teko_TpetraOperatorWrapper.hpp"
#include "Teko_ConfigDefs.hpp"

namespace Teko {
namespace TpetraHelpers {

class TpetraStridedMappingStrategy : public MappingStrategy {
 public:
  //! \name Constructors
  //@{

  /** Creates a strided mapping strategy. This class is useful
   * for breaking up nodally ordered matrices (i.e. the unknowns
   * in a FEM problem are ordered [u0,v0,p0,u1,v1,p1,...]). Current
   * implimentation only supports a fixed number of variables
   *
   * \param[in]      vars   Vector describing the blocking of variables
   * \param[in]      map    original Tpetra::Map<LO,GO,NT>  to be broken up
   * \param[in]      comm   Teuchos::Comm<int> object related to the map
   */
  TpetraStridedMappingStrategy(const std::vector<int>& vars,
                               const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >& map,
                               const Teuchos::Comm<int>& comm);
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
  virtual void copyTpetraIntoThyra(const Tpetra::MultiVector<ST, LO, GO, NT>& tpetra_X,
                                   const Teuchos::Ptr<Thyra::MultiVectorBase<ST> >& thyra_X) const;

  /** Virtual function defined in MappingStrategy.  This copies
   * an Tpetra::MultiVector<ST,LO,GO,NT> into a Thyra::MultiVectorBase with
   * blocking handled by the strides defined in the constructor.
   *
   * \param[in]     thyra_Y  source Thyra::MultiVectorBase
   * \param[in,out]     tpetra_Y destination Tpetra::MultiVector<ST,LO,GO,NT>
   */
  virtual void copyThyraIntoTpetra(const Teuchos::RCP<const Thyra::MultiVectorBase<ST> >& thyra_Y,
                                   Tpetra::MultiVector<ST, LO, GO, NT>& tpetra_Y) const;

  /** Returns the domain and range maps used by this class.
   * This faciliates building an Tpetra_Operator around this
   * class with its core functionality being a Thyra::LinearOpBase
   * operator
   *
   * \returns Range map corresponding to this class
   */
  virtual const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> > domainMap() const {
    return domainMap_;
  }

  /** Returns the domain and range maps used by this class.
   * This faciliates building an Tpetra_Operator around this
   * class with its core functionality being a Thyra::LinearOpBase
   * operator
   *
   * \returns Range map corresponding to this class
   */
  virtual const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> > rangeMap() const { return rangeMap_; }

  /** A function for my sanity
   *
   * \returns String with description of this class
   */
  virtual std::string toString() const { return std::string("TpetraStridedMappingStrategy"); }

  //@}

  //! \name Locally useful functions
  //@{

  /** \brief This is the core routine that builds the maps
   * and importers/exporters.
   *
   * This is the core routine that builds the maps
   * and importers/exporters neccessary for all the
   * transfers. Currently it simply calls out to the
   * interlaced tpetra functions. (Comment: this
   * routine should probably be private or protected
   * ... it is basically the meat of the constructor)
   *
   * \param[in]      vars      Vector describing the blocking of variables
   * \param[in]      baseMap   basic map to use in the transfers
   * \param[in]      comm      Teuchos::Comm<int> object
   */
  void buildBlockTransferData(const std::vector<int>& vars,
                              const Teuchos::RCP<const Tpetra::Map<LO, GO, NT> >& baseMap,
                              const Teuchos::Comm<int>& comm);

  /** \brief  Get the individual block maps underlying that
   * make up a strided vector/operator.
   *
   * Get the individual block maps underlying that
   * make up a strided vector/operator. These are
   * useful if you want to tear something ... i.e.
   * a matrix ... apart.
   *
   * \returns Return a vector of block maps
   *             created for this strided operator
   */
  const std::vector<std::pair<int, Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > >& getMaps() const {
    return blockMaps_;
  }

  /** Builds a blocked Thyra operator that uses the strided
   * mapping strategy to define sub blocks.
   *
   * \param[in] mat Tpetra::CrsMatrix<ST,LO,GO,NT>  with FillComplete called, this
   *             matrix is assumed to be square, with the same
   *             range and domain maps
   * \param[in] label String to be used to prefix the sub block's labels
   *
   * \returns Blocked Thyra linear operator with sub blocks
   *          defined by this mapping strategy
   */
  const Teuchos::RCP<Thyra::BlockedLinearOpBase<ST> > buildBlockedThyraOp(
      const Teuchos::RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> >& mat,
      const std::string& label = "<ANYM>") const;

  /** Rebuilds a block Thyra operator using the strided mapping
   * strategy to define sub blocks.
   *
   * \param[in] mat Tpetra::CrsMatrix<ST,LO,GO,NT>  with FillComplete called, this
   *             matrix is assumed to be square, with the same
   *             range and domain maps
   * \param[in] A Destination block linear op composed of blocks of
   *            Tpetra::CrsMatrix<ST,LO,GO,NT>  at all relevant locations
   */
  void rebuildBlockedThyraOp(const RCP<const Tpetra::CrsMatrix<ST, LO, GO, NT> >& mat,
                             const RCP<Thyra::BlockedLinearOpBase<ST> >& A) const;

  //@}

 protected:
  // member variables

  //! \name storage for sanity
  //@{
  Teuchos::RCP<const Tpetra::Map<LO, GO, NT> > domainMap_;
  Teuchos::RCP<const Tpetra::Map<LO, GO, NT> > rangeMap_;
  //@}

  //! \name block transfer data
  //@{
  std::vector<std::pair<int, Teuchos::RCP<Tpetra::Map<LO, GO, NT> > > > blockMaps_;
  std::vector<Teuchos::RCP<Tpetra::Export<LO, GO, NT> > > blockExport_;
  std::vector<Teuchos::RCP<Tpetra::Import<LO, GO, NT> > > blockImport_;
  //@}
};

}  // end namespace TpetraHelpers
}  // end namespace Teko

#endif
