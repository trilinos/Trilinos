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

#ifndef __Teko_BlockedMappingStrategy_hpp__
#define __Teko_BlockedMappingStrategy_hpp__

// stl includes
#include <vector>

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"

// Epetra includes
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Import.h"
#include "Epetra_Export.h"
#include "Epetra_CrsMatrix.h"

// Teko includes
#include "Teko_EpetraOperatorWrapper.hpp"
#include "Teko_BlockingEpetra.hpp"

namespace Teko {
namespace Epetra {

class BlockedMappingStrategy : public MappingStrategy {
 public:
  //! \name Constructors
  //@{

  /** Creates a strided mapping strategy. This class is useful
   * for breaking up nodally ordered matrices (i.e. the unknowns
   * in a FEM problem are ordered [u0,v0,p0,u1,v1,p1,...]). Current
   * implimentation only supports a fixed number of variables
   *
   * \param[in]      vars   Vector describing the blocking of variables
   * \param[in]      baseMap    original Epetra_Map to be broken up
   * \param[in]      comm   Epetra_Comm object related to the map
   */
  BlockedMappingStrategy(const std::vector<std::vector<int> >& vars,
                         const Teuchos::RCP<const Epetra_Map>& baseMap, const Epetra_Comm& comm);
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
  virtual std::string toString() const { return std::string("BlockedMappingStrategy"); }

  //@}

  //! \name Locally useful functions
  //@{

  /** \brief This is the core routine that builds the maps
   * and importers/exporters.
   *
   * This is the core routine that builds the maps
   * and importers/exporters neccessary for all the
   * transfers. Currently it simply calls out to the
   * interlaced epetra functions. (Comment: this
   * routine should probably be private or protected
   * ... it is basically the meat of the constructor)
   *
   * \param[in]      vars      Vector describing the blocking of variables
   * \param[in]      baseMap   basic map to use in the transfers
   * \param[in]      comm      Epetra_Comm object
   */
  void buildBlockTransferData(const std::vector<std::vector<int> >& vars,
                              const Teuchos::RCP<const Epetra_Map>& baseMap,
                              const Epetra_Comm& comm);

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
  const std::vector<Blocking::MapPair>& getMaps() const { return blockMaps_; }

  /** Builds a blocked Thyra operator that uses the strided
   * mapping strategy to define sub blocks.
   *
   * \param[in] mat Epetra_CrsMatrix with FillComplete called, this
   *             matrix is assumed to be square, with the same
   *             range and domain maps
   * \param[in] label String to be used to prefix the sub block's labels
   *
   * \returns Blocked Thyra linear operator with sub blocks
   *          defined by this mapping strategy
   */
  const Teuchos::RCP<Thyra::BlockedLinearOpBase<double> > buildBlockedThyraOp(
      const Teuchos::RCP<const Epetra_CrsMatrix>& mat, const std::string& label = "<ANYM>") const;

  /** Rebuilds a block Thyra operator using the strided mapping
   * strategy to define sub blocks.
   *
   * \param[in] mat Epetra_CrsMatrix with FillComplete called, this
   *             matrix is assumed to be square, with the same
   *             range and domain maps
   * \param[in] A Destination block linear op composed of blocks of
   *            Epetra_CrsMatrix at all relevant locations
   */
  void rebuildBlockedThyraOp(const RCP<const Epetra_CrsMatrix>& mat,
                             const RCP<Thyra::BlockedLinearOpBase<double> >& A) const;

  //@}

 protected:
  // member variables

  //! \name storage for sanity
  //@{
  Teuchos::RCP<const Epetra_Map> domainMap_;
  Teuchos::RCP<const Epetra_Map> rangeMap_;
  //@}

  //! \name block transfer data
  //@{
  std::vector<Blocking::MapPair> blockMaps_;
  std::vector<Teuchos::RCP<Epetra_Import> > blockImport_;
  std::vector<Teuchos::RCP<Epetra_Export> > blockExport_;
  //@}
};

}  // end namespace Epetra
}  // end namespace Teko

#endif
