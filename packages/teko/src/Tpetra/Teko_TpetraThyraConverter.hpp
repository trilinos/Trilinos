// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_TpetraThyraConverter_hpp__
#define __Teko_TpetraThyraConverter_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Tpetra includes
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Map.hpp"

// Thyra includes
#include "Thyra_VectorSpaceBase.hpp"

#include "Teko_ConfigDefs.hpp"

namespace Teko {
namespace TpetraHelpers {

/** \brief Convert a Epetra_MultiVector with assumed block structure dictated by the
  *        vector space into a Thyra::MultiVectorBase object.
  *
  * Converts a Epetra_MultiVector to a Thyra::MultiVectorBase object. The structure of
  * the Epetra_MultiVector is assumed to be described by a map produced from the vector space
  * containing the destination vector. For instance if there is a vector \f$v = [ [v_0, v_1] ,
  v_2]^T\f$
  * where \f$v_0\f$, \f$v_1\f$ and \f$v_2\f$ are all subvectors then the Epetra_Map
  * corresponding to \f$v\f$ would have indicies

    \f$v_0^0, v_0^1,\ldots,v_0^{n_0},
       v_1^0, v_1^1,\ldots,v_1^{n_1}
       v_2^0, v_2^1,\ldots,v_2^{n_2} \f$.

  * That is the each of the subvectors are stacked on top of each other. The possibly recursive
  * block structure is then dictated by the vector space. The resulting Thyra::MultiVectorBase
  * object will be in the vectors space, with the contents of the Epetra_MultiVector copied to it.
  *
  * \param[in]     epetraX  Source Epetra_MultiVector object to be converted. See assumptions in
  Preconditions section.
  * \param[in,out] thyraX   Destination Thyra::MultiVectorBase. See assumptions in Preconditions
  section
  *
  * <b>Preconditions</b><ul>
  * <li> [mv.Map()==thyraVSToEpetraMap(thyraX->space())] This method assumes that the map used to
  create
  *      <code>epetraX</code> was created using the thyraVSToEpetraMap function using the vector
  space defining <code>thyraX</code>
  * <li> All subvectors are of type Thyra::SpmdMultiVectorBase or Thyra::ProductMultiVectorBase
  * </ul>
  *
  * <b>Postconditions</b><ul>
  * <li> [<code>thryaX==epetraX</code>] Contents of <code>epetraX</code> are copied into
  <code>thyraX</code>
  * </ul>
  *
  * \note Due to a Thyra issue with a somewhat incomplete inheritance hierarchy surrounding
  *       the <code>SpmdMultiVectorBase</code> and <code>SpmdVectorBase</code> interfaces
  <code>thyraX</code> must be of
  *       type <code>SpmdMultiVectorBase</code>, <code>ProductMultiVectorBase</code>, or
  <code>ProductVectorBase</code>.
  *       Notice that this does not include the <code>SpmdVectorBase</code> class. A fix of this
  might involve a more
  *       general implementation and use of <code>DetachedSpmdMultiVectorView</code>.
  */
void blockTpetraToThyra(const Tpetra::MultiVector<ST, LO, GO, NT>& tpetraX,
                        const Teuchos::Ptr<Thyra::MultiVectorBase<ST> >& thyraX);

/** \brief Convert a Thyra::MultiVectorBase object to a Epetra_MultiVector object with
  *        the map defined by the Epetra_Map.
  *
  * Converts a Thyra::MultiVectorBase object to a Epetra_MultiVector object. The Epetra_MultiVector
  * should have been created to be compatiable with the Thyra multi-vector, i.e. using the
  corresponding
  * map create by a call to thyraVSToEpetraMap. For example, for a Thyra::MultiVectorBase object
  * \f$v = [ [v_0, v_1] , v_2]^T\f$ where \f$v_0\f$, \f$v_1\f$ and \f$v_2\f$ are all subvectors,
  * the Epetra_MultiVector object will have global indicies

    \f$v_0^0, v_0^1,\ldots,v_0^{n_0},
       v_1^0, v_1^1,\ldots,v_1^{n_1}
       v_2^0, v_2^1,\ldots,v_2^{n_2} \f$.

  * That is the each of the subvectors are stacked on top of each other. The possibly recursive
  * block structure is then dictated by the Thyra::MultiVectorBase object. The resulting
  Epetra_MultiVector
  * object will be a copy of the Thyra::MultiVectorBase object.
  *
  * \param[in]     thyraX  Source Thyra::MultiVectorBase object to be converted. See assumptions in
  Preconditions section.
  * \param[in,out] epetraX Destination Epetra_MultiVector. See assumptions in Preconditions section
  *
  * <b>Preconditions</b><ul>
  * <li> [<code>epetraX.Map()==thyraVSToEpetraMap(thyraX->space())</code>] This method assumes that
  the map used to create
  *      <code>epetraX</code> was created using the thyraVSToEpetraMap function using the vector
  space defining <code>thyraX</code>
  * <li> All subvectors are of type Thyra::SpmdMultiVectorBase or Thyra::ProductMultiVectorBase
  * </ul>
  *
  * <b>Postconditions</b><ul>
  * <li> [<code>thryaX==epetraX</code>] Contents of <code>epetraX</code> are copied into
  <code>thyraX</code>
  * </ul>
  *
  * \note Due to a Thyra issue with a somewhat incomplete inheritance hierarchy surrounding
  *       the <code>SpmdMultiVectorBase</code> and <code>SpmdVectorBase</code> interfaces
  <code>thyraX</code> must be of
  *       type <code>SpmdMultiVectorBase</code>, <code>ProductMultiVectorBase</code>, or
  <code>ProductVectorBase</code>.
  *       Notice that this does not include the <code>SpmdVectorBase</code> class. A fix of this
  might involve a more
  *       general implementation and use of <code>ConstDetachedSpmdMultiVectorView</code>.
  */
void blockThyraToTpetra(const Teuchos::RCP<const Thyra::MultiVectorBase<ST> >& thyraX,
                        Tpetra::MultiVector<ST, LO, GO, NT>& tpetraX);

/** \brief From a Thyra vector space create a compatable Epetra_Map
  *
  * Build a distributed Epetra_Map from a Thyra::VectorSpaceBase object. This vector space
  * should only be composed of ProductVectorSpaceBase and SpmdVectorSpaceBase objects. Elements
  * stored locally on this process will map to global indicies stored locally. So the parallel
  * distribution should be the same for both the vector space and the map.
  *
  * As an example consider the vector space that describes vectors of the form
  * \f$v = [ [v_0, v_1] , v_2]^T\f$ where \f$v_0\f$, \f$v_1\f$ and \f$v_2\f$ are all subvectors.
  * The Epetra_Map created from this vector space will have global indicies

    \f$v_0^0, v_0^1,\ldots,v_0^{n_0},
       v_1^0, v_1^1,\ldots,v_1^{n_1}
       v_2^0, v_2^1,\ldots,v_2^{n_2} \f$.

  * That is the each of the subvectors are stacked on top of each other.
  *
  * \param[in] vs   A vector space object from which the Epetra_Map is to be created
  * \param[in] comm The Epetra_Comm object used to create the Epetra_Map
  *
  * \return An Epetra_Map is returned whose global indicies are distributed in a
  *         way that mirrors the structure of the vector space <code>vs</code>.
  *
  * <b>Preconditions</b><ul>
  * <li> Assumes that the vector space is composed of Thyra::ProductMultiVectorBase and
  *      Thyra::SpmdMultiVectorBase objects.
  * </ul>
  */
const Teuchos::RCP<Tpetra::Map<LO, GO, NT> > thyraVSToTpetraMap(
    const Thyra::VectorSpaceBase<ST>& vs,
    const Teuchos::RCP<const Teuchos::Comm<Thyra::Ordinal> >& comm);

}  // namespace TpetraHelpers
}  // end namespace Teko

#endif
