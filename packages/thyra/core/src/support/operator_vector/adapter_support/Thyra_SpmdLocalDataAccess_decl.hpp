// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_SPMD_LOCAL_DATA_ACCESS_DECL_HPP
#define THYRA_SPMD_LOCAL_DATA_ACCESS_DECL_HPP


#include "Thyra_OperatorVectorTypes.hpp"


namespace Thyra {


/** \defgroup Thyra_Op_Vec_adapters_spmd_local_data_access_grp Official utilities for accessing local data in SPMD vectors and multi-vectors.
 *
 * \ingroup Thyra_Op_Vec_spmd_adapters_grp
 *
 * \brief These non-member helper functions provide the standard way by which
 * clients can get access to local data in VectorBase and MultiVector base
 * objects for the standard SPMD objects.  They are the primary
 * interoperability mechanism for clients to get at local SPMD VectorBase and
 * MultiVectorBase data.
 *
 * All of these views follow the constraints of the Generalized View Design
 * Pattern defined
 * <a href="../../../teuchos/doc/html/TeuchosMemoryManagementSAND.pdf">here</a>.
 * That means a few things:
 * <ul>
 *
 * <li>Unlimited const views are allowed but only a single non-const view is
 * allowed.</li>
 *
 * <li>Once a single non-const view is created, no other const or non-const
 * views are allowed.</li>
 *
 * <li>If any view is active, a non-const view can not be created.</li>
 *
 * <li>Changes to an object through a non-const view are not guaranteed to be
 *     reflected in the parent object until the view is released.</li>
 *
 * </ul>
 *
 * ToDo: Finish Documentation!
 */


/** \brief Return a contiguous non-const semi-persisting view of the local
 * process data of a VectorBase object.
 *
 * \ingroup Thyra_Op_Vec_adapters_spmd_local_data_access_grp
 */
template<class Scalar>
RTOpPack::SubVectorView<Scalar>
getNonconstLocalSubVectorView(
  const RCP<VectorBase<Scalar> > &vec);


/** \brief Return a contiguous const semi-persisting view of the local process
 * data of a VectorBase object.
 *
 * \ingroup Thyra_Op_Vec_adapters_spmd_local_data_access_grp
 */
template<class Scalar>
RTOpPack::ConstSubVectorView<Scalar>
getLocalSubVectorView(
  const RCP<const VectorBase<Scalar> > &vec);


/** \brief Return a contiguous non-const semi-persisting view of the local
 * process data of a MultiVectorBase object.
 *
 * \ingroup Thyra_Op_Vec_adapters_spmd_local_data_access_grp
 */
template<class Scalar>
RTOpPack::SubMultiVectorView<Scalar>
getNonconstLocalSubMultiVectorView(
  const RCP<MultiVectorBase<Scalar> > &multivec);


/** \brief Return a contiguous const semi-persisting view of the local process
 * data of a MultiVectorBase object.
 *
 * \ingroup Thyra_Op_Vec_adapters_spmd_local_data_access_grp
 */
template<class Scalar>
RTOpPack::ConstSubMultiVectorView<Scalar>
getLocalSubMultiVectorView(
  const RCP<const MultiVectorBase<Scalar> > &multivec);


} // end namespace Thyra


#endif // THYRA_SPMD_LOCAL_DATA_ACCESS_DECL_HPP
