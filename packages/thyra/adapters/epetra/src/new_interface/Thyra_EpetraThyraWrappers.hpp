// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef THYRA_EPETRA_THYRA_WRAPPERS_HPP
#define THYRA_EPETRA_THYRA_WRAPPERS_HPP

// Not directly needed in this file, but this way we made the
// macro HAVE_THYRA_EPETRA_REFACTOR available to files that include
// this header. This way, they do not need to include the config.h
// header manually. That's nice, because in the future we may deprecate
// and then remove the old interface, making the config.h file pointless.
// If that happens, we may remove it, and at that point all files including
// it would have to be updated. This was, only the adapters headers need to
// be updated.
#include "ThyraEpetraAdapters_config.h"

#include "Thyra_LinearOpBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "Teuchos_Comm.hpp"

class Epetra_Comm;
class Epetra_Map;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_Operator;

namespace Teuchos { template<class Ordinal> class Comm; }


namespace Thyra {

/** \brief Given an Epetra <tt>Teuchos::Comm<int></tt> object, return an
 * equivalent <tt>Teuchos::Comm<Ordinal></tt> object.
 *
 * Will throw if conversion is not successful.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<const Teuchos::Comm<Ordinal> >
convertEpetraToThyraComm( const Epetra_Comm& epetraComm );


/** \brief Create a Thyra::VectorSpaceBase object given a Epetra_Map.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<const VectorSpaceBase<double> >
createVectorSpace(const RCP<const Epetra_Map>& epetraMap);


/** \brief .
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<VectorBase<double> >
createVector(
  const RCP<Epetra_Vector>& epetraVector,
  const RCP<const VectorSpaceBase<double> > space = Teuchos::null
  );


/** \brief .
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<const VectorBase<double> >
createConstVector(
  const RCP<const Epetra_Vector>& epetraVector,
  const RCP<const VectorSpaceBase<double> > space = Teuchos::null
  );


/** \brief .
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<MultiVectorBase<double> >
createMultiVector(
  const RCP<Epetra_MultiVector>& epetraMultiVector,
  const RCP<const VectorSpaceBase<double> > rangeSpace = Teuchos::null,
  const RCP<const VectorSpaceBase<double> > domainSpace = Teuchos::null
  );


/** \brief .
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<const MultiVectorBase<double> >
createConstMultiVector(
  const RCP<const Epetra_MultiVector>& epetraMultiVector,
  const RCP<const VectorSpaceBase<double> > rangeSpace = Teuchos::null,
  const RCP<const VectorSpaceBase<double> > domainSpace = Teuchos::null
  );


/** \brief .
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<LinearOpBase<double> >
createLinearOp(
  const RCP<Epetra_Operator>& epetraOperator,
  const RCP<const VectorSpaceBase<double> > rangeSpace = Teuchos::null,
  const RCP<const VectorSpaceBase<double> > domainSpace = Teuchos::null
  );


/** \brief .
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
RCP<const LinearOpBase<double> >
createConstLinearOp(
  const RCP<const Epetra_Operator>& epetraOperator,
  const RCP<const VectorSpaceBase<double> > rangeSpace = Teuchos::null,
  const RCP<const VectorSpaceBase<double> > domainSpace = Teuchos::null
  );


/** \brief Traits class that enables the extraction of Epetra operator/vector
 * objects wrapped in Thyra operator/vector objects.
 *
 * Example usage:

 \code

  typedef Thyra::EpetraObjectExtraction TOE;
  typedef Epetra_MultiVector EpetraMultiVector_t;

  RCP<EpetraMultiVector_t> epetraMv = TOE::getEpetraMultiVector(thyraMv);
  RCP<EpetraVector_t> epetraV = TOE::getEpetraVector(thyraV);

 \endcode

 *
 * \todo Finish documentation
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
class EpetraOperatorVectorExtraction {
public:

  /** \brief Get a non-const Epetra_Vector from a non-const
   * Thyra::VectorBase object.
   */
  static RCP<Epetra_Vector>
  getEpetraVector(const RCP<VectorBase<double> > &v);

  /** \brief Get a const Epetra_Vector from a const
   * Thyra::VectorBase object.
   */
  static RCP<const Epetra_Vector>
  getConstEpetraVector(const RCP<const VectorBase<double> > &v);

  /** \brief Get a non-const Epetra_MultiVector from a non-const
   * Thyra::MultiVectorBase object.
   */
  static RCP<Epetra_MultiVector>
  getEpetraMultiVector(const RCP<MultiVectorBase<double> > &mv);

  /** \brief Get a const Epetra_MultiVector from a const
   * Thyra::MultiVectorBase object.
   */
  static RCP<const Epetra_MultiVector>
  getConstEpetraMultiVector(const RCP<const MultiVectorBase<double> > &mv);

  /** \brief Get a non-const Epetra_Operator from a non-const
   * Thyra::LinearOpBase object.
   */
  static RCP<Epetra_Operator>
  getEpetraOperator(const RCP<LinearOpBase<double> > &op);

  /** \brief Get a const Epetra_Operator from a const
   * Thyra::LinearOpBase object.
   */
  static RCP<const Epetra_Operator>
  getConstEpetraOperator(const RCP<const LinearOpBase<double> > &op);

  ///////////////////////////////////////////////////
  //             getOrCreate methods               //
  ///////////////////////////////////////////////////
  /** \brief Get an Epetra_* object from a generic
   * Thyra one. If the Thyra object is of concrete
   * type EpetraXYZ, simply get the Epetra object
   * out from it. If not, try to build an Epetra
   * object from the input.
   */

  // Given a vector space, retrieves the Epetra_Map (if the vs was of type
  // EpetraVectorSpace) or creates it from scratch, in case it is another
  // type of vector space.
  static RCP<const Epetra_Map>
  getOrCreateEpetraMap (const RCP<const VectorSpaceBase<double>> &vs);

  static RCP<Epetra_Vector>
  getOrCreateEpetraVector(const RCP<VectorBase<double> > &v);

  static RCP<const Epetra_Vector>
  getOrCreateConstEpetraVector(const RCP<const VectorBase<double> > &v);

  static RCP<Epetra_MultiVector>
  getOrCreateEpetraMultiVector(const RCP<MultiVectorBase<double> > &mv);

  static RCP<const Epetra_MultiVector>
  getOrCreateConstEpetraMultiVector(const RCP<const MultiVectorBase<double> > &mv);

  static RCP<Epetra_Operator>
  getOrCreateEpetraOperator(const RCP<LinearOpBase<double> > &op);

  static RCP<const Epetra_Operator>
  getOrCreateConstEpetraOperator(const RCP<const LinearOpBase<double> > &op);

  static RCP<const Epetra_Comm>
  createEpetraComm(const RCP<const Teuchos::Comm<Ordinal>> comm);

};

} // namespace Thyra

#endif // THYRA_EPETRA_THYRA_WRAPPERS_HPP
