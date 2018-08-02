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


#ifndef THYRA_EPETRA_VECTOR_SPACE_HPP
#define THYRA_EPETRA_VECTOR_SPACE_HPP

// Not directly needed in this file, but this way we made the
// macro HAVE_THYRA_EPETRA_REFACTOR available to files that include
// this header. This way, they do not need to include the config.h
// header manually. That's nice, because in the future we may deprecate
// and then remove the old interface, making the config.h file pointless.
// If that happens, we may remove it, and at that point all files including
// it would have to be updated. This was, only the adapters headers need to
// be updated.
#include "ThyraEpetraAdapters_config.h"

#include "Thyra_SpmdVectorSpaceDefaultBase.hpp"

class Epetra_Map;
class Epetra_Vector;
class Epetra_MultiVector;

namespace Thyra {


/** \brief Concrete implementation of an SPMD vector space for Epetra.
 *
 * \ingroup Epetra_Thyra_Op_Vec_adapters_grp
 */
class EpetraVectorSpace : public SpmdVectorSpaceDefaultBase<double>
{
public:

  /** @name Constructors and initializers */
  //@{

  /** \brief Create with weak ownership to self. */
  static RCP<EpetraVectorSpace> create();

  /** \brief Initialize a serial space. */
  void initialize(const RCP<const Epetra_Map> &epetraMap);

  //@}

  /** @name Public overridden from VectorSpaceBase */
  //@{
  /** \brief Returns true if all the elements in <tt>rng</tt> are in this
   * process.
   */
  bool hasInCoreView(
    const Range1D& rng, const EViewType viewType, const EStrideType strideType
    ) const;
  /** \brief . */
  RCP< const VectorSpaceBase<double> > clone() const;
  //@}

  // This method is helpful so that Epetra(Multi)Vector constructor/init methods
  // that accept an Epetra_(Multi)Vector and a VectorSpaceBase can check that
  //   1) the vector space is an EpetraVectorSpace
  //   2) the EpetraVectorSpace stores the same map of the Epetra_(Multi)Vector
  RCP<const Epetra_Map> getEpetraMap() const;

protected:

  /** @name Protected overridden from VectorSpaceBase */
  //@{

  /** \brief . */
  RCP<VectorBase<double>> createMember() const;
  /** \brief . */
  RCP<MultiVectorBase<double>> createMembers(int numMembers) const;

  //@}

public:

  /** @name Public overridden from SpmdVectorSpaceBase */
  //@{

  /** \brief . */
  RCP<const Teuchos::Comm<Ordinal> > getComm() const;
  /** \brief . */
  Ordinal localSubDim() const;

  //@}

private:

  // //////////////////////////////////////
  // Private data members

  RCP<const Epetra_Map> epetraMap_;
  // The only reason Thyra needs this comm_ object is because Thyra
  // uses Ordinal as the Comm template parameter, while Epetra uses
  // int.  Ordinal is some 64-bit type, which doesn't make any sense,
  // given that MPI implementations currently only allow 32-bit
  // process ranks.  This is why Thyra does not just use the Map's
  // stored communicator.
  RCP<const Teuchos::Comm<Ordinal> > comm_;
  RCP<EpetraVectorSpace> weakSelfPtr_;

  // /////////////////////////////////////
  // Private member functions

  EpetraVectorSpace();

}; // end class EpetraVectorSpace


/** \brief Nonmember consturctor that creats a serial vector space.
 *
 * \relates EpetraVectorSpace
 */
inline
RCP<EpetraVectorSpace>
epetraVectorSpace(const RCP<const Epetra_Map> &epetraMap)
{
  RCP<EpetraVectorSpace> vs = EpetraVectorSpace::create();
  vs->initialize(epetraMap);
  return vs;
}

} // end namespace Thyra

#endif // THYRA_EPETRA_VECTOR_SPACE_HPP
