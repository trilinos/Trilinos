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

#ifndef __Panzer_ResponseMESupportBase_hpp__
#define __Panzer_ResponseMESupportBase_hpp__

#include <string>

#include "Teuchos_RCP.hpp"

#include "PanzerDiscFE_config.hpp"
#ifdef PANZER_HAVE_EPETRA
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#endif

#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"

#include "Panzer_ResponseBase.hpp"

namespace panzer {

template <typename EvalT>
class ResponseMESupportBase : public ResponseBase {
public:
   ResponseMESupportBase(const std::string & responseName)
     : ResponseBase(responseName) {}

   virtual ~ResponseMESupportBase() {}

#ifdef PANZER_HAVE_EPETRA
   // This is the epetra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the <code>Epetra_Map</code> for this response, map is constructed lazily.
   virtual Teuchos::RCP<const Epetra_Map> getMap() const = 0;

   /** Set the vector (to be filled) for this response. This must be
     * constructed from the vector space returned by <code>getMap</code>.
     */
   virtual void setVector(const Teuchos::RCP<Epetra_Vector> & destVec) = 0;
#endif

   // This is the Thyra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the vector space for this response, vector space is constructed lazily.
   virtual Teuchos::RCP<const Thyra::VectorSpaceBase<double> > getVectorSpace() const = 0;

   /** Set the vector (to be filled) for this response. This must be
     * constructed from the vector space returned by <code>getVectorSpace</code>.
     */
   virtual void setVector(const Teuchos::RCP<Thyra::VectorBase<double> > & destVec) = 0;

private:
   // hide these methods
   ResponseMESupportBase();
   ResponseMESupportBase(const ResponseMESupportBase<EvalT> &);
};

template < >
class ResponseMESupportBase<panzer::Traits::Jacobian> : public ResponseBase {
public:
   ResponseMESupportBase(const std::string & responseName)
     : ResponseBase(responseName) {}

   virtual ~ResponseMESupportBase() {}

   //! Does this response support derivative evaluation?
   virtual bool supportsDerivative() const = 0;

#ifdef PANZER_HAVE_EPETRA
   // This is the epetra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the <code>Epetra_MultiVector</code> for this response, map is constructed lazily.
   virtual Teuchos::RCP<Epetra_MultiVector> buildEpetraDerivative() const = 0;

   /** Set the derivative (to be filled) for this response.
     */
   virtual void setDerivative(const Teuchos::RCP<Epetra_MultiVector> & derivative) = 0;
#endif

   // This is the Thyra view of the world
   ///////////////////////////////////////////////////////////

   //! Get an <code>Epetra_Operator</code> for this response, map is constructed lazily.
   virtual Teuchos::RCP<Thyra::MultiVectorBase<double> > buildDerivative() const = 0;

   /** Set the derivative (to be filled) for this response. This must be
     * constructed from the vector space returned by <code>getMap</code>.
     */
   virtual void setDerivative(const Teuchos::RCP<Thyra::MultiVectorBase<double> > & derivative) = 0;

private:
   // hide these methods
   ResponseMESupportBase();
   ResponseMESupportBase(const ResponseMESupportBase<panzer::Traits::Jacobian> &);
};

template < >
class ResponseMESupportBase<panzer::Traits::Tangent> : public ResponseBase {
public:
   ResponseMESupportBase(const std::string & responseName)
     : ResponseBase(responseName) {}

   virtual ~ResponseMESupportBase() {}

#ifdef PANZER_HAVE_EPETRA
   // This is the epetra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the <code>Epetra_Map</code> for this response, map is constructed lazily.
   virtual Teuchos::RCP<const Epetra_Map> getMap() const = 0;

   /** Set the vector (to be filled) for this response. This must be
     * constructed from the vector space returned by <code>getMap</code>.
     */
   virtual void setVector(const Teuchos::RCP<Epetra_MultiVector> & destVec) = 0;
#endif

   // This is the Thyra view of the world
   ///////////////////////////////////////////////////////////

   //! Get the vector space for this response, vector space is constructed lazily.
   virtual Teuchos::RCP<const Thyra::VectorSpaceBase<double> > getVectorSpace() const = 0;

   /** Set the vector (to be filled) for this response. This must be
     * constructed from the vector space returned by <code>getVectorSpace</code>.
     */
   virtual void setVector(const Teuchos::RCP<Thyra::MultiVectorBase<double> > & destVec) = 0;

private:
   // hide these methods
   ResponseMESupportBase();
   ResponseMESupportBase(const ResponseMESupportBase<panzer::Traits::Tangent> &);
};

#ifdef Panzer_BUILD_HESSIAN_SUPPORT
template < >
class ResponseMESupportBase<panzer::Traits::Hessian> : public ResponseBase {
public:
   ResponseMESupportBase(const std::string & responseName)
     : ResponseBase(responseName) {}

   virtual ~ResponseMESupportBase() {}

   //! Does this response support derivative evaluation?
   virtual bool supportsDerivative() const = 0;

   // This is the Thyra view of the world
   ///////////////////////////////////////////////////////////

   //! Get an <code>Epetra_Operator</code> for this response, map is constructed lazily.
   virtual Teuchos::RCP<Thyra::MultiVectorBase<double> > buildDerivative() const = 0;

   /** Set the derivative (to be filled) for this response. This must be
     * constructed from the vector space returned by <code>getMap</code>.
     */
   virtual void setDerivative(const Teuchos::RCP<Thyra::MultiVectorBase<double> > & derivative) = 0;

private:
   // hide these methods
   ResponseMESupportBase();
   ResponseMESupportBase(const ResponseMESupportBase<panzer::Traits::Hessian> &);
};

#endif

}

#endif
