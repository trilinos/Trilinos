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

#ifndef __Panzer_ResponseMESupport_Default_impl_hpp__
#define __Panzer_ResponseMESupport_Default_impl_hpp__

namespace panzer {

#ifdef PANZER_HAVE_EPETRA
template <typename EvalT>
Epetra_Vector & ResponseMESupport_Default<EvalT>::
getEpetraVector() const
{
   TEUCHOS_ASSERT(useEpetra());

   return *eVector_;
}
#endif

template <typename EvalT>
Thyra::ArrayRCP<double> ResponseMESupport_Default<EvalT>::
getThyraVector() const
{
   TEUCHOS_ASSERT(useThyra());

   Teuchos::ArrayRCP<double> data;
   Teuchos::rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(tVector_,true)->getNonconstLocalData(Teuchos::outArg(data));

   return data;
}

#ifdef PANZER_HAVE_EPETRA
template <typename EvalT>
Teuchos::RCP<const Epetra_Map> ResponseMESupport_Default<EvalT>::
getMap() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(useThyra_,std::logic_error,
                             "Reponse field \"" << this->getName() << "\" has previously been initialized as a "
                             "Thyra object, now trying to initalize as a Epetra! Error!");

  // lazily construct the map only as needed
  if(map_==Teuchos::null) {
    if(this->vectorIsDistributed())
      map_ = Teuchos::rcp(new Epetra_Map(-1,(int) this->localSizeRequired(),0,eComm_));
    else
      map_ = Teuchos::rcp(new Epetra_LocalMap((int) this->localSizeRequired(),0,eComm_));
  }

  return map_;
}

template <typename EvalT>
void ResponseMESupport_Default<EvalT>::
setVector(const Teuchos::RCP<Epetra_Vector> & destVec)
{
  TEUCHOS_TEST_FOR_EXCEPTION(useThyra_,std::logic_error,
                             "Reponse field \"" << this->getName() << "\" has previously been initialized as a "
                             "Thyra object, now trying to initalize as a Epetra! Error!");

  eVector_ = destVec;

  useEpetra_ = true;
}
#endif

template <typename EvalT>
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > ResponseMESupport_Default<EvalT>::
getVectorSpace() const
{
  TEUCHOS_TEST_FOR_EXCEPTION(useEpetra_,std::logic_error,
                             "Reponse field \"" << this->getName() << "\" has previously been initialized as an "
                             "Epetra object, now trying to initalize as a Thyra object! Error!");

  // lazily build the space and return it
  if(vSpace_==Teuchos::null) {
    if(this->vectorIsDistributed())
      vSpace_ = Thyra::defaultSpmdVectorSpace<double>(tComm_,this->localSizeRequired(),-1);
    else
      vSpace_ = Thyra::locallyReplicatedDefaultSpmdVectorSpace<double>(tComm_,this->localSizeRequired());
  }

  return vSpace_;
}

template <typename EvalT>
void ResponseMESupport_Default<EvalT>::
setVector(const Teuchos::RCP<Thyra::VectorBase<double> > & destVec)
{
  TEUCHOS_TEST_FOR_EXCEPTION(useEpetra_,std::logic_error,
                             "Reponse field \"" << this->getName() << "\" has previously been initialized as an "
                             "Epetra object, now trying to initalize as a Thyra object! Error!");

  tVector_ = destVec;

  useThyra_ = true;
}

}

#endif
