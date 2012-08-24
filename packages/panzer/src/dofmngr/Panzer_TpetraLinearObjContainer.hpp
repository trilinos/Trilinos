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

#ifndef __Panzer_TpetraLinearObjContainer_hpp__
#define __Panzer_TpetraLinearObjContainer_hpp__

#include "Panzer_config.hpp"

#include <map>

// Tpetra includes
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "Thyra_TpetraThyraWrappers.hpp"

#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_ThyraObjContainer.hpp"

#include "Teuchos_RCP.hpp"

namespace panzer {

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT=Kokkos::DefaultNode::DefaultNodeType>
class TpetraLinearObjContainer : public LinearObjContainer
                               , public ThyraObjContainer<ScalarT> {
   TpetraLinearObjContainer();

public:
   typedef LinearObjContainer::Members Members;

   typedef Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> VectorType;
   typedef Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> CrsMatrixType;
   typedef Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> CrsGraphType;
   typedef Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> MapType;
   typedef Tpetra::Import<LocalOrdinalT,GlobalOrdinalT,NodeT> ImportType;
   typedef Tpetra::Export<LocalOrdinalT,GlobalOrdinalT,NodeT> ExportType;

   TpetraLinearObjContainer(const Teuchos::RCP<const Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > & domain,
                            const Teuchos::RCP<const Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> > & range)
   {
      domainSpace = Thyra::createVectorSpace<ScalarT>(domain);
      rangeSpace = Thyra::createVectorSpace<ScalarT>(range);
   }

   virtual void initialize() 
   {
      if(get_x()!=Teuchos::null) get_x()->putScalar(0.0);
      if(get_dxdt()!=Teuchos::null) get_dxdt()->putScalar(0.0);
      if(get_f()!=Teuchos::null) get_f()->putScalar(0.0);
      if(get_A()!=Teuchos::null) get_A()->setAllToScalar(0.0);
   }

   //! Wipe out stored data.
   void clear()
   {
      set_x(Teuchos::null);
      set_dxdt(Teuchos::null);
      set_f(Teuchos::null);
      set_A(Teuchos::null);
   }

   inline void set_x(const Teuchos::RCP<VectorType> & in) { x = in; } 
   inline const Teuchos::RCP<VectorType> get_x() const { return x; }

   inline void set_dxdt(const Teuchos::RCP<VectorType> & in) { dxdt = in; } 
   inline const Teuchos::RCP<VectorType> get_dxdt() const { return dxdt; }

   inline void set_f(const Teuchos::RCP<VectorType> & in) { f = in; } 
   inline const Teuchos::RCP<VectorType> get_f() const { return f; }

   inline void set_A(const Teuchos::RCP<CrsMatrixType> & in) { A = in; } 
   inline const Teuchos::RCP<CrsMatrixType> get_A() const { return A; }

   void initializeMatrix(ScalarT value)
   { A->setAllToScalar(value); }

   virtual void set_x_th(const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & in) 
   { x = (in==Teuchos::null) ? Teuchos::null : TOE::getTpetraVector(in); } 
   virtual Teuchos::RCP<Thyra::VectorBase<ScalarT> > get_x_th() const 
   { return (x==Teuchos::null) ? Teuchos::null : Thyra::createVector(x,domainSpace); }

   virtual void set_dxdt_th(const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & in)
   { dxdt = (in==Teuchos::null) ? Teuchos::null : TOE::getTpetraVector(in); } 
   virtual Teuchos::RCP<Thyra::VectorBase<ScalarT> > get_dxdt_th() const 
   { return (dxdt==Teuchos::null) ? Teuchos::null : Thyra::createVector(dxdt,domainSpace); }

   virtual void set_f_th(const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & in)
   { f = (in==Teuchos::null) ? Teuchos::null : TOE::getTpetraVector(in); } 
   virtual Teuchos::RCP<Thyra::VectorBase<ScalarT> > get_f_th() const 
   { return (f==Teuchos::null) ? Teuchos::null : Thyra::createVector(f,rangeSpace); }

   virtual void set_A_th(const Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > & in) 
   { A = (in==Teuchos::null) ? Teuchos::null : Teuchos::rcp_dynamic_cast<CrsMatrixType>(TOE::getTpetraOperator(in),true); }
   virtual Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > get_A_th() const
   { return (A==Teuchos::null) ? Teuchos::null : Thyra::createLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(A,rangeSpace,domainSpace); }
    
private:
   typedef Thyra::TpetraOperatorVectorExtraction<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> TOE;

   Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > domainSpace;
   Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > rangeSpace;

   Teuchos::RCP<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > x, dxdt, f;
   Teuchos::RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > A;
};

}

#endif
