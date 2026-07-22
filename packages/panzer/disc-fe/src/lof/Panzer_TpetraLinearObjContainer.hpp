// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_TpetraLinearObjContainer_hpp__
#define __Panzer_TpetraLinearObjContainer_hpp__

#include "PanzerDiscFE_config.hpp"

#include <map>

// Tpetra includes
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"

#include "Thyra_TpetraThyraWrappers.hpp"

#include "Panzer_LinearObjFactory.hpp"
#include "Panzer_ThyraObjContainer.hpp"
#include "Panzer_NodeType.hpp"

#include "Teuchos_RCP.hpp"

namespace panzer {

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT=panzer::TpetraNodeType>
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
      if(get_A()!=Teuchos::null) {
        Teuchos::RCP<CrsMatrixType> mat = get_A(); 
        mat->setAllToScalar(0.0);
      }
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
   {  
     A->setAllToScalar(value); 
   }

   virtual void set_x_th(const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & in) 
   { 
     if(in==Teuchos::null) { x = Teuchos::null; return; }

     Teuchos::RCP<const Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > x_const 
         = TOE::getConstTpetraVector(in);
     x = Teuchos::rcp_const_cast<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(x_const); 
   } 
   virtual Teuchos::RCP<Thyra::VectorBase<ScalarT> > get_x_th() const 
   { return (x==Teuchos::null) ? Teuchos::null : Thyra::createVector(x,domainSpace); }

   virtual void set_dxdt_th(const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & in)
   { 
     if(in==Teuchos::null) { dxdt = Teuchos::null; return; }

     Teuchos::RCP<const Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > dxdt_const 
         = TOE::getConstTpetraVector(in);
     dxdt = Teuchos::rcp_const_cast<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(dxdt_const); 
   } 
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
