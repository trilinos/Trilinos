// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Thyra_VectorStdOps.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_TpetraLinearOp.hpp"

#include "Tpetra_CrsMatrix.hpp"

namespace panzer {

//! Make sure row and column spaces match up
template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
bool BlockedTpetraLinearObjContainer<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
checkCompatibility() const
{
   using Thyra::VectorSpaceBase;
   using Teuchos::RCP;
   using Teuchos::null;

   bool x_matches=false, f_matches=false, dxdt_matches=false;

   if(get_A()!=null) {
      RCP<const VectorSpaceBase<ScalarT> > range  = get_A()->range();
      RCP<const VectorSpaceBase<ScalarT> > domain = get_A()->domain();

      if(get_x()!=null)
         x_matches = range->isCompatible(*get_x()->space());
      else
         x_matches = true; // nothing to compare

      if(get_dxdt()!=null)
         dxdt_matches = range->isCompatible(*get_dxdt()->space());
      else
         dxdt_matches = true; // nothing to compare

      if(get_f()!=null)
         f_matches = range->isCompatible(*get_f()->space());
      else
         f_matches = true; // nothing to compare
   }
   else if(get_x()!=null && get_dxdt()!=null) {
      f_matches = true; // nothing to compare f to
      x_matches = get_x()->space()->isCompatible(*get_dxdt()->space());  // dxdt and x are in the same space
      dxdt_matches = x_matches;
   }
   else {
      f_matches = x_matches = dxdt_matches = true; // nothing to compare to
   }

   return x_matches && dxdt_matches && f_matches;
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjContainer<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
initialize()
{
   using Thyra::LinearOpBase;
   using Thyra::PhysicallyBlockedLinearOpBase;
   using Thyra::ProductVectorSpaceBase;
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   if(get_x()!=Teuchos::null)    Thyra::assign<ScalarT>(x.ptr(),0.0);
   if(get_dxdt()!=Teuchos::null) Thyra::assign<ScalarT>(get_dxdt().ptr(),0.0);
   if(get_f()!=Teuchos::null)    Thyra::assign<ScalarT>(get_f().ptr(),0.0);
   if(get_A()!=Teuchos::null) {
      RCP<PhysicallyBlockedLinearOpBase<ScalarT> > Amat
            = rcp_dynamic_cast<PhysicallyBlockedLinearOpBase<ScalarT> >(get_A(),true);
      RCP<const ProductVectorSpaceBase<ScalarT> > range = Amat->productRange();
      RCP<const ProductVectorSpaceBase<ScalarT> > domain = Amat->productDomain();

      // loop over block entries
      for(int i=0;i<range->numBlocks();i++) {
         for(int j=0;j<domain->numBlocks();j++) {
            RCP<LinearOpBase<ScalarT> > block = Amat->getNonconstBlock(i,j);
            if(block!=Teuchos::null) {
               RCP<Tpetra::Operator<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > t_block =
                   rcp_dynamic_cast<Thyra::TpetraLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(block,true)->getTpetraOperator();

               RCP<const MapType> map_i = t_block->getRangeMap();
               RCP<const MapType> map_j = t_block->getDomainMap();

               RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > mat =
                   rcp_dynamic_cast<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(t_block,true);

               mat->resumeFill();
               mat->setAllToScalar(0.0);
               mat->fillComplete(map_j,map_i);
            }
         }
      }
   }
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjContainer<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
initializeMatrix(ScalarT value)
{
   using Thyra::LinearOpBase;
   using Thyra::PhysicallyBlockedLinearOpBase;
   using Thyra::ProductVectorSpaceBase;
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   if(get_A()!=Teuchos::null) {
      RCP<PhysicallyBlockedLinearOpBase<ScalarT> > Amat
            = rcp_dynamic_cast<PhysicallyBlockedLinearOpBase<ScalarT> >(get_A(),true);
      RCP<const ProductVectorSpaceBase<ScalarT> > range = Amat->productRange();
      RCP<const ProductVectorSpaceBase<ScalarT> > domain = Amat->productDomain();

      // loop over block entries
      for(int i=0;i<range->numBlocks();i++) {
         for(int j=0;j<domain->numBlocks();j++) {
            RCP<LinearOpBase<ScalarT> > block = Amat->getNonconstBlock(i,j);
            if(block!=Teuchos::null) {
               RCP<Tpetra::Operator<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > t_block =
                   rcp_dynamic_cast<Thyra::TpetraLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(block,true)->getTpetraOperator();

               // why do I have to do this?
               RCP<const MapType> map_i = t_block->getRangeMap();
               RCP<const MapType> map_j = t_block->getDomainMap();

               RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > mat =
                   rcp_dynamic_cast<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(t_block,true);

               mat->resumeFill();
               mat->setAllToScalar(value);
               mat->fillComplete(map_j,map_i);
            }
         }
      }
   }
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjContainer<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
clear()
{
   set_x(Teuchos::null);
   set_dxdt(Teuchos::null);
   set_f(Teuchos::null);
   set_A(Teuchos::null);
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjContainer<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
beginFill()
{
   using Thyra::LinearOpBase;
   using Thyra::PhysicallyBlockedLinearOpBase;
   using Thyra::ProductVectorSpaceBase;
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   if(get_A()!=Teuchos::null) {
      RCP<PhysicallyBlockedLinearOpBase<ScalarT> > Amat
            = rcp_dynamic_cast<PhysicallyBlockedLinearOpBase<ScalarT> >(get_A(),true);
      RCP<const ProductVectorSpaceBase<ScalarT> > range = Amat->productRange();
      RCP<const ProductVectorSpaceBase<ScalarT> > domain = Amat->productDomain();

      // loop over block entries
      for(int i=0;i<range->numBlocks();i++) {
         for(int j=0;j<domain->numBlocks();j++) {
            RCP<LinearOpBase<ScalarT> > block = Amat->getNonconstBlock(i,j);
            if(block!=Teuchos::null) {
               RCP<Tpetra::Operator<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > t_block =
                   rcp_dynamic_cast<Thyra::TpetraLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(block,true)->getTpetraOperator();

               RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > mat =
                   rcp_dynamic_cast<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(t_block,true);

               mat->resumeFill();
            }
         }
      }
   }
}

template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjContainer<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
endFill()
{
   using Thyra::LinearOpBase;
   using Thyra::PhysicallyBlockedLinearOpBase;
   using Thyra::ProductVectorSpaceBase;
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   if(get_A()!=Teuchos::null) {
      RCP<PhysicallyBlockedLinearOpBase<ScalarT> > Amat
            = rcp_dynamic_cast<PhysicallyBlockedLinearOpBase<ScalarT> >(get_A(),true);
      RCP<const ProductVectorSpaceBase<ScalarT> > range = Amat->productRange();
      RCP<const ProductVectorSpaceBase<ScalarT> > domain = Amat->productDomain();

      // loop over block entries
      for(int i=0;i<range->numBlocks();i++) {
         for(int j=0;j<domain->numBlocks();j++) {
            RCP<LinearOpBase<ScalarT> > block = Amat->getNonconstBlock(i,j);
            if(block!=Teuchos::null) {
               RCP<Tpetra::Operator<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > t_block =
                   rcp_dynamic_cast<Thyra::TpetraLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(block,true)->getTpetraOperator();

               RCP<const MapType> map_i = t_block->getRangeMap();
               RCP<const MapType> map_j = t_block->getDomainMap();

               RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > mat =
                   rcp_dynamic_cast<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(t_block,true);

               mat->fillComplete(map_j,map_i);
            }
         }
      }
   }
}

}
