// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_BlockedEpetraLinearObjContainer.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_ProductVectorSpaceBase.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

#include "Epetra_CrsMatrix.h"

namespace panzer {

//! Make sure row and column spaces match up
bool BlockedEpetraLinearObjContainer::
checkCompatibility() const
{
   using Thyra::VectorSpaceBase;
   using Teuchos::RCP;
   using Teuchos::null;

   bool x_matches=false, f_matches=false, dxdt_matches=false;

   if(get_A()!=null) {
      RCP<const VectorSpaceBase<double> > range  = get_A()->range();
      RCP<const VectorSpaceBase<double> > domain = get_A()->domain();

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

void BlockedEpetraLinearObjContainer::
initialize()
{
   using Thyra::LinearOpBase;
   using Thyra::PhysicallyBlockedLinearOpBase;
   using Thyra::ProductVectorSpaceBase;
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   if(get_x()!=Teuchos::null)    Thyra::assign<double>(x.ptr(),0.0);
   if(get_dxdt()!=Teuchos::null) Thyra::assign<double>(get_dxdt().ptr(),0.0);
   if(get_f()!=Teuchos::null)    Thyra::assign<double>(get_f().ptr(),0.0);
   if(get_A()!=Teuchos::null) {
      RCP<PhysicallyBlockedLinearOpBase<double> > Amat
            = rcp_dynamic_cast<PhysicallyBlockedLinearOpBase<double> >(get_A(),true);
      RCP<const ProductVectorSpaceBase<double> > range = Amat->productRange();
      RCP<const ProductVectorSpaceBase<double> > domain = Amat->productDomain();

      // loop over block entries
      for(int i=0;i<range->numBlocks();i++) {
         for(int j=0;j<domain->numBlocks();j++) {
            RCP<LinearOpBase<double> > block = Amat->getNonconstBlock(i,j);
            if(block!=Teuchos::null) {
               RCP<Epetra_Operator> e_block = Thyra::get_Epetra_Operator(*block);
               rcp_dynamic_cast<Epetra_CrsMatrix>(e_block,true)->PutScalar(0.0);
            }
         }
      }
   }
}

void BlockedEpetraLinearObjContainer::
initializeMatrix(double value)
{
   using Thyra::LinearOpBase;
   using Thyra::PhysicallyBlockedLinearOpBase;
   using Thyra::ProductVectorSpaceBase;
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   if(get_A()!=Teuchos::null) {
      RCP<PhysicallyBlockedLinearOpBase<double> > Amat
            = rcp_dynamic_cast<PhysicallyBlockedLinearOpBase<double> >(get_A(),true);
      RCP<const ProductVectorSpaceBase<double> > range = Amat->productRange();
      RCP<const ProductVectorSpaceBase<double> > domain = Amat->productDomain();

      // loop over block entries
      for(int i=0;i<range->numBlocks();i++) {
         for(int j=0;j<domain->numBlocks();j++) {
            RCP<LinearOpBase<double> > block = Amat->getNonconstBlock(i,j);
            if(block!=Teuchos::null) {
               RCP<Epetra_Operator> e_block = Thyra::get_Epetra_Operator(*block);
               rcp_dynamic_cast<Epetra_CrsMatrix>(e_block,true)->PutScalar(value);
            }
         }
      }
   }
}

void BlockedEpetraLinearObjContainer::
clear()
{
   set_x(Teuchos::null);
   set_dxdt(Teuchos::null);
   set_f(Teuchos::null);
   set_A(Teuchos::null);
}

}