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
