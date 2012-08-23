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

#ifndef __Panzer_EpetraLinearObjContainer_hpp__
#define __Panzer_EpetraLinearObjContainer_hpp__

#include "Panzer_config.hpp"

#include <map>

// Epetra includes
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "Panzer_LinearObjFactory.hpp" 
#include "Panzer_ThyraObjContainer.hpp"

#include "Teuchos_RCP.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

namespace panzer {

class EpetraLinearObjContainer : public LinearObjContainer
                               , public ThyraObjContainer<double> {

   EpetraLinearObjContainer();

public:
   typedef LinearObjContainer::Members Members;

   typedef Epetra_Vector VectorType;
   typedef Epetra_CrsMatrix CrsMatrixType;

   EpetraLinearObjContainer(const Teuchos::RCP<const Epetra_Map> & domain,
                            const Teuchos::RCP<const Epetra_Map> & range)
      : domainMap(domain), rangeMap(range) 
   {
      domainSpace = Thyra::create_VectorSpace(domainMap);
      rangeSpace = Thyra::create_VectorSpace(rangeMap);
   }

   virtual void initialize() 
   {
      if(get_x()!=Teuchos::null) get_x()->PutScalar(0.0);
      if(get_dxdt()!=Teuchos::null) get_dxdt()->PutScalar(0.0);
      if(get_f()!=Teuchos::null) get_f()->PutScalar(0.0);
      if(get_A()!=Teuchos::null) get_A()->PutScalar(0.0);
   }

   //! Wipe out stored data.
   void clear()
   {
      set_x(Teuchos::null);
      set_dxdt(Teuchos::null);
      set_f(Teuchos::null);
      set_A(Teuchos::null);
   }

   inline void set_x(const Teuchos::RCP<Epetra_Vector> & in) { x = in; } 
   inline const Teuchos::RCP<Epetra_Vector> get_x() const { return x; }

   inline void set_dxdt(const Teuchos::RCP<Epetra_Vector> & in) { dxdt = in; } 
   inline const Teuchos::RCP<Epetra_Vector> get_dxdt() const { return dxdt; }

   inline void set_f(const Teuchos::RCP<Epetra_Vector> & in) { f = in; } 
   inline const Teuchos::RCP<Epetra_Vector> get_f() const { return f; }

   inline void set_A(const Teuchos::RCP<Epetra_CrsMatrix> & in) { A = in; } 
   inline const Teuchos::RCP<Epetra_CrsMatrix> get_A() const { return A; }

   void initializeMatrix(double value)
   { A->PutScalar(value); }

   virtual void set_x_th(const Teuchos::RCP<Thyra::VectorBase<double> > & in) 
   { x = (in==Teuchos::null) ? Teuchos::null : Thyra::get_Epetra_Vector(*domainMap,in); }
   virtual Teuchos::RCP<Thyra::VectorBase<double> > get_x_th() const 
   { return (x==Teuchos::null) ? Teuchos::null : Thyra::create_Vector(x,domainSpace); }

   virtual void set_dxdt_th(const Teuchos::RCP<Thyra::VectorBase<double> > & in)
   { dxdt = (in==Teuchos::null) ? Teuchos::null : Thyra::get_Epetra_Vector(*domainMap,in); }
   virtual Teuchos::RCP<Thyra::VectorBase<double> > get_dxdt_th() const 
   { return (dxdt==Teuchos::null) ? Teuchos::null : Thyra::create_Vector(dxdt,domainSpace); }

   virtual void set_f_th(const Teuchos::RCP<Thyra::VectorBase<double> > & in)
   { f = (in==Teuchos::null) ? Teuchos::null : Thyra::get_Epetra_Vector(*rangeMap,in); }
   virtual Teuchos::RCP<Thyra::VectorBase<double> > get_f_th() const 
   { return (f==Teuchos::null) ? Teuchos::null : Thyra::create_Vector(f,rangeSpace); }

   virtual void set_A_th(const Teuchos::RCP<Thyra::LinearOpBase<double> > & in) 
   { A = (in==Teuchos::null) ? Teuchos::null : Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(Thyra::get_Epetra_Operator(*in),true); }
   virtual Teuchos::RCP<Thyra::LinearOpBase<double> > get_A_th() const
   { return (A==Teuchos::null) ? Teuchos::null : Thyra::nonconstEpetraLinearOp(A); }

private:
   Teuchos::RCP<const Epetra_Map> domainMap;
   Teuchos::RCP<const Epetra_Map> rangeMap;
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > domainSpace;
   Teuchos::RCP<const Thyra::VectorSpaceBase<double> > rangeSpace;
   Teuchos::RCP<Epetra_Vector> x, dxdt, f;
   Teuchos::RCP<Epetra_CrsMatrix> A;
};

}

#endif
