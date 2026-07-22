// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ThyraObjContainer_hpp__
#define __Panzer_ThyraObjContainer_hpp__

#include "PanzerDiscFE_config.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace Thyra {
   template <typename ScalarT> class VectorBase;
   template <typename ScalarT> class LinearOpBase;
}

namespace panzer {

template <typename ScalarT>
class ThyraObjContainer {
public:
   virtual ~ThyraObjContainer() {}

   //! Put a particular scalar in the matrix
   virtual void initializeMatrix(ScalarT value) = 0;

   virtual void set_x_th(const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & in) = 0;
   virtual Teuchos::RCP<Thyra::VectorBase<ScalarT> > get_x_th() const = 0;

   virtual void set_dxdt_th(const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & in) = 0;
   virtual Teuchos::RCP<Thyra::VectorBase<ScalarT> > get_dxdt_th() const = 0;

   virtual void set_f_th(const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & in) = 0;
   virtual Teuchos::RCP<Thyra::VectorBase<ScalarT> > get_f_th() const = 0;

   virtual void set_A_th(const Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > & in) = 0;
   virtual Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > get_A_th() const = 0;

   void clear()
   {
     set_x_th(Teuchos::null);
     set_dxdt_th(Teuchos::null);
     set_f_th(Teuchos::null);
     set_A_th(Teuchos::null);
   }
};

}

#endif
