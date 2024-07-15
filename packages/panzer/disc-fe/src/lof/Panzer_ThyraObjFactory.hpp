// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_ThyraObjFactory_hpp__
#define __Panzer_ThyraObjFactory_hpp__

#include "PanzerDiscFE_config.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_dyn_cast.hpp"

namespace Thyra {
   template <typename ScalarT> class VectorSpaceBase;
   template <typename ScalarT> class LinearOpBase;
}

namespace panzer {

template <typename ScalarT>
class ThyraObjFactory {
public:
   virtual ~ThyraObjFactory() {}
   
   //! Get the domain space
   virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > getThyraDomainSpace() const = 0;

   //! Get the range space
   virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > getThyraRangeSpace() const = 0;

   //! Get a matrix operator
   virtual Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > getThyraMatrix() const = 0;

};

}

#endif
