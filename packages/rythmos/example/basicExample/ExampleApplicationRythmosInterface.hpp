//
// @HEADER
// ***********************************************************************
// 
//                           Rythmos Package
//                 Copyright (2005) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef Rythmos_EXAMPLE_APPLICATION_RYTHMOS_INTERFACE_H
#define Rythmos_EXAMPLE_APPLICATION_RYTHMOS_INTERFACE_H

#include "Epetra_Map.h"
//#include "Epetra_Vector.h"
//#include "Rythmos_ConfigDefs.h"

#include "ExampleApplication.hpp"

#include "InOutArgs.hpp"
#include "NonlinearModel.hpp"

#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"

#include "Teuchos_RefCountPtr.hpp"

//-----------------------------------------------------------------------------
// Class         : ExampleApplicationRythmosInterface
// Purpose       : Interface code to link EampleApplication to Rythmos
// Special Notes :
// Creator       : Todd Coffey, SNL
// Creation Date : 05/17/05
//-----------------------------------------------------------------------------
class ExampleApplicationRythmosInterface : public Rythmos::NonlinearModel<double>
{
  public:

    ExampleApplicationRythmosInterface();

    ~ExampleApplicationRythmosInterface();

    int evalModel(const Rythmos::InArgs<double> & inargs, const Rythmos::OutArgs<double> & outargs) const;

    Teuchos::RefCountPtr<Thyra::VectorBase<double> > get_vector() const;

    const Teuchos::RefCountPtr<const Epetra_Map> &get_Epetra_Map() const;

  protected:
    Teuchos::RefCountPtr<ExampleApplication> problem_;
    Teuchos::RefCountPtr<const Epetra_Map> epetra_map_;
    Teuchos::RefCountPtr<const Thyra::MPIVectorSpaceBase<double> > thyra_vs_;

};


#endif // Rythmos_EXAMPLE_APPLICATION_RYTHMOS_INTERFACE_H
