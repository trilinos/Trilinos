// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//               PyTrilinos.NOX: Python Interface to NOX
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

%module(package="PyTrilinos.NOX") Epetra

%{
// Epetra includes
#include "Epetra_BLAS.h"
#include "Epetra_Object.h"
#include "Epetra_CompObject.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_DistObject.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_MapColoring.h"

// NOX includes
#include "NOX_Abstract_Group.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Epetra_Interface_Required.H"
#include "NOX_Epetra_Interface_Jacobian.H"
#include "NOX_Epetra_Interface_Preconditioner.H"
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_FiniteDifference.H"
#include "NOX_Epetra_FiniteDifferenceColoring.H"
#include "NOX_Epetra_Vector.H"

// Local includes
#include "NumPyImporter.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_VectorHelper.h"
#include "Callback.h"
#include "PyInterface.h"

// Namespace flattening
using namespace NOX          ;
using namespace NOX::Abstract;
using namespace NOX::Epetra  ;
%}

// Ignore directives
%ignore *::print(ostream &, int) const;
%ignore Epetra_CompObject::operator=(const Epetra_CompObject &);
%ignore Epetra_IntVector::operator=(const Epetra_IntVector &);
%ignore Epetra_IntVector::operator[](int);
%ignore Epetra_IntVector::operator[](int) const;
%ignore Epetra_MultiVector::operator=(const Epetra_MultiVector &);
%ignore Epetra_MultiVector::operator[](int);
%ignore Epetra_MultiVector::operator[](int) const;
%ignore Epetra_CrsGraph::operator[](int);
%ignore Epetra_CrsGraph::operator[](int) const;
%ignore Epetra_CrsGraph::operator=(const Epetra_CrsGraph&);
%ignore Epetra_MapColoring::operator[](int);
%ignore Epetra_MapColoring::operator[](int) const;
%ignore NOX::Abstract::Group::operator=(const NOX::Abstract::Group&);
%ignore NOX::Abstract::Vector::operator=(const Epetra_Vector&);
%ignore NOX::Abstract::Vector::operator=(const NOX::Epetra::Vector&);
%ignore NOX::Abstract::Vector::operator=(const NOX::Abstract::Vector&);
%ignore NOX::Abstract::Vector::print() const;
%ignore NOX::Epetra::Group::operator=(const NOX::Epetra::Group&);
%ignore NOX::Epetra::Group::operator=(const NOX::Abstract::Group&);
%ignore NOX::Epetra::Vector(Epetra_Vector&, NOX::CopyType, bool);
%ignore NOX::Epetra::Vector::getEpetraVector() const;
%ignore Callback::getFunction() const;

// Rename directives
%rename(Group_None) NOX::Epetra::Group::None;

// SWIG library includes
%include "std_vector.i"

// Epetra, EpetraExt and NOX::Abstract imports
%import "Epetra_BLAS.h"
%import "Epetra_Object.h"
%import "Epetra_CompObject.h"
%import "Epetra_SrcDistObject.h"
%import "Epetra_DistObject.h"
%import "Epetra_IntVector.h"
%import "Epetra_MultiVector.h"
%import "Epetra_Vector.h"
%import "Epetra_Operator.h"
%import "Epetra_RowMatrix.h"
%import "Epetra_CrsGraph.h"
%import "Epetra_MapColoring.h"
%import "NOX_Abstract_Group.H"
%import "NOX_Abstract_Vector.H"
%import "Epetra_NumPyVector.h"

// NOX interface includes
using namespace std;
%include "NOX_Epetra_Interface_Required.H"
// %include "NOX_Epetra_Interface_Jacobian.H"
// %include "NOX_Epetra_Interface_Preconditioner.H"
// using namespace NOX::Epetra::Interface;
// %include "NOX_Epetra_Group.H"
// %include "NOX_Epetra_FiniteDifference.H"
// %include "NOX_Epetra_FiniteDifferenceColoring.H"
// %include "NOX_Epetra_Vector.H"

// // Local interface includes
// %include "Callback.h"
// %include "PyInterface.h"

// // Extensions
// %extend NOX::Epetra::Group {
//   void getSoln(PyObject * p_pyObject) {
//     Epetra_VectorHelper::unloadViaCopy(& const_cast<Epetra_Vector &> 
//                                     ((dynamic_cast<const NOX::Epetra::Vector &>
//                                       (self->getX())).getEpetraVector()),
//                                     p_pyObject);
//   }
// }
