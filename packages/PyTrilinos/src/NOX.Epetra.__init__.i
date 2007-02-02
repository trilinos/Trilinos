// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//              PyTrilinos: Python Interface to Trilinos
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
// Questions? Contact Bill Spotz (wfspotz@sandia.gov)
//
// ***********************************************************************
// @HEADER

%module(package      = "PyTrilinos.NOX.Epetra",
	directors    = "1",
	autodoc      = "1",
	implicitconv = "1") __init__

%{
// System includes
#include <vector>

// Teuchos includes
#include "Teuchos_PythonParameter.h"

// Epetra includes
#include "Epetra_BLAS.h"
#include "Epetra_Object.h"
#include "Epetra_CompObject.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_DistObject.h"
#include "Epetra_IntVector.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_FEVector.h"
#include "Epetra_Operator.h"
#include "Epetra_InvOperator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_MapColoring.h"
#include "Epetra_JadMatrix.h"

// EpetraExt includes
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"

// NOX includes
#include "NOX_Abstract_Group.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_Vector.H"
#include "NOX_Epetra_FiniteDifference.H"
#include "NOX_Epetra_FiniteDifferenceColoring.H"
#include "NOX_Epetra_MatrixFree.H"
#include "NOX_Epetra_Scaling.H"
#include "NOX_Epetra_LinearSystem.H"
#include "NOX_Epetra_LinearSystem_AztecOO.H"

// Local includes
#include "NumPyImporter.h"
#include "Epetra_NumPyMultiVector.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_NumPyIntVector.h"

// Namespace flattening
using namespace NOX;
using namespace NOX::Abstract;
using namespace NOX::Epetra;
%}

// Define macro for handling exceptions thrown by NOX.Epetra methods and
// constructors
%define %nox_epetra_exception(className,methodName)
  %exception NOX::Epetra::className::methodName {
  try {
    $action
    if (PyErr_Occurred()) SWIG_fail;
  }
  catch(Teuchos::Exceptions::InvalidParameterType e) {
    PyErr_SetString(PyExc_TypeError, e.what());
    SWIG_fail;
  }
  catch(Teuchos::Exceptions::InvalidParameter e) {
    PyErr_SetString(PyExc_KeyError, e.what());
    SWIG_fail;
  }
  catch(std::runtime_error e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    SWIG_fail;
  }
  catch(std::logic_error e) {
    PyErr_SetString(PyExc_RuntimeError, e.what());
    SWIG_fail;
  }
}
%enddef

// General ignore directives
%ignore *::print(std::ostream &) const;
%ignore *::print(std::ostream &, int) const;
%ignore *::operator=;
%ignore *::operator<<;

// SWIG library includes
using namespace std;
%include "stl.i"

// Trilinos interface support
%import "Teuchos.i"

// Support for Teuchos::RefCountPtrs
%teuchos_rcp_typemaps(NOX::Epetra::LinearSystem)
%teuchos_rcp_typemaps(NOX::Epetra::Scaling)
%teuchos_rcp_typemaps(Epetra_CrsGraph)
%teuchos_rcp_typemaps(Epetra_MapColoring)
%teuchos_rcp_typemaps(Epetra_Operator)

//////////////
// Typemaps //
//////////////

// Make Epetra_Vector and NOX::Epetra::Vector input arguments
// interchangeable
%typemap(in) NOX::Epetra::Vector & (void* argp=0, int res=0) {
  res = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
  if (!SWIG_IsOK(res)) {
    res = SWIG_ConvertPtr($input, &argp, $descriptor(Epetra_Vector*), %convertptr_flags);
    if (!SWIG_IsOK(res)) {
      %argument_fail(res, "$type", $symname, $argnum);
    }
    if (!argp) { %argument_nullref("$type", $symname, $argnum); }
    $1 = new NOX::Epetra::Vector(*%reinterpret_cast( argp, Epetra_Vector*));
  } else {
    $1 = %reinterpret_cast(argp, NOX::Epetra::Vector*);
  }
}
%typecheck(1190) NOX::Epetra::Vector & {
  static void * argp = 0;
  $1 = SWIG_CheckState(SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags)) ? 1 : 0;
  if (!$1)
    $1 = SWIG_CheckState(SWIG_ConvertPtr($input, &argp, $descriptor(Epetra_Vector*),
                                         %convertptr_flags)) ? 1 : 0;
}
%typemap(freearg) NOX::Epetra::Vector {
  delete $1;
}

// Convert NOX::Abstract::Vector return arguments to Epetra.Vectors
%typemap(out) NOX::Abstract::Vector & {
  static swig_type_info * swig_ENPV_ptr = SWIG_TypeQuery("Epetra_NumPyVector*");
  NOX::Epetra::Vector * nevResult = dynamic_cast<NOX::Epetra::Vector*>($1);
  if (nevResult == NULL) {
    // If we can't upcast, then return the NOX::Abstract::Vector
    $result = SWIG_NewPointerObj((void*)&$1, $descriptor, 1);
  } else {
    Epetra_NumPyVector * enpvResult = new Epetra_NumPyVector(View, nevResult->getEpetraVector(), 0);
    $result = SWIG_NewPointerObj((void*)enpvResult, swig_ENPV_ptr, 1);
  }
}

// Convert Epetra_Vector return arguments to Epetra.Vectors
%typemap(out) Epetra_Vector & {
  static swig_type_info * swig_ENPV_ptr = SWIG_TypeQuery("Epetra_NumPyVector*");
  Epetra_NumPyVector * enpvResult = new Epetra_NumPyVector(View, *$1, 0);
  $result = SWIG_NewPointerObj((void*)enpvResult, swig_ENPV_ptr, 1);
}

// Epetra imports
%import "Epetra_SrcDistObject.h"
%import "Epetra_Operator.h"
%import "Epetra_RowMatrix.h"

// NOX imports
%import "NOX_Abstract_Group.H"
%import "NOX_Abstract_Vector.H"

// NOX::Epetra::Interface imports
%import "NOX.Epetra.Interface.i"

//////////////////////////////
// NOX.Epetra.Group support //
//////////////////////////////
%nox_epetra_exception(Group,Group)
%rename(Group_None) NOX::Epetra::Group::None;
%include "NOX_Epetra_Group.H"

///////////////////////////////
// NOX.Epetra.Vector support //
///////////////////////////////
%ignore NOX::Epetra::Vector(Epetra_Vector&, NOX::CopyType, bool);
%ignore NOX::Epetra::Vector::getEpetraVector() const;
%include "NOX_Epetra_Vector.H"

/////////////////////////////////////////
// NOX.Epetra.FiniteDifference support //
/////////////////////////////////////////
%nox_epetra_exception(FiniteDifference,FiniteDifference)
%include "NOX_Epetra_FiniteDifference.H"

/////////////////////////////////////////////////
// NOX.Epetra.FiniteDifferenceColoring support //
/////////////////////////////////////////////////
%nox_epetra_exception(FiniteDifferenceColoring,FiniteDifferenceColoring)
namespace NOX {
  namespace Epetra {
    %extend FiniteDifferenceColoring {
      FiniteDifferenceColoring(Teuchos::ParameterList & printingParams,
			       Interface::Required & i,
			       const NOX::Epetra::Vector & initialGuess,
			       Epetra_CrsGraph & rawGraph,
			       bool parallelColoring = false,
			       bool distance1 = false,
			       double beta = 1.0e-6, double alpha = 1.0e-4) {
	// Wrap the interface and CRS graph in reference counters
	Teuchos::RefCountPtr<Interface::Required> i_ptr     = Teuchos::rcp(&i,        false);
	Teuchos::RefCountPtr<Epetra_CrsGraph>     graph_ptr = Teuchos::rcp(&rawGraph, false);
	// Construct the coloring algorithm functor and generate the color map
	EpetraExt::CrsGraph_MapColoring          *mapColor =
	  new EpetraExt::CrsGraph_MapColoring();
	Teuchos::RefCountPtr<Epetra_MapColoring> colorMap  = 
	  Teuchos::rcp(&(*mapColor)(rawGraph));
	// Construct the color index functor and generate the column indexes
	EpetraExt::CrsGraph_MapColoringIndex                 *mapColorIndex =
	  new EpetraExt::CrsGraph_MapColoringIndex(*colorMap);
	Teuchos::RefCountPtr<std::vector<Epetra_IntVector> > columns        =
	  Teuchos::rcp(&(*mapColorIndex)(rawGraph));
	// Construct the FiniteDifferenceColoring object
	FiniteDifferenceColoring *fdc = new FiniteDifferenceColoring(printingParams,
								     i_ptr,
								     initialGuess,
								     graph_ptr,
								     colorMap, columns,
								     parallelColoring,
								     distance1, beta,
								     alpha);
	// Delete the temporary functors
	delete mapColor;
	delete mapColorIndex;
	// Return the pointer to FiniteDifferenceColoring object
	return fdc;
      }
    }
    %ignore FiniteDifferenceColoring::FiniteDifferenceColoring;
  }
}
%include "NOX_Epetra_FiniteDifferenceColoring.H"

///////////////////////////////////
// NOX.Epetra.MatrixFree support //
///////////////////////////////////
%nox_epetra_exception(MatrixFree,MatrixFree)
%include "NOX_Epetra_MatrixFree.H"

////////////////////////////////
// NOX.Epetra.Scaling support //
////////////////////////////////
%ignore operator<<(ostream&, NOX::Epetra::Scaling&);
%rename(Scaling_None) NOX::Epetra::Scaling::None;
%include "NOX_Epetra_Scaling.H"

/////////////////////////////////////
// NOX.Epetra.LinearSystem support //
/////////////////////////////////////
%nox_epetra_exception(LinearSystem,LinearSystem)
%feature("director") NOX::Epetra::LinearSystem;
// The following #define is to change the name of NOX method
// arguments that conflict with a SWIG director method argument
#define result nox_result
%include "NOX_Epetra_LinearSystem.H"

////////////////////////////////////////////
// NOX.Epetra.LinearSystemAztecOO support //
////////////////////////////////////////////
%nox_epetra_exception(LinearSystemAztecOO,LinearSystemAztecOO)
%include "NOX_Epetra_LinearSystem_AztecOO.H"
