// -*- c++ -*-

%module(package="NOX") Epetra

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
#include "NOX_Epetra_Interface.H"
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_FiniteDifference.H"
#include "NOX_Epetra_FiniteDifferenceColoring.H"
#include "NOX_Epetra_Vector.H"

// Local includes
#include "Epetra_NumPyVector.h"
#include "Epetra_VectorHelper.h"
#include "NumPyWrapper.h"
#include "Callback.H"
#include "PyInterface.H"
%}

// Ignore directives
%ignore *::print(ostream &, int) const;
%ignore operator=;
%ignore operator[];
%ignore NOX::Abstract::Vector::print() const;
%ignore NOX::Epetra::Vector(Epetra_Vector&, NOX::CopyType, bool);
%ignore NOX::Epetra::Vector::getEpetraVector() const;
%ignore Callback::getFunction() const;

// Rename directives
%rename(Group_None) NOX::Epetra::Group::None;

// SWIG library includes
%include "std_vector.i"

// Epetra, EpetraExt and NOX::Abstract imports
%import "RawEpetra.i"
%import "EpetraExt.i"
%import "NOX_Abstract.i"

// NOX interface includes
using namespace std;
%include "NOX_Epetra_Interface.H"
%include "NOX_Epetra_Group.H"
%include "NOX_Epetra_FiniteDifference.H"
%include "NOX_Epetra_FiniteDifferenceColoring.H"
%include "NOX_Epetra_Vector.H"

// Local interface includes
%include "Callback.H"
%include "PyInterface.H"

// Extensions
%extend NOX::Epetra::Group {
  void getSoln(PyObject * p_pyObject) {
    Epetra_VectorHelper::unloadViaCopy(& const_cast<Epetra_Vector &> 
                                    ((dynamic_cast<const NOX::Epetra::Vector &>
                                      (self->getX())).getEpetraVector()),
                                    p_pyObject);
  }
}
