// -*- c++ -*-

%module NOX_Epetra

%{
// Epetra includes
#include "Epetra_LocalMap.h"
#include "Epetra_SerialComm.h"
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_RowMatrix.h"

// NOX includes
#include "NOX_Abstract_Group.H"
#include "NOX_Abstract_Vector.H"
#include "NOX_Epetra_Interface.H"
#include "NOX_Epetra_Group.H"
#include "NOX_Epetra_FiniteDifference.H"
#include "NOX_Epetra_FiniteDifferenceColoring.H"
#include "NOX_Epetra_Vector.H"

// Local includes
#include "Callback.h"
#include "Epetra_NumPyVector.h"
#include "Epetra_VectorHelper.h"
#include "NumPyWrapper.h"
#include "PyInterface.h"

// Namespace flattening
using namespace NOX          ;
using namespace NOX::Abstract;
using namespace NOX::Epetra  ;
%}

// Ignore directives
%ignore *::print(ostream &, int) const;
%ignore Epetra_CompObject::operator=(const Epetra_CompObject &);
%ignore NOX::Abstract::Group;
%ignore NOX::Abstract::Group::operator=(const NOX::Abstract::Group&);
%ignore NOX::Abstract::Group::operator=(const NOX::Epetra::Group&);
%ignore NOX::Abstract::Vector;
%ignore NOX::Abstract::Vector(Epetra_Vector, NOX::CopyType, bool);
%ignore NOX::Abstract::Vector::operator=(const Epetra_Vector&);
%ignore NOX::Abstract::Vector::operator=(const NOX::Epetra::Vector&);
%ignore NOX::Abstract::Vector::operator=(const NOX::Abstract::Vector&);
%ignore NOX::Abstract::Vector::print() const;
%ignore NOX::Epetra::Group::operator=(const NOX::Epetra::Group&);
%ignore NOX::Epetra::Group::operator=(const NOX::Abstract::Group&);
%ignore NOX::Epetra::Vector(Epetra_Vector&, NOX::CopyType, bool);

// Rename directives
%rename(Group_None) NOX::Epetra::Group::None;

// Epetra interface import
%import "RawEpetra.i"

// NOX interface includes
using namespace std;
%include "NOX_Abstract_Group.H"
%include "NOX_Abstract_Vector.H"
%include "NOX_Epetra_Interface.H"
%include "NOX_Epetra_Group.H"
%include "NOX_Epetra_FiniteDifference.H"
%include "NOX_Epetra_FiniteDifferenceColoring.H"
%include "NOX_Epetra_Vector.H"

// Local interface includes
%include "Callback.h"
%include "PyInterface.h"

// Extensions
%extend NOX::Epetra::Group {
  void getSoln(PyObject * p_pyObject) {
    Epetra_VectorHelper::unloadViaCopy(& const_cast<Epetra_Vector &> 
				       ((dynamic_cast<const NOX::Epetra::Vector &>
					 (self->getX())).getEpetraVector()),
				       p_pyObject);
  }
}
