// -*- c++ -*-

%module EpetraExt

%{
// System includes
#include <vector>

// Epetra includes
#include "Epetra_Object.h"
#include "Epetra_SrcDistObject.h"
#include "Epetra_DistObject.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_IntVector.h"
#include "Epetra_MapColoring.h"

// EpetraExt includes
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
%}

// Ignore directives
%ignore Epetra_CrsGraph::operator[](int);
%ignore Epetra_CrsGraph::operator[](int) const;
%ignore Epetra_CrsGraph::operator=(const Epetra_CrsGraph &);
%ignore Epetra_IntVector::operator=(const Epetra_IntVector &);
%ignore Epetra_IntVector::operator[](int);
%ignore Epetra_IntVector::operator[](int) const;
%ignore Epetra_MapColoring::operator[](int);
%ignore Epetra_MapColoring::operator[](int) const;

// SWIG library includes
%include "std_vector.i"

// Epetra interface import
//%import "RawEpetra.i"
%import "Epetra_Object.h"
%import "Epetra_SrcDistObject.h"
%import "Epetra_DistObject.h"
%import "Epetra_CrsGraph.h"
%import "Epetra_IntVector.h"
%import "Epetra_MapColoring.h"

// Epetra interface includes
%include "EpetraExt_Transform.h"
%template (Transform_CrsGraph_MapColoring) EpetraExt::Transform<
  Epetra_CrsGraph, Epetra_MapColoring>;
%template (Transform_CrsGraph_MapColoringIndex) EpetraExt::Transform<
  Epetra_CrsGraph, std::vector<Epetra_IntVector> >;
%template (StructuralTransform_CrsGraph_MapColoring) EpetraExt::StructuralTransform<
  Epetra_CrsGraph, Epetra_MapColoring>;
%template (StructuralTransform_CrsGraph_MapColoringIndex)
  EpetraExt::StructuralTransform<Epetra_CrsGraph, std::vector<Epetra_IntVector> >;

%include "EpetraExt_MapColoring.h"
%include "EpetraExt_MapColoringIndex.h"
