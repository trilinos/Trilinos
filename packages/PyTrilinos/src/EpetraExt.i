// -*- c++ -*-

%module EpetraExt

%{
// System includes
#include <vector>

// Epetra includes
#include "Epetra_CrsGraph.h"
#include "Epetra_MapColoring.h"

// Epetraext includes
#include "EpetraExt_MapColoring.h"
#include "EpetraExt_MapColoringIndex.h"
%}

// SWIG library includes
%include "std_vector.i"

// Epetra interface import
%import "Epetra.i"

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
