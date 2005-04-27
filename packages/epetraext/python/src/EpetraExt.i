// -*- c++ -*-

// @HEADER
// ***********************************************************************
//
//         PyTrilinos.EpetraExt: Python Interface to EpetraExt
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

%module(package="PyTrilinos") EpetraExt

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
#include "EpetraExt_Version.h"
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

// Epetra interface import
%import "Epetra_Object.h"
%import "Epetra_SrcDistObject.h"
%import "Epetra_DistObject.h"
%import "Epetra_CrsGraph.h"
%import "Epetra_IntVector.h"
%import "Epetra_MapColoring.h"

// C++ STL support
%include "std_string.i"
%include "std_vector.i"

// EpetraExt interface includes
%include "EpetraExt_Version.h"
%include "EpetraExt_Transform.h"
%template () std::vector<Epetra_IntVector>;
%template () EpetraExt::Transform<Epetra_CrsGraph, Epetra_MapColoring>;
%template () EpetraExt::Transform<Epetra_CrsGraph, std::vector<Epetra_IntVector,
							       std::allocator<Epetra_IntVector> > >;
//%template () EpetraExt::Transform<Epetra_CrsGraph, std::vector<Epetra_IntVector> >;
%template () EpetraExt::StructuralTransform<Epetra_CrsGraph, Epetra_MapColoring>;
%template () EpetraExt::StructuralTransform<Epetra_CrsGraph, std::vector<Epetra_IntVector> >;

%include "EpetraExt_MapColoring.h"
%include "EpetraExt_MapColoringIndex.h"

// Python code.
%pythoncode %{

  __version__ = EpetraExt_Version().split()[2]

%}
