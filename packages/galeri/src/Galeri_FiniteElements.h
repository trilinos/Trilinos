// @HEADER
// ************************************************************************
//
//                  Galeri Matrix Generation Package
//                 Copyright (2006) Sandia Corporation
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
// ************************************************************************
// @HEADER

#include "Galeri_ConfigDefs.h"
#include "Galeri_Exception.h"
#include "FiniteElements/Galeri_Workspace.h"
// Finite element includes:
// - abstract files
#include "FiniteElements/Galeri_AbstractGrid.h"
#include "FiniteElements/Galeri_AbstractQuadrature.h"
#include "FiniteElements/Galeri_AbstractProblem.h"
#include "FiniteElements/Galeri_AbstractVariational.h"
// - grid files;
#include "FiniteElements/Galeri_FileGrid.h"
#include "FiniteElements/Galeri_TriangleRectangleGrid.h"
#include "FiniteElements/Galeri_QuadRectangleGrid.h"
#include "FiniteElements/Galeri_TetCubeGrid.h"
#include "FiniteElements/Galeri_HexCubeGrid.h"
// - quadrature files;
#include "FiniteElements/Galeri_AbstractQuadrature.h"
#include "FiniteElements/Galeri_QuadRectangleGrid.h"
#include "FiniteElements/Galeri_TriangleQuadrature.h"
#include "FiniteElements/Galeri_QuadQuadrature.h"
#include "FiniteElements/Galeri_TetQuadrature.h"
#include "FiniteElements/Galeri_HexQuadrature.h"
// - variational files
#include "FiniteElements/Galeri_SUPGVariational.h"
#include "FiniteElements/Galeri_GalerkinVariational.h"
// - problem files
#include "FiniteElements/Galeri_LinearProblem.h"
// - other files
#include "FiniteElements/Galeri_MEDITInterface.h"
