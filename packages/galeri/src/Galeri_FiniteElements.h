// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
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
// Questions about Galeri? Contact Marzio Sala (marzio.sala _AT_ gmail.com)
//
// ************************************************************************
// @HEADER

#include "Galeri_ConfigDefs.h"
#include "Galeri_Exception.h"
#include "Galeri_Workspace.h"
// Finite element includes:
// - abstract files
#include "Galeri_AbstractGrid.h"
#include "Galeri_AbstractQuadrature.h"
#include "Galeri_AbstractProblem.h"
#include "Galeri_AbstractVariational.h"
// - grid files;
#include "Galeri_FileGrid.h"
#include "Galeri_TriangleRectangleGrid.h"
#include "Galeri_QuadRectangleGrid.h"
#include "Galeri_TetCubeGrid.h"
#include "Galeri_HexCubeGrid.h"
// - quadrature files;
#include "Galeri_AbstractQuadrature.h"
#include "Galeri_QuadRectangleGrid.h"
#include "Galeri_TriangleQuadrature.h"
#include "Galeri_QuadQuadrature.h"
#include "Galeri_TetQuadrature.h"
#include "Galeri_HexQuadrature.h"
// - variational files
#include "Galeri_SUPGVariational.h"
#include "Galeri_GalerkinVariational.h"
// - problem files
#include "Galeri_LinearProblem.h"
// - other files
#include "Galeri_MEDITInterface.h"
