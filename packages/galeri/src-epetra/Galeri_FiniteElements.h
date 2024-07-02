// @HEADER
// *****************************************************************************
//           Galeri: Finite Element and Matrix Generation Package
//
// Copyright 2006 ETHZ/NTESS and the Galeri contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
