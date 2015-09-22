// @HEADER
// ************************************************************************
//
//           Galeri: Finite Element and Matrix Generation Package
//                 Copyright (2006) ETHZ/Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
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
