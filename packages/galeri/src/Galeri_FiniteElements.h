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
