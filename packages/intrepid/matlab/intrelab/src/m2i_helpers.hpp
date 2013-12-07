#include "mex.h"
#include "Intrepid_CellTools.hpp"
#include "Intrepid_HGRAD_LINE_C1_FEM.hpp"
#include "Intrepid_HGRAD_LINE_Cn_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C1_FEM.hpp"
#include "Intrepid_HGRAD_TRI_C2_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C1_FEM.hpp"
#include "Intrepid_HGRAD_QUAD_C2_FEM.hpp"
#include "Intrepid_HGRAD_TET_C1_FEM.hpp"
#include "Intrepid_HGRAD_TET_C2_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C1_FEM.hpp"
#include "Intrepid_HGRAD_HEX_C2_FEM.hpp"

using namespace Intrepid;

void m2iGetArrayDims(Teuchos::Array<int> &dim_array, const mxArray* a, bool strip_last=false);

void m2iGetCellTopo(Teuchos::RCP<shards::CellTopology> &cellTopo, const std::string &cell_type, std::string &descriptor);

EOperator m2iGetDiffOperator(const std::string &op, std::string &descriptor);

void m2iGetBasis(Teuchos::RCP<Intrepid::Basis<double,FieldContainer<double> > > &basis,
                 const std::string &cell_type,
                 int degree,
                 std::string &descriptor);

