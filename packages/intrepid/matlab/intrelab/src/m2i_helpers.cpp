#include "mex.h"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void m2iGetArrayDims(Teuchos::Array<int> &dim_array, const mxArray* a, bool strip_last)
{
    // Get the size of the input mxArray
    int number_of_dims = mxGetNumberOfDimensions(a);
    // Get pointer to dimensions array
    //const int* dim_ptr = mxGetDimensions(a);
    const mwSize* dim_ptr = mxGetDimensions(a);
    // Resize and fill Teuchos array
    dim_array.resize(number_of_dims);
    for(int i = 0; i < number_of_dims; i++)
    {
        dim_array[number_of_dims - i - 1] = *dim_ptr++;
    }
    if(strip_last && (number_of_dims == 2) && (dim_array[1] == 1))
    {
        dim_array.resize(1);
    }
}

void m2iGetCellTopo(Teuchos::RCP<shards::CellTopology> &cellTopo, const std::string &cell_type, std::string &descriptor)
{
    // Create cell topology
    if(cell_type == "Line")
    {
        cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData<shards::Line<> >()));
    }
    else if(cell_type == "Triangle")
    {
        cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData<shards::Triangle<> >()));
    }
    else if(cell_type == "Quadrilateral")
    {
        cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<> >()));
    }
    else if(cell_type == "Tetrahedron")
    {
        cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<> >()));
    }
    else if(cell_type == "Hexahedron")
    {
        cellTopo = Teuchos::rcp(new shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<> >()));
    }
    else
    {
        std::string ioError = descriptor + "Unknown cell type!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
}

EOperator m2iGetDiffOperator(const std::string &op, std::string &descriptor)
{
    EOperator retOp = OPERATOR_MAX;
    // Create Intrepid differential operator
    if(op == "OPERATOR_VALUE")
    {
        retOp = OPERATOR_VALUE;
    }
    else if(op == "OPERATOR_GRAD")
    {
        retOp = OPERATOR_GRAD;
    }
    else
    {
        std::string ioError = descriptor + "Unknown differential operator!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
    return retOp;
}

void m2iGetBasis(Teuchos::RCP<Intrepid::Basis<double, FieldContainer<double> > > &basis,
const std::string &cell_type,
int order,
std::string &descriptor)
{
    // Create basis
    if(cell_type == "Line")
    {
        if(order == 1)
        {
            basis = Teuchos::rcp(new Basis_HGRAD_LINE_C1_FEM<double, FieldContainer<double> >);
        }
        else if(order == 2)
        {
            basis = Teuchos::rcp(new Basis_HGRAD_LINE_Cn_FEM<double, FieldContainer<double> >(2, POINTTYPE_EQUISPACED));
        }
    }
    else if(cell_type == "Triangle")
    {
        if(order == 1)
        {
            basis = Teuchos::rcp(new Basis_HGRAD_TRI_C1_FEM<double, FieldContainer<double> >);
        }
        else if(order == 2)
        {
            basis = Teuchos::rcp(new Basis_HGRAD_TRI_C2_FEM<double, FieldContainer<double> >);
        }
    }
    else if(cell_type == "Quadrilateral")
    {
        if(order == 1)
        {
            basis = Teuchos::rcp(new Basis_HGRAD_QUAD_C1_FEM<double, FieldContainer<double> >);
        }
        else if(order == 2)
        {
            basis = Teuchos::rcp(new Basis_HGRAD_QUAD_C2_FEM<double, FieldContainer<double> >);
        }
    }
    else if(cell_type == "Tetrahedron")
    {
        if(order == 1)
        {
            basis = Teuchos::rcp(new Basis_HGRAD_TET_C1_FEM<double, FieldContainer<double> >);
        }
        else if(order == 2)
        {
            basis = Teuchos::rcp(new Basis_HGRAD_TET_C2_FEM<double, FieldContainer<double> >);
        }
    }
    else if(cell_type == "Hexahedron")
    {
        if(order == 1)
        {
            basis = Teuchos::rcp(new Basis_HGRAD_HEX_C1_FEM<double, FieldContainer<double> >);
        }
        else if(order == 2)
        {
            basis = Teuchos::rcp(new Basis_HGRAD_HEX_C2_FEM<double, FieldContainer<double> >);
        }
    }
    else
    {
        std::string ioError = descriptor + "Unknown cell type or invalid basis order!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
}

