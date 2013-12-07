#include "mex.h"
#include "Intrepid_CellTools.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_mapToReferenceSubcell ..... MEX interface for the Intrepid (Trilinos) function Intrepid::CellTools::mapToReferenceSubcell.\n\n"
                    "\tintrepid_mapToReferenceSubCell(refSubcellPoints,paramPoints,subcellDim,subcellOrd,parentCell)\n\n"
                    "\t<1-in/out> refSubcellPoints = 2D array of size [spaceDim x #evalPoints] with image of parameter space points\n"
                    "\t<2-in>     paramPoints = 2D array of size [spaceDim x #evalPoints] with points in 1D or 2D parameter domain\n"
                    "\t<3-in>     subcellDim = dimension of subcell where points are mapped to (int)\n"
                    "\t<4-in>     subcellOrd = subcell ordinal (int)\n"
                    "\t<5-in>     parentCellType = 'Triangle' | 'Quadrilateral' | 'Tetrahedron' | 'Hexahedron' (string)\n\n");

    // Check the number of input arguments
    if(nInput != 5)
    {
        std::string ioError = descriptor + "Incorrect number of input arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
    if(nOutput != 0)
    {
        std::string ioError = descriptor + "There can be no output arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }

    // Get the dimensions of the reference subcell point array
    Teuchos::Array<int> ref_subcell_points_dims;
    m2iGetArrayDims(ref_subcell_points_dims, pInput[0]);
    // Get the dimensions of the parameter point array
    Teuchos::Array<int> param_pts_dims;
    m2iGetArrayDims(param_pts_dims, pInput[1]);

    // Get the (pointers to) data
    double* ref_subcell_points_raw = mxGetPr(pInput[0]);
    double* param_pts_raw = mxGetPr(pInput[1]);
    const int subcellDim = (int) mxGetScalar(pInput[2]);
    const int subcellOrd = (int) mxGetScalar(pInput[3]);
    const std::string cell_type(mxArrayToString(pInput[4]));

    // Create cell topology
    Teuchos::RCP<shards::CellTopology> cellTopo;
    m2iGetCellTopo(cellTopo, cell_type, descriptor);

    FieldContainer<double> refSubcellPoints(ref_subcell_points_dims, ref_subcell_points_raw);
    FieldContainer<double> paramPoints(param_pts_dims, param_pts_raw);

    try
    {
        CellTools<double>::mapToReferenceSubcell(refSubcellPoints, paramPoints, subcellDim, subcellOrd,
                *(cellTopo.get()));
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_mapToReferenceSubcell\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
