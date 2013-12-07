#include "mex.h"
#include "Intrepid_CellTools.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_mapToPhysicalFrame ..... MEX interface for the Intrepid (Trilinos) function Intrepid::CellTools::mapToPhysicalFrame.\n\n"
                    "\tintrepid_mapToPhysicalFrame(physPoints,refPoints,cellWorkset,cellType,whichCell=-1)\n\n"
                    "\t<1-in/out> physPoints = Images of the reference points (3D array of size [spaceDim x #evalPoints x #cells] or 2D array of size [spaceDim x #evalPoints])\n"
                    "\t<2-in>     refPoints = Points in reference frame (3D array of size [spaceDim x #evalPoints x #cells] or 2D array of size [spaceDim x #evalPoints])\n"
                    "\t<3-in>     cellWorkset = Cell nodes (3D array of size [spaceDim x #nodes x #cells])\n"
                    "\t<4-in>     cellType = 'Line' | 'Triangle' | 'Quadrilateral' | 'Tetrahedron' | 'Hexahedron' (string)\n"
                    "\t<5-in/opt> whichCell = Cell number (default: -1 means ALL cells) (int)\n\n");

    // Check the number of input arguments
    if((nInput != 4) && (nInput != 5))
    {
        std::string ioError = descriptor + "Incorrect number of input arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
    if(nOutput != 0)
    {
        std::string ioError = descriptor + "There can be no output arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }

    // Get the dimensions of the physical point array
    Teuchos::Array<int> phys_pts_dims;
    m2iGetArrayDims(phys_pts_dims, pInput[0]);
    // Get the dimensions of the reference point array
    Teuchos::Array<int> ref_points_dims;
    m2iGetArrayDims(ref_points_dims, pInput[1]);
    // Get the dimensions of the cell workset
    Teuchos::Array<int> workset_dims;
    m2iGetArrayDims(workset_dims, pInput[2], true);

    // Get the (pointers to) data
    double* phys_pts_raw = mxGetPr(pInput[0]);
    double* ref_points_raw = mxGetPr(pInput[1]);
    double* workset_raw = mxGetPr(pInput[2]);
    const std::string cell_type(mxArrayToString(pInput[3]));
    int whichCell = -1;
    if(nInput == 5)
    {
        whichCell = (int) mxGetScalar(pInput[4]);
    }

    // Create cell topology
    Teuchos::RCP<shards::CellTopology> cellTopo;
    m2iGetCellTopo(cellTopo, cell_type, descriptor);

    FieldContainer<double> physPoints(phys_pts_dims, phys_pts_raw);
    FieldContainer<double> refPoints(ref_points_dims, ref_points_raw);
    FieldContainer<double> cellWorkset(workset_dims, workset_raw);

    try
    {
        CellTools<double>::mapToPhysicalFrame(physPoints, refPoints, cellWorkset, *(cellTopo.get()), whichCell);
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_mapToPhysicalFrame.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
