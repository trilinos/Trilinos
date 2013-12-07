#include "mex.h"
#include "Intrepid_CellTools.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_getPhysicalSideNormals ..... MEX interface for the Intrepid (Trilinos) function Intrepid::CellTools::getPhysicalSideNormals.\n\n"
                    "\tintrepid_getPhysicalSideNormals(sideNormals,Jacobians,subcellOrd,parentCell)\n\n"
                    "\t<1-in/out> sideNormals = 3D array of size [spaceDim x #cubPoints x #cells] normals at cell sides\n"
                    "\t<2-in>     Jacobians = 4D array of size [spaceDim x spaceDim x #cubPoints x #cells] with Jacobians at reference side points\n"
                    "\t<3-in>     subcellOrd = subcell ordinal (int)\n"
                    "\t<4-in>     parentCellType = 'Triangle' | 'Quadrilateral' | 'Tetrahedron' | 'Hexahedron' (string)\n\n");

    // Check the number of input arguments
    if(nInput != 4)
    {
        std::string ioError = descriptor + "Incorrect number of input arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
    if(nOutput != 0)
    {
        std::string ioError = descriptor + "There can be no output arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }

    // Get the dimensions of the reference workset side normals array
    Teuchos::Array<int> ref_side_normals_dims;
    m2iGetArrayDims(ref_side_normals_dims, pInput[0]);
    // Get the dimensions of the workset jacobians array
    Teuchos::Array<int> jacobians_dims;
    m2iGetArrayDims(jacobians_dims, pInput[1]);

    // Get the (pointers to) data
    double* ref_side_normals_raw = mxGetPr(pInput[0]);
    double* jacobians_raw = mxGetPr(pInput[1]);
    const int subcellOrd = (int) mxGetScalar(pInput[2]);
    const std::string cell_type(mxArrayToString(pInput[3]));

    // Create cell topology
    Teuchos::RCP<shards::CellTopology> cellTopo;
    m2iGetCellTopo(cellTopo, cell_type, descriptor);

    FieldContainer<double> refSideNormals(ref_side_normals_dims, ref_side_normals_raw);
    FieldContainer<double> jacobians(jacobians_dims, jacobians_raw);

    try
    {
        CellTools<double>::getPhysicalSideNormals(refSideNormals, jacobians, subcellOrd, *(cellTopo.get()));
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_mapToReferenceSubcell\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
