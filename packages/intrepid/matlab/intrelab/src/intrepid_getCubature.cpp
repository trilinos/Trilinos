#include "mex.h"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_getCubature ..... MEX interface for the Intrepid (Trilinos) function Intrepid::Cubature::getCubature.\n\n"
                    "\tintrepid_getCubature(cubPoints,cubWeights,cellType,degree)\n\n"
                    "\t<1-in/out> cubPoints = Cubature points (2D array of size [spaceDim x #cubPoints])\n"
                    "\t<2-in/out> cubWeights = Cubature weights (1D array of size [#cubPoints])\n"
                    "\t<3-in>     cellType = 'Line' | 'Triangle' | 'Quadrilateral' | 'Tetrahedron' | 'Hexahedron' (string)\n"
                    "\t<4-in>     degree = Degree of polynomials to be integrated exactly (integer)\n\n");

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

    // Get the dimensions of the cubature point array
    Teuchos::Array<int> points_dims;
    m2iGetArrayDims(points_dims, pInput[0]);
    // Get the dimensions of the weights array
    Teuchos::Array<int> weights_dims;
    m2iGetArrayDims(weights_dims, pInput[1], true);

    // Get the (pointers to) data
    const std::string cell_type(mxArrayToString(pInput[2]));
    const int cubDegree = (int) mxGetScalar(pInput[3]);
    double* cub_points_raw = mxGetPr(pInput[0]);
    double* cub_weights_raw = mxGetPr(pInput[1]);

    // Create cell topology
    Teuchos::RCP<shards::CellTopology> cellTopo;
    m2iGetCellTopo(cellTopo, cell_type, descriptor);

    // Create the cubature factory
    DefaultCubatureFactory<double> cubFactory;
    Teuchos::RCP<Cubature<double> > cellCub = cubFactory.create(*(cellTopo.get()), cubDegree);

    FieldContainer<double> cub_points(points_dims, cub_points_raw);
    FieldContainer<double> cub_weights(weights_dims, cub_weights_raw);

    try
    {
        cellCub->getCubature(cub_points, cub_weights);
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_getCubature.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
