#include "mex.h"
#include "Intrepid_DefaultCubatureFactory.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_getNumCubaturePoints ..... MEX interface for the Intrepid (Trilinos) function Intrepid::Cubature::getNumPoints.\n\n"
                    "\t[numPts] = intrepid_getNumCubaturePoints(cellType,degree)\n\n"
                    "\t<return> numPts = Number of cubature points (integer)\n"
                    "\t<1-in>   cellType = 'Line' | 'Triangle' | 'Quadrilateral' | 'Tetrahedron' | 'Hexahedron' (string)\n"
                    "\t<2-in>   degree = Degree of polynomials to be integrated exactly (integer)\n\n");

    // Check the number of input arguments
    if(nInput != 2)
    {
        std::string ioError = descriptor + "Incorrect number of input arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
    if(nOutput > 1)
    {
        std::string ioError = descriptor + "Incorrect number of output arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }

    // Get the (pointers to) data
    const std::string cell_type(mxArrayToString(pInput[0]));
    const int cubDegree = (int) mxGetScalar(pInput[1]);

    // Create cell topology
    Teuchos::RCP<shards::CellTopology> cellTopo;
    m2iGetCellTopo(cellTopo, cell_type, descriptor);

    // Create the cubature factory
    DefaultCubatureFactory<double> cubFactory;
    Teuchos::RCP<Cubature<double> > cellCub = cubFactory.create(*(cellTopo.get()), cubDegree);

    pOutput[0] = mxCreateNumericMatrix(1, 1, mxUINT32_CLASS, mxREAL);
    int* numPoints = (int*) mxGetPr(pOutput[0]);

    try
    {
        numPoints[0] = cellCub->getNumPoints();
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_getNumCubaturePoints.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
