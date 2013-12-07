#include "mex.h"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_computeEdgeMeasure ..... MEX interface for the Intrepid (Trilinos) function Intrepid::FunctionSpaceTools::computeEdgeMeasure.\n\n"
                    "\tintrepid_computeEdgeMeasure(outMeasures,inJacobians,inWeights,whichEdge,parentCell)\n\n"
                    "\t<1-in/out> outMeasures = Weighted edge measures (2D array of size [#cubPoints x #cells])\n"
                    "\t<2-in>     inJacobians = Cell Jacobians at cubature points (4D array of size [#spaceDim x #spaceDim x #cubPoints x #cells])\n"
                    "\t<3-in>     inWeights = Cubature weights (1D array of size [#cubPoints])\n"
                    "\t<4-in>     whichEdge = Index of the edge subcell relative to the parent cell; defines the domain of integration (int)\n"
                    "\t<5-in>     parentCell = Parent cell topology (string)\n\n");

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

    // Get the dimensions of the measures array
    Teuchos::Array<int> measures_dims;
    m2iGetArrayDims(measures_dims, pInput[0]);
    // Get the dimensions of the jacobians array
    Teuchos::Array<int> jac_dims;
    m2iGetArrayDims(jac_dims, pInput[1]);
    // Get the dimensions of the weights array
    Teuchos::Array<int> weights_dims;
    m2iGetArrayDims(weights_dims, pInput[2], true);

    // Get the (pointers to) data
    double* measures_raw = mxGetPr(pInput[0]);
    double* jac_raw = mxGetPr(pInput[1]);
    double* weights_raw = mxGetPr(pInput[2]);
    const int whichEdge = (int) mxGetScalar(pInput[3]);
    const std::string cell_type(mxArrayToString(pInput[4]));

    // Create cell topology
    Teuchos::RCP<shards::CellTopology> cellTopo;
    m2iGetCellTopo(cellTopo, cell_type, descriptor);

    FieldContainer<double> cellMeasures(measures_dims, measures_raw);
    FieldContainer<double> cellJac(jac_dims, jac_raw);
    FieldContainer<double> cubWeights(weights_dims, weights_raw);

    try
    {
        FunctionSpaceTools::computeEdgeMeasure<double>(cellMeasures, cellJac, cubWeights, whichEdge, *(cellTopo.get()));
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_computeEdgeMeasure.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
