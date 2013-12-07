#include "mex.h"
#include "Intrepid_FunctionSpaceTools.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_computeCellMeasure ..... MEX interface for the Intrepid (Trilinos) function Intrepid::FunctionSpaceTools::computeCellMeasure.\n\n"
                    "\tintrepid_computeCellMeasure(outMeasures,inDets,inWeights)\n\n"
                    "\t<1-in/out> outMeasures = Weighted measures (2D array of size [#cubPoints x #cells])\n"
                    "\t<2-in>     inDets = Determinants of cell Jacobians at cubature points (2D array of size [#cubPoints x #cells])\n"
                    "\t<3-in>     inWeights = Cubature weights (1D array of size [#cubPoints])\n\n");

    // Check the number of input arguments
    if(nInput != 3)
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
    // Get the dimensions of the determinant array
    Teuchos::Array<int> dets_dims;
    m2iGetArrayDims(dets_dims, pInput[1]);
    // Get the dimensions of the weights array
    Teuchos::Array<int> weights_dims;
    m2iGetArrayDims(weights_dims, pInput[2], true);

    // Get the (pointers to) data
    double* measures_raw = mxGetPr(pInput[0]);
    double* dets_raw = mxGetPr(pInput[1]);
    double* weights_raw = mxGetPr(pInput[2]);

    FieldContainer<double> cellMeasures(measures_dims, measures_raw);
    FieldContainer<double> cellDets(dets_dims, dets_raw);
    FieldContainer<double> cubWeights(weights_dims, weights_raw);

    try
    {
        FunctionSpaceTools::computeCellMeasure<double>(cellMeasures, cellDets, cubWeights);
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_computeCellMeasure.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
