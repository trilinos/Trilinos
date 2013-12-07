#include "mex.h"
#include "Intrepid_CellTools.hpp"
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_setJacobianInv ..... MEX interface for the Intrepid (Trilinos) function Intrepid::CellTools::setJacobianInv.\n\n"
                    "\tintrepid_setJacobian(jacobianInv,jacobian)\n\n"
                    "\t<1-in/out> jacobianInv = Inverse Jacobians (4D array of size [spaceDim x spaceDim x #evalPoints x #cells] or 3D array of size [spaceDim x spaceDim x #evalPoints])\n"
                    "\t<2-in>     jacobian = Jacobians of reference-to-physical map (4D array of size [spaceDim x spaceDim x #evalPoints x #cells] or 3D array of size [spaceDim x spaceDim x #evalPoints])\n\n");

    // Check the number of input arguments
    if(nInput != 2)
    {
        std::string ioError = descriptor + "Incorrect number of input arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }
    if(nOutput != 0)
    {
        std::string ioError = descriptor + "There can be no output arguments!!!\n";
        mexErrMsgTxt(ioError.c_str());
    }

    // Get the dimensions of the Jacobian arrays
    Teuchos::Array<int> jacs_dims;
    m2iGetArrayDims(jacs_dims, pInput[0]);

    // Get the (pointers to) data
    double* jacinv_raw = mxGetPr(pInput[0]);
    double* jacs_raw = mxGetPr(pInput[1]);

    FieldContainer<double> cellJacInv(jacs_dims, jacinv_raw);
    FieldContainer<double> cellJacs(jacs_dims, jacs_raw);

    try
    {
        CellTools<double>::setJacobianInv(cellJacInv, cellJacs);
    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_setJacobianInv.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
