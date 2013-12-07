#include "mex.h"
#include <Intrepid_FieldContainer.hpp>
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_multiplyVectorDataField ..... MEX interface for function intrepid_multiplyVectorDataField.\n\n"
                    "\tintrepid_multiplyVectorDataField(outputFields,inputData,inputFields)\n\n"
                    "\t<1-in/out> outputFields = Output fields array (4D array of size [#cubPoints x #fields x #cells x spaceDim])\n"
                    "\t<2-in> 	  inputData = Input data array (3D array of size [#cubPoints x #cells x spaceDim])\n"
                    "\t<3-in> 	  inputFields = Input fields array (5D array of size [#cubPoints x #fields x #cells x spaceDim x spaceDim] or 4D array of size [#cubPoints x #fields x #cells x spaceDim] or 3D array of size [#cubPoints x #fields x #cells])\n");

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

    // Get the dimensions of the output field values array
    Teuchos::Array<int> oFields_dims;
    m2iGetArrayDims(oFields_dims, pInput[0]);
    // Get the dimensions of the input vector data values array
    Teuchos::Array<int> iData_dims;
    m2iGetArrayDims(iData_dims, pInput[1]);
    // Get the dimensions of the input field values array
    Teuchos::Array<int> iFields_dims;
    m2iGetArrayDims(iFields_dims, pInput[2]);

    // Get the (pointers to) data
    double* oFields_raw = mxGetPr(pInput[0]);
    double* iData_raw = mxGetPr(pInput[1]);
    double* iFields_raw = mxGetPr(pInput[2]);

    FieldContainer<double> outputFields(oFields_dims, oFields_raw);
    FieldContainer<double> inputData(iData_dims, iData_raw);
    FieldContainer<double> inputFields(iFields_dims, iFields_raw);

    // get sizes
    int invalRank = inputFields.rank();
    int outvalRank = outputFields.rank();
    int numCells = outputFields.dimension(0);
    int numFields = outputFields.dimension(1);
    int numPoints = outputFields.dimension(2);
    int spaceDim = outputFields.dimension(3);

    try
    {
        for(int cl = 0; cl < numCells; cl++)
        {
            for(int bf = 0; bf < numFields; bf++)
            {
                for(int pt = 0; pt < numPoints; pt++)
                {
                    for(int iVec = 0; iVec < spaceDim; iVec++)
                    {
                        outputFields(cl, bf, pt, iVec) = inputFields(cl, bf, pt) * inputData(cl, pt, iVec);
                    } // D-loop
                } // P-loop
            } // F-loop
        } // C-loop

        TEUCHOS_TEST_FOR_EXCEPTION(!((invalRank == 3)), std::invalid_argument,
                ">>> ERROR: This method is defined only for rank-3 input containers.");

    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_scalarMultiplyDataField.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}
