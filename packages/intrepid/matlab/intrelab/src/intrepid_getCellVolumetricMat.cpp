#include "mex.h"
#include <Intrepid_FieldContainer.hpp>
#include "m2i_helpers.hpp"

using namespace Intrepid;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\nintrepid_getCellVolumetricMat ..... MEX interface for function intrepid_getCellVolumetricMat.\n\n"
                    "\tintrepid_getCellVolumetricMat(outputFields)\n\n"
                    "\t<1-in/out> outputCellMat = Output fields array (3D array of size [#cells x #fields x #fields])\n"
                    "\t<2-in>     spaceDim = spatial dimensions (int)\n\n");

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

    // Get the dimensions of the output field values array
    Teuchos::Array<int> oFields_dims;
    m2iGetArrayDims(oFields_dims, pInput[0]);

    // Get the (pointers to) data
    double* oFields_raw = mxGetPr(pInput[0]);
    const int spaceDim = (int) mxGetScalar(pInput[1]);

    FieldContainer<double> outputFields(oFields_dims, oFields_raw);

    // get sizes
    int outvalRank = outputFields.rank();
    int numCells = outputFields.dimension(0);
    int numPoints = outputFields.dimension(1);

    try
    {

        for(int cell = 0; cell < numCells; cell++)
        {
            for(int pt = 0; pt < numPoints; pt++)
            {
                switch(spaceDim)
                {
                    case 2:
                    {
                        outputFields(cell, pt, 0, 0) = 1.0;
                        outputFields(cell, pt, 0, 1) = 1.0;
                        outputFields(cell, pt, 1, 0) = 1.0;
                        outputFields(cell, pt, 1, 1) = 1.0;
                        break;
                    }
                    case 3:
                    {
                        outputFields(cell, pt, 0, 0) = 1.0;
                        outputFields(cell, pt, 0, 1) = 1.0;
                        outputFields(cell, pt, 0, 2) = 1.0;
                        outputFields(cell, pt, 1, 0) = 1.0;
                        outputFields(cell, pt, 1, 1) = 1.0;
                        outputFields(cell, pt, 1, 2) = 1.0;
                        outputFields(cell, pt, 2, 0) = 1.0;
                        outputFields(cell, pt, 2, 1) = 1.0;
                        outputFields(cell, pt, 2, 2) = 1.0;
                        break;
                    }
                }
            } // P-loop
        } // C-loop
        TEUCHOS_TEST_FOR_EXCEPTION(!((outvalRank == 4)), std::invalid_argument,
                ">>> ERROR: This method is defined only for rank-4 output containers.");

    } catch(const std::exception &e)
    {
        std::string intrepiderr = e.what();
        std::string matlaberr = "------------------------------------------------------------\n"
                + ("MATLAB returned:  Invalid arguments in the call to intrepid_getCellVolumetricMat.\n"
                        + ("Intrepid (Trilinos) returned:\n" + intrepiderr));
        mexErrMsgTxt(matlaberr.c_str());
    }
}

