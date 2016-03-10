// @HEADER
// ************************************************************************
//
//                           Intrepid2 Package
//                 Copyright (2007) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Kyungjoo Kim  (kyukim@sandia.gov), or
//                    Mauro Perego  (mperego@sandia.gov)
//
// ************************************************************************
// @HEADER

#include "mex.h"
#include <Intrepid2_FieldContainer.hpp>
#include "m2i_helpers.hpp"

using namespace Intrepid2;

void mexFunction(int nOutput, mxArray *pOutput[], int nInput, const mxArray *pInput[])
{

    std::string descriptor =
            ("\n evaluateVectorGradField ..... MEX interface for function evaluateVectorGradField.\n\n"
                    "\t evaluateVectorGradField(outputFields,inputData,inputFields)\n\n"
                    "\t<1-in/out> outputFields = Output fields array (4D array of size [#spaceDim x #spaceDim x #cubPoints x #cells])\n"
                    "\t<2-in> 	   inputData = Input data array (3D array of size [#spaceDim x #numFields x #cells])\n"
                    "\t<3-in> 	   inputFields = Input fields array (3D array of size [#cubPoints x #fields x #cells])\n"
                    "\t<4-in>     transpose = If 'T', use transposed tensor; if 'N', no transpose. Default: 'N'\n\n");

    // Check the number of input arguments
    if((nInput != 3) && (nInput != 4))
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

    char transpose = 'N';
    if(nInput == 3)
    {
        transpose = 'N';
    }
    else if(nInput == 4)
    {
        const std::string si_str(mxArrayToString(pInput[4]));
        if((si_str == "T") || (si_str == "t"))
        {
            transpose = 'T';
        }
        else
        {
            std::string ioError = descriptor + ("ACCEPTABLE INPUT ARGUMENTS: EITHER 'T' OR 't' = TRANSPOSE.\n"
                    "AND EITHER 'N' OR 'n' = NO TRANSPOSE.\n\n");
            mexErrMsgTxt(ioError.c_str());
        }
    }

    FieldContainer<double> outputFields(oFields_dims, oFields_raw);
    FieldContainer<double> inputData(iData_dims, iData_raw);
    FieldContainer<double> inputFields(iFields_dims, iFields_raw);

    // get sizes
    int dataRank = inputData.rank();
    int numDataPts = inputData.dimension(1);
    int inRank = inputFields.rank();
    int numCells = outputFields.dimension(0);
    int numFields = outputFields.dimension(1);
    int numPoints = outputFields.dimension(2);
    int matDim = outputFields.dimension(3);

    /*********************************************************************************************
     *                              inputFields is (C,F,P,D)                                     *
     *********************************************************************************************/
    if(inRank == 4)
    {
        if(numDataPts != 1)
        {
            // non-constant data
            switch(dataRank)
            {
                case 2:
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int field = 0; field < numFields; field++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputFields(cell, field, point, row) = inputData(cell, point)
                                            * inputFields(cell, field, point, row);
                                } // Row-loop
                            } // P-loop
                        } // F-loop
                    } // C-loop
                    break;
                case 3:
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int field = 0; field < numFields; field++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputFields(cell, field, point, row) = inputData(cell, point, row)
                                            * inputFields(cell, field, point, row);
                                } // Row-loop
                            } // P-loop
                        } // F-loop
                    } // C-loop
                    break;
                case 4:
                    if((transpose == 'n') || (transpose == 'N'))
                    {
                        for(int cell = 0; cell < numCells; cell++)
                        {
                            for(int field = 0; field < numFields; field++)
                            {
                                for(int point = 0; point < numPoints; point++)
                                {
                                    for(int row = 0; row < matDim; row++)
                                    {
                                        outputFields(cell, field, point, row) = 0.0;
                                        for(int col = 0; col < matDim; col++)
                                        {
                                            outputFields(cell, field, point, row) += inputData(cell, point, row, col)
                                                    * inputFields(cell, field, point, col);
                                        } // col
                                    } //row
                                } // point
                            } // field
                        } // cell
                    } // no transpose
                    else if((transpose == 't') || (transpose == 'T'))
                    {
                        for(int cell = 0; cell < numCells; cell++)
                        {
                            for(int field = 0; field < numFields; field++)
                            {
                                for(int point = 0; point < numPoints; point++)
                                {
                                    for(int row = 0; row < matDim; row++)
                                    {
                                        outputFields(cell, field, point, row) = 0.0;
                                        for(int col = 0; col < matDim; col++)
                                        {
                                            outputFields(cell, field, point, row) += inputData(cell, point, col, row)
                                                    * inputFields(cell, field, point, col);
                                        } // col
                                    } //row
                                } // point
                            } // field
                        } // cell
                    } //transpose
                    else
                    {
                        TEUCHOS_TEST_FOR_EXCEPTION(
                                !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ),
                                std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataField): The transpose flag must be 'n', 'N', 't' or 'T'.");
                    }
                    break;

                default:
                    TEUCHOS_TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ),
                            std::invalid_argument,
                            ">>> ERROR (ArrayTools::matvecProductDataField): inputData rank 2, 3 or 4 required.")
            } // switch inputData rank
        }
        else
        { // constant data case
            switch(dataRank)
            {
                case 2:
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int field = 0; field < numFields; field++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputFields(cell, field, point, row) = inputData(cell, 0)
                                            * inputFields(cell, field, point, row);
                                } // Row-loop
                            } // P-loop
                        } // F-loop
                    } // C-loop
                    break;

                case 3:
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int field = 0; field < numFields; field++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputFields(cell, field, point, row) = inputData(cell, 0, row)
                                            * inputFields(cell, field, point, row);
                                } // Row-loop
                            } // P-loop
                        } // F-loop
                    } // C-loop
                    break;

                case 4:
                    if((transpose == 'n') || (transpose == 'N'))
                    {
                        for(int cell = 0; cell < numCells; cell++)
                        {
                            for(int field = 0; field < numFields; field++)
                            {
                                for(int point = 0; point < numPoints; point++)
                                {
                                    for(int row = 0; row < matDim; row++)
                                    {
                                        outputFields(cell, field, point, row) = 0.0;
                                        for(int col = 0; col < matDim; col++)
                                        {
                                            outputFields(cell, field, point, row) += inputData(cell, 0, row, col)
                                                    * inputFields(cell, field, point, col);
                                        } // col
                                    } //row
                                } // point
                            } // field
                        } // cell
                    } // no transpose
                    else if((transpose == 't') || (transpose == 'T'))
                    {
                        for(int cell = 0; cell < numCells; cell++)
                        {
                            for(int field = 0; field < numFields; field++)
                            {
                                for(int point = 0; point < numPoints; point++)
                                {
                                    for(int row = 0; row < matDim; row++)
                                    {
                                        outputFields(cell, field, point, row) = 0.0;
                                        for(int col = 0; col < matDim; col++)
                                        {
                                            outputFields(cell, field, point, row) += inputData(cell, 0, col, row)
                                                    * inputFields(cell, field, point, col);
                                        } // col
                                    } //row
                                } // point
                            } // field
                        } // cell
                    } //transpose
                    else
                    {
                        TEUCHOS_TEST_FOR_EXCEPTION(
                                !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ),
                                std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataField): The transpose flag must be 'n', 'N', 't' or 'T'.");
                    }
                    break;

                default:
                    TEUCHOS_TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ),
                            std::invalid_argument,
                            ">>> ERROR (ArrayTools::matvecProductDataField): inputData rank 2, 3 or 4 required.")
            } // switch inputData rank
        } // end constant data case
    } // inputFields rank 4
    /*********************************************************************************************
     *                              inputFields is (F,P,D)                                       *
     *********************************************************************************************/
    else if(inRank == 3)
    {
        if(numDataPts != 1)
        { // non-constant data

            switch(dataRank)
            {
                case 2:
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int field = 0; field < numFields; field++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputFields(cell, field, point, row) = inputData(cell, point)
                                            * inputFields(field, point, row);
                                } // Row-loop
                            } // P-loop
                        } // F-loop
                    } // C-loop
                    break;

                case 3:
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int field = 0; field < numFields; field++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputFields(cell, field, point, row) = inputData(cell, point, row)
                                            * inputFields(field, point, row);
                                } // Row-loop
                            } // P-loop
                        } // F-loop
                    } // C-loop
                    break;

                case 4:
                    if((transpose == 'n') || (transpose == 'N'))
                    {
                        for(int cell = 0; cell < numCells; cell++)
                        {
                            for(int field = 0; field < numFields; field++)
                            {
                                for(int point = 0; point < numPoints; point++)
                                {
                                    for(int row = 0; row < matDim; row++)
                                    {
                                        outputFields(cell, field, point, row) = 0.0;
                                        for(int col = 0; col < matDim; col++)
                                        {
                                            outputFields(cell, field, point, row) += inputData(cell, point, row, col)
                                                    * inputFields(field, point, col);
                                        } // col
                                    } //row
                                } // point
                            } // field
                        } // cell
                    } // no transpose
                    else if((transpose == 't') || (transpose == 'T'))
                    {
                        for(int cell = 0; cell < numCells; cell++)
                        {
                            for(int field = 0; field < numFields; field++)
                            {
                                for(int point = 0; point < numPoints; point++)
                                {
                                    for(int row = 0; row < matDim; row++)
                                    {
                                        outputFields(cell, field, point, row) = 0.0;
                                        for(int col = 0; col < matDim; col++)
                                        {
                                            outputFields(cell, field, point, row) += inputData(cell, point, col, row)
                                                    * inputFields(field, point, col);
                                        } // col
                                    } //row
                                } // point
                            } // field
                        } // cell
                    } //transpose
                    else
                    {
                        TEUCHOS_TEST_FOR_EXCEPTION(
                                !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ),
                                std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataField): The transpose flag must be 'n', 'N', 't' or 'T'.");
                    }
                    break;

                default:
                    TEUCHOS_TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ),
                            std::invalid_argument,
                            ">>> ERROR (ArrayTools::matvecProductDataField): inputData rank 2, 3 or 4 required.")
            } // switch inputData rank
        }
        else
        { // constant data case
            switch(dataRank)
            {
                case 2:
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int field = 0; field < numFields; field++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputFields(cell, field, point, row) = inputData(cell, 0)
                                            * inputFields(field, point, row);
                                } // Row-loop
                            } // P-loop
                        } // F-loop
                    } // C-loop
                    break;

                case 3:
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int field = 0; field < numFields; field++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputFields(cell, field, point, row) = inputData(cell, 0, row)
                                            * inputFields(field, point, row);
                                } // Row-loop
                            } // P-loop
                        } // F-loop
                    } // C-loop
                    break;

                case 4:
                    if((transpose == 'n') || (transpose == 'N'))
                    {
                        for(int cell = 0; cell < numCells; cell++)
                        {
                            for(int field = 0; field < numFields; field++)
                            {
                                for(int point = 0; point < numPoints; point++)
                                {
                                    for(int row = 0; row < matDim; row++)
                                    {
                                        outputFields(cell, field, point, row) = 0.0;
                                        for(int col = 0; col < matDim; col++)
                                        {
                                            outputFields(cell, field, point, row) += inputData(cell, 0, row, col)
                                                    * inputFields(field, point, col);
                                        } // col
                                    } //row
                                } // point
                            } // field
                        } // cell
                    } // no transpose
                    else if((transpose == 't') || (transpose == 'T'))
                    {
                        for(int cell = 0; cell < numCells; cell++)
                        {
                            for(int field = 0; field < numFields; field++)
                            {
                                for(int point = 0; point < numPoints; point++)
                                {
                                    for(int row = 0; row < matDim; row++)
                                    {
                                        outputFields(cell, field, point, row) = 0.0;
                                        for(int col = 0; col < matDim; col++)
                                        {
                                            outputFields(cell, field, point, row) += inputData(cell, 0, col, row)
                                                    * inputFields(field, point, col);
                                        } // col
                                    } //row
                                } // point
                            } // field
                        } // cell
                    } //transpose
                    else
                    {
                        TEUCHOS_TEST_FOR_EXCEPTION(
                                !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ),
                                std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataField): The transpose flag must be 'n', 'N', 't' or 'T'.");
                    }
                    break;

                default:
                    TEUCHOS_TEST_FOR_EXCEPTION( !( (dataRank == 2) || (dataRank == 3) || (dataRank == 4) ),
                            std::invalid_argument,
                            ">>> ERROR (ArrayTools::matvecProductDataField): inputData rank 2, 3 or 4 required.")
            } // switch inputData rank
        } // end constant data case
    } // inputFields rank 3
    else
    {
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                ">>> ERROR (ArrayTools::matvecProductDataField): inputFields rank 3 or 4 required.")
    } // rank error
}
