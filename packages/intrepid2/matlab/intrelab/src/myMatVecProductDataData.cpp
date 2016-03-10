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
            ("\n myMatVecProductDataData ..... MEX interface for function myMatVecProductDataData.\n\n"
                    "\t evaluateVectorGradField(outputFields,inputData,inputFields)\n\n"
                    "\t<1-in/out> outputFields = Output fields array (4D array of size [#spaceDim x #spaceDim x #cubPoints x #cells])\n"
                    "\t<2-in>      inputData = Input data array (3D array of size [#spaceDim x #numFields x #cells])\n"
                    "\t<3-in>      inputFields = Input fields array (3D array of size [#cubPoints x #fields x #cells])\n"
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
    Teuchos::Array<int> oData_dims;
    m2iGetArrayDims(oData_dims, pInput[0]);
    // Get the dimensions of the input vector data values array
    Teuchos::Array<int> iDataLeft_dims;
    m2iGetArrayDims(iDataLeft_dims, pInput[1]);
    // Get the dimensions of the input field values array
    Teuchos::Array<int> iDataRight_dims;
    m2iGetArrayDims(iDataRight_dims, pInput[2]);

    // Get the (pointers to) data
    double* oData_raw = mxGetPr(pInput[0]);
    double* iDataLeft_raw = mxGetPr(pInput[1]);
    double* iDataRight_raw = mxGetPr(pInput[2]);

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

    FieldContainer<double> outputData(oData_dims, oData_raw);
    FieldContainer<double> inputDataLeft(iDataLeft_dims, iDataLeft_raw);
    FieldContainer<double> inputDataRight(iDataRight_dims, iDataRight_raw);

    // get sizes
    int dataLeftRank = inputDataLeft.rank();
    int numDataLeftPts = inputDataLeft.dimension(1);
    int dataRightRank = inputDataRight.rank();
    int numCells = outputData.dimension(0);
    int numPoints = outputData.dimension(1);
    int matDim = outputData.dimension(2);

    /*********************************************************************************************
     *                              inputDataRight is (C,P,D)                                   *
     *********************************************************************************************/
    if(dataRightRank == 3)
    {
        if(numDataLeftPts != 1)
        { // non-constant left data
            switch(dataLeftRank)
            {
                case 2:
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int point = 0; point < numPoints; point++)
                        {
                            for(int row = 0; row < matDim; row++)
                            {
                                outputData(cell, point, row) = inputDataLeft(cell, point)
                                        * inputDataRight(cell, point, row);
                            } // Row-loop
                        } // P-loop
                    } // C-loop
                    break;
                case 3:
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int point = 0; point < numPoints; point++)
                        {
                            for(int row = 0; row < matDim; row++)
                            {
                                outputData(cell, point, row) = inputDataLeft(cell, point, row)
                                        * inputDataRight(cell, point, row);
                            } // Row-loop
                        } // P-loop
                    } // C-loop
                    break;
                case 4:
                    if((transpose == 'n') || (transpose == 'N'))
                    {
                        for(int cell = 0; cell < numCells; cell++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputData(cell, point, row) = 0.0;
                                    for(int col = 0; col < matDim; col++)
                                    {
                                        outputData(cell, point, row) += inputDataLeft(cell, point, row, col)
                                                * inputDataRight(cell, point, col);
                                    } // col
                                } //row
                            } // point
                        } // cell
                    } // no transpose
                    else if((transpose == 't') || (transpose == 'T'))
                    {
                        for(int cell = 0; cell < numCells; cell++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputData(cell, point, row) = 0.0;
                                    for(int col = 0; col < matDim; col++)
                                    {
                                        outputData(cell, point, row) += inputDataLeft(cell, point, col, row)
                                                * inputDataRight(cell, point, col);
                                    } // col
                                } //row
                            } // point
                        } // cell
                    } //transpose
                    else
                    {
                        TEUCHOS_TEST_FOR_EXCEPTION(
                                !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ),
                                std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataData): The transpose flag must be 'n', 'N', 't' or 'T'.");
                    }
                    break;
                default:
                    TEUCHOS_TEST_FOR_EXCEPTION( !( (dataLeftRank == 2) || (dataLeftRank == 3) || (dataLeftRank == 4) ),
                            std::invalid_argument,
                            ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft rank 2, 3 or 4 required.")
            } // switch inputDataLeft rank
        }
        else
        { // constant data case
            switch(dataLeftRank)
            {
                case 2:
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int point = 0; point < numPoints; point++)
                        {
                            for(int row = 0; row < matDim; row++)
                            {
                                outputData(cell, point, row) = inputDataLeft(cell, 0)
                                        * inputDataRight(cell, point, row);
                            } // Row-loop
                        } // F-loop
                    } // C-loop
                    break;
                case 3:
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int point = 0; point < numPoints; point++)
                        {
                            for(int row = 0; row < matDim; row++)
                            {
                                outputData(cell, point, row) = inputDataLeft(cell, 0, row)
                                        * inputDataRight(cell, point, row);
                            } // Row-loop
                        } // P-loop
                    } // C-loop
                    break;
                case 4:
                    if((transpose == 'n') || (transpose == 'N'))
                    {
                        for(int cell = 0; cell < numCells; cell++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputData(cell, point, row) = 0.0;
                                    for(int col = 0; col < matDim; col++)
                                    {
                                        outputData(cell, point, row) += inputDataLeft(cell, 0, row, col)
                                                * inputDataRight(cell, point, col);
                                    } // col
                                } //row
                            } // point
                        } // cell
                    } // no transpose
                    else if((transpose == 't') || (transpose == 'T'))
                    {
                        for(int cell = 0; cell < numCells; cell++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputData(cell, point, row) = 0.0;
                                    for(int col = 0; col < matDim; col++)
                                    {
                                        outputData(cell, point, row) += inputDataLeft(cell, 0, col, row)
                                                * inputDataRight(cell, point, col);
                                    } // col
                                } //row
                            } // point
                        } // cell
                    } //transpose
                    else
                    {
                        TEUCHOS_TEST_FOR_EXCEPTION(
                                !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ),
                                std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataData): The transpose flag must be 'n', 'N', 't' or 'T'.");
                    }
                    break;
                default:
                    TEUCHOS_TEST_FOR_EXCEPTION( !( (dataLeftRank == 2) || (dataLeftRank == 3) || (dataLeftRank == 4) ),
                            std::invalid_argument,
                            ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft rank 2, 3 or 4 required.")
            } // switch inputDataLeft rank
        } // end constant data case
    } // inputDataRight rank 4
    /*********************************************************************************************
     *                              inputDataRight is (P,D)                                     *
     *********************************************************************************************/
    else if(dataRightRank == 2)
    {
        if(numDataLeftPts != 1)
        { // non-constant data
            switch(dataLeftRank)
            {
                case 2:
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int point = 0; point < numPoints; point++)
                        {
                            for(int row = 0; row < matDim; row++)
                            {
                                outputData(cell, point, row) = inputDataLeft(cell, point) * inputDataRight(point, row);
                            } // Row-loop
                        } // P-loop
                    } // C-loop
                    break;
                case 3:
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int point = 0; point < numPoints; point++)
                        {
                            for(int row = 0; row < matDim; row++)
                            {
                                outputData(cell, point, row) = inputDataLeft(cell, point, row)
                                        * inputDataRight(point, row);
                            } // Row-loop
                        } // P-loop
                    } // C-loop
                    break;
                case 4:
                    if((transpose == 'n') || (transpose == 'N'))
                    {
                        for(int cell = 0; cell < numCells; cell++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputData(cell, point, row) = 0.0;
                                    for(int col = 0; col < matDim; col++)
                                    {
                                        outputData(cell, point, row) += inputDataLeft(cell, point, row, col)
                                                * inputDataRight(point, col);
                                    } // col
                                } //row
                            } // point
                        } // cell
                    } // no transpose
                    else if((transpose == 't') || (transpose == 'T'))
                    {
                        for(int cell = 0; cell < numCells; cell++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputData(cell, point, row) = 0.0;
                                    for(int col = 0; col < matDim; col++)
                                    {
                                        outputData(cell, point, row) += inputDataLeft(cell, point, col, row)
                                                * inputDataRight(point, col);
                                    } // col
                                } //row
                            } // point
                        } // cell
                    } //transpose
                    else
                    {
                        TEUCHOS_TEST_FOR_EXCEPTION(
                                !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ),
                                std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataData): The transpose flag must be 'n', 'N', 't' or 'T'.");
                    }
                    break;
                default:
                    TEUCHOS_TEST_FOR_EXCEPTION( !( (dataLeftRank == 2) || (dataLeftRank == 3) || (dataLeftRank == 4) ),
                            std::invalid_argument,
                            ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft rank 2, 3 or 4 required.")
            } // switch inputDataLeft rank
        }
        else
        { // constant data case
            switch(dataLeftRank)
            {
                case 2:
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int point = 0; point < numPoints; point++)
                        {
                            for(int row = 0; row < matDim; row++)
                            {
                                outputData(cell, point, row) = inputDataLeft(cell, 0) * inputDataRight(point, row);
                            } // Row-loop
                        } // P-loop
                    } // C-loop
                    break;
                case 3:
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int point = 0; point < numPoints; point++)
                        {
                            for(int row = 0; row < matDim; row++)
                            {
                                outputData(cell, point, row) = inputDataLeft(cell, 0, row) * inputDataRight(point, row);
                            } // Row-loop
                        } // P-loop
                    } // C-loop
                    break;
                case 4:
                    if((transpose == 'n') || (transpose == 'N'))
                    {
                        for(int cell = 0; cell < numCells; cell++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputData(cell, point, row) = 0.0;
                                    for(int col = 0; col < matDim; col++)
                                    {
                                        outputData(cell, point, row) += inputDataLeft(cell, 0, row, col)
                                                * inputDataRight(point, col);
                                    } // col
                                } //row
                            } // point
                        } // cell
                    } // no transpose
                    else if((transpose == 't') || (transpose == 'T'))
                    {
                        for(int cell = 0; cell < numCells; cell++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputData(cell, point, row) = 0.0;
                                    for(int col = 0; col < matDim; col++)
                                    {
                                        outputData(cell, point, row) += inputDataLeft(cell, 0, col, row)
                                                * inputDataRight(point, col);
                                    } // col
                                } //row
                            } // point
                        } // cell
                    } //transpose
                    else
                    {
                        TEUCHOS_TEST_FOR_EXCEPTION(
                                !( (transpose == 'n') || (transpose == 'N') || (transpose == 't') || (transpose == 'T') ),
                                std::invalid_argument,
                                ">>> ERROR (ArrayTools::matvecProductDataData): The transpose flag must be 'n', 'N', 't' or 'T'.");
                    }
                    break;
                default:
                    TEUCHOS_TEST_FOR_EXCEPTION( !( (dataLeftRank == 2) || (dataLeftRank == 3) || (dataLeftRank == 4) ),
                            std::invalid_argument,
                            ">>> ERROR (ArrayTools::matvecProductDataData): inputDataLeft rank 2, 3 or 4 required.")
            } // switch inputDataLeft rank
        } // end constant inputDataLeft case
    } // inputDataRight rank 2
    else
    {
        TEUCHOS_TEST_FOR_EXCEPTION( true, std::invalid_argument,
                ">>> ERROR (ArrayTools::matvecProductDataData): inputDataRight rank 2 or 3 required.")
    } // rank error
}
