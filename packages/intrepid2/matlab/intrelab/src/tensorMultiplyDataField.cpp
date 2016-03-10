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
            ("\n tensorMultiplyDataField ..... MEX interface for function tensorMultiplyDataField.\n\n"
                    "\t tensorMultiplyDataField(outputFields,inputData,inputFields)\n\n"
                    "\t<1-in/out> outputFields = Output fields array (4D array of size [matrixDim x #dof x #cubPoints x #cells])\n"
                    "\t<2-in>     inputData = Input data array (3D array size [#spaceDim x #cubPoints x # cells] or 4D array size [#spaceDim x #spaceDim x #cubPoints x # cells])\n"
                    "\t<3-in> 	   inputFields = Input fields array (4D array of size [#matrixDim x #dof x cubPoints x #cells] or 5D array of size [#matrixDim x #matrixDim x #dof x cubPoints x #cells])\n\n");

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
    // Get the dimensions of the input data values array
    Teuchos::Array<int> iData_dims;
    m2iGetArrayDims(iData_dims, pInput[1]);
    // Get the dimensions of the input field values array
    Teuchos::Array<int> iFields_dims;
    m2iGetArrayDims(iFields_dims, pInput[2]);

    // Get the (pointers to) data
    double* oFields_raw = mxGetPr(pInput[0]);
    double* iData_raw = mxGetPr(pInput[1]);
    double* iFields_raw = mxGetPr(pInput[2]);

    bool transpose = false;
    if(nInput == 3)
    {
        transpose = false;
    }
    else if(nInput == 4)
    {
        const std::string tr_str(mxArrayToString(pInput[3]));
        if(tr_str == "true")
        {
            transpose = true;
        }
    }

    FieldContainer<double> outputFields(oFields_dims, oFields_raw);
    FieldContainer<double> inputData(iData_dims, iData_raw);
    FieldContainer<double> inputFields(iFields_dims, iFields_raw);

    // get sizes
    int dataRank = inputData.rank();
    int fieldRank = inputFields.rank();
    int numCells = inputFields.dimension(0);
    int numDofs = inputFields.dimension(1);
    int numPoints = inputFields.dimension(2);
    int matDim = inputFields.dimension(3);

    /***********************************************************************
     *                    inputFields is (C,F,P,D)                         *
     ************************************************************************/
    if(fieldRank == 4)
    {
        switch(dataRank)
        {
            case 4:
                if(!transpose)
                {
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int dof = 0; dof < numDofs; dof++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputFields(cell, dof, point, row) = 0.0;
                                    for(int col = 0; col < matDim; col++)
                                    {
                                        outputFields(cell, dof, point, row) += inputData(cell, point, row, col)
                                                * inputFields(cell, dof, point, col);
                                    } // col
                                } //row
                            } // point
                        } // dof
                    } // cell
                } // no transpose
                else
                {
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int dof = 0; dof < numDofs; dof++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputFields(cell, dof, point, row) = 0.0;
                                    for(int col = 0; col < matDim; col++)
                                    {
                                        outputFields(cell, dof, point, row) += inputData(cell, point, col, row)
                                                * inputFields(cell, dof, point, col);
                                    } // col
                                } //row
                            } // point
                        } // dof
                    } // cell
                } // transpose
                break;

        } // switch data rank
    } // if statement

    /***********************************************************************
     *                    inputFields is (C,F,P,D,D)                         *
     ************************************************************************/

    if(fieldRank == 5)
    {
        switch(dataRank)
        {
            case 3:
                if(!transpose)
                {
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int dof = 0; dof < numDofs; dof++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputFields(cell, dof, point, row) = 0.0;
                                    for(int col = 0; col < matDim; col++)
                                    {
                                        outputFields(cell, dof, point, row) += inputFields(cell, dof, point, row, col)
                                                * inputData(cell, point, col);
                                    } // col
                                } //row
                            } // point
                        } // dof
                    } // cell
                } // no transpose
                else
                {
                    for(int cell = 0; cell < numCells; cell++)
                    {
                        for(int dof = 0; dof < numDofs; dof++)
                        {
                            for(int point = 0; point < numPoints; point++)
                            {
                                for(int row = 0; row < matDim; row++)
                                {
                                    outputFields(cell, dof, point, row) = 0.0;
                                    for(int col = 0; col < matDim; col++)
                                    {
                                        outputFields(cell, dof, point, row) += inputFields(cell, dof, point, col, row)
                                                * inputData(cell, point, col);
                                    } // col
                                } //row
                            } // point
                        } // dof
                    } // cell
                } // transpose
                break;

        } // switch data rank
    } // if statement
}
