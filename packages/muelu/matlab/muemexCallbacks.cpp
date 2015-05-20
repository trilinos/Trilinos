// @HEADER
//
// ***********************************************************************
//
//                MueLu: A package for multigrid based preconditioning
//                                      Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact
//                                        Jonathan Hu           (jhu@sandia.gov)
//                                        Andrey Prokopenko (aprokop@sandia.gov)
//                                        Ray Tuminaro          (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "muemexCallbacks.h"
#include "muemexTypes_def.hpp"
#include "MueLu_Utilities.hpp"

using namespace Teuchos;

void MuemexCallback::callMatlabNoArgs(std::string function)
{
  int result = mexEvalString(function.c_str());
  if(result != 0)
    mexPrintf("An error occurred while running a MATLAB command.");\
}

std::vector<RCP<MuemexArg>> MuemexCallback::callMatlab(std::string function, int numOutputs, std::vector<RCP<MuemexArg>> args)
{
  mxArray** matlabArgs = new mxArray*[args.size()];
  mxArray** matlabOutput = new mxArray*[numOutputs];
  std::vector<RCP<MuemexArg>> output;
  for(int i = 0; i < int(args.size()); i++)
    {
      try
        {
          switch(args[i]->type)
            {
            case INT:
              matlabArgs[i] = rcp_static_cast<MuemexData<int>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case DOUBLE:
              matlabArgs[i] = rcp_static_cast<MuemexData<double>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case STRING:
              matlabArgs[i] = rcp_static_cast<MuemexData<string>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case COMPLEX:
              matlabArgs[i] = rcp_static_cast<MuemexData<complex_t>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case XPETRA_ORDINAL_VECTOR:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Xpetra_ordinal_vector>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case TPETRA_MULTIVECTOR_DOUBLE:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case TPETRA_MULTIVECTOR_COMPLEX:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case TPETRA_MATRIX_DOUBLE:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Tpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case TPETRA_MATRIX_COMPLEX:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case XPETRA_MATRIX_DOUBLE:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Xpetra_Matrix_double>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case XPETRA_MATRIX_COMPLEX:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Xpetra_Matrix_complex>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case XPETRA_MULTIVECTOR_DOUBLE:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Xpetra_MultiVector_double>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case XPETRA_MULTIVECTOR_COMPLEX:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Xpetra_MultiVector_complex>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case EPETRA_CRSMATRIX:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Epetra_CrsMatrix>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            case EPETRA_MULTIVECTOR:
              matlabArgs[i] = rcp_static_cast<MuemexData<RCP<Epetra_MultiVector>>, MuemexArg>(args[i])->convertToMatlab();
              break;
            }
        }
      catch (std::exception& e)
        {
          mexPrintf("An error occurred while converting arg #%d to MATLAB:\n", i);
	  std::cout << e.what() << std::endl;
          mexPrintf("Passing 0 instead.\n");
          matlabArgs[i] = mxCreateDoubleScalar(0);
        }
    }
  //now matlabArgs is populated with MATLAB data types
  int result = mexCallMATLAB(numOutputs, matlabOutput, args.size(), matlabArgs, function.c_str());
  if(result != 0)
    mexPrintf("Matlab encountered an error while running command through muemexCallbacks.\n");
  //now, if all went well, matlabOutput contains all the output to return to user
  for(int i = 0; i < numOutputs; i++)
    {
      try
        {
          //Identify the type of each output, and put into output vector
          mxArray* item = matlabOutput[i];
          switch(mxGetClassID(item))
            {
            case mxCHAR_CLASS:
              //string
              output.push_back(rcp(new MuemexData<string>(item)));
              break;
            case mxINT32_CLASS:
              if(mxGetM(item) == 1 && mxGetN(item) == 1)
                //single int
                output.push_back(rcp(new MuemexData<int>(item)));
              else if(mxGetM(item) != 1 || mxGetN(item) != 1)
                //ordinal vector
                output.push_back(rcp(new MuemexData<RCP<Xpetra_ordinal_vector>>(item)));
              else
                throw std::runtime_error("Error: Don't know what to do with integer array.\n");
              break;
            case mxDOUBLE_CLASS:
              if(mxGetM(item) == 1 && mxGetN(item) == 1)
                {
                  if(mxIsComplex(item))
                    //single double (scalar, real)
                    output.push_back(rcp(new MuemexData<double>(item)));
                  else
                    //single complex scalar
                    output.push_back(rcp(new MuemexData<complex_t>(item)));
                }
              else if(mxIsSparse(item))
                {
                  //Default to Tpetra matrix for this
                  if(mxIsComplex(item))
                    //complex Tpetra matrix (sparse)
                    output.push_back(rcp(new MuemexData<RCP<Tpetra_CrsMatrix_double>>(item)));
                  else
                    //real Tpetra matrix
                    output.push_back(rcp(new MuemexData<RCP<Tpetra_CrsMatrix_complex>>(item)));
                }
              else
                {
                  //Default to Tpetra multivector for this case
                  if(mxIsComplex(item))
                    output.push_back(rcp(new MuemexData<RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>(item)));
                  else
                    output.push_back(rcp(new MuemexData<RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>>>(item)));
                }
              break;
            default:
	        throw std::runtime_error("MATLAB returned an unsupported type as a function output.\n");
            }
        }
	      catch(std::exception& e)
        {
          mexPrintf("An error occurred while converting output #%d from MATLAB:\n", i);
	  std::cout << e.what() << std::endl;
        }
    }
  return output;
}
