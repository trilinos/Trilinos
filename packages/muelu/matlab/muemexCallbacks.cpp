#include "muemexCallbacks.h"

using namespace std;
using namespace Teuchos;

void MuemexCallback::callMatlabNoArgs(std::string function)
{
	int result = mexEvalString(function.c_str());
	if(result != 0)
		mexPrintf("An error occurred while running a MATLAB command.");\
}

vector<RCP<MuemexArg>> MuemexCallback::callMatlab(std::string function, int numOutputs, vector<RCP<MuemexArg>> args)
{
	mxArray** matlabArgs = new mxArray*[args.size()];
	mxArray** matlabOutput = new mxArray*[numOutputs];
	vector<RCP<MuemexArg>> output;
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
		catch (exception& e)
		{
			mexPrintf("An error occurred while converting arg #%d to MATLAB:\n", i);
			cout << e.what() << endl;
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
						throw runtime_error("Error: Don't know what to do with integer array.\n");
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
					throw runtime_error("MATLAB returned an unsupported type as a function output.\n");
			}
		}
		catch(exception& e)
		{
			mexPrintf("An error occurred while converting output #%d from MATLAB:\n", i);
			cout << e.what() << endl;
		}
	}
	return output;
}
