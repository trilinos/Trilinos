// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

#ifndef MUEMEX_H
#define MUEMEX_H

#include <cstdio>
#include <cstring>
#include <string>
#include <vector>
#include <complex>
#include <stdexcept>

#include "Teuchos_ParameterList.hpp"
#include "MueLu_config.hpp"
#include "MueLu.hpp"
#include "MueLu_EpetraOperator.hpp"
#include "MueLu_TpetraOperator.hpp"
#include "MueLu_Hierarchy_decl.hpp"
#include "MueLu_CreateEpetraPreconditioner.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Tpetra_CrsMatrix_decl.hpp"
#include "Xpetra_EpetraCrsMatrix.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosTpetraAdapter.hpp"
#include "BelosMueLuAdapter.hpp"
#include "mex.h"

//Chris, put the complex scalar support def here (I don't know what it's called)
//I will test HAVE_COMPLEX_SCALARS in my code to control what gets defined/allowed

//For now, always allow them
#define HAVE_TRILINOS_COMPLEX_SUPPORT

#ifdef HAVE_TRILINOS_COMPLEX_SUPPORT
#define HAVE_COMPLEX_SCALARS
#endif

typedef enum
{
    EPETRA,
    TPETRA,
    TPETRA_COMPLEX
} DataPackType;

//Mode value passed to MATLAB muemex function as 1st arg
typedef enum
{
    MODE_SETUP,	//0
    MODE_SOLVE, //1
    MODE_CLEANUP, //2
    MODE_STATUS, //3
    MODE_AGGREGATE, //4
    MODE_SETUP_MAXWELL, //5
    MODE_SOLVE_NEWMATRIX, //6
    MODE_ERROR, //7
	MODE_GET //8
} MODE_TYPE;

typedef enum
{
	MATRIX,
	MULTIVECTOR,
	LOVECTOR,
	SCALAR,
	UNKNOWN
} HierAttribType;

typedef Tpetra::Vector<>::node_type mm_node_t;
typedef Tpetra::Vector<>::local_ordinal_type mm_LocalOrd;
typedef Tpetra::Vector<>::global_ordinal_type mm_GlobalOrd;
typedef Tpetra::Map<> muemex_map_type;
typedef Tpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_CrsMatrix_double;
typedef std::complex<double> complex_t;
typedef Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_CrsMatrix_complex;
typedef Xpetra::Vector<mm_LocalOrd, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_ordinal_vector;
typedef Xpetra::Matrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_Matrix_double;
typedef Xpetra::Matrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_Matrix_complex;
typedef Xpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_MultiVector_double;
typedef Xpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Xpetra_MultiVector_complex;
typedef MueLu::Hierarchy<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Hierarchy_double;
typedef MueLu::Hierarchy<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Hierarchy_complex;

class muelu_data_pack
{
    public:
        muelu_data_pack(DataPackType type);	//type is one of EPETRA, TPETRA or TPETRA_COMPLEX
        virtual ~muelu_data_pack() = 0;
		virtual int status() = 0;
        virtual int setup(const mxArray* mxa) = 0;
        int id;
        Teuchos::RCP<Teuchos::ParameterList> List;
		DataPackType type;
		mxArray* getHierarchyData(std::string dataName, HierAttribType dataType, int levelID); //Works for all dp types
};

class muelu_epetra_data_pack : public muelu_data_pack
{
    public:
        muelu_epetra_data_pack();
        ~muelu_epetra_data_pack();
        int setup(const mxArray* mxa);
        int status();
        int solve(Teuchos::RCP<Teuchos::ParameterList> TPL, Teuchos::RCP<Epetra_CrsMatrix> Amat, double* b, double* x, int numVecs, int &iters);
        Teuchos::RCP<Epetra_CrsMatrix> GetMatrix()
        {
            return A;
        }
        Teuchos::RCP<Epetra_Operator> GetPrec()
        {
            return prec;
        }
        int NumGlobalRows()
        {
            return A->NumGlobalRows();
        }
        int NumMyCols()
        {
            return A->NumGlobalCols();
        }
        double operatorComplexity;
		Teuchos::RCP<Hierarchy_double> getHierarchy();
    private:
        Teuchos::RCP<Epetra_CrsMatrix> A;
        Teuchos::RCP<Epetra_Operator> prec;
};

//Scalar can be double or std::complex<double> (complex_t)
//Note: DataPackType is either TPETRA or TPETRA_COMPLEX
class muelu_tpetra_double_data_pack : public muelu_data_pack
{
    public:
        muelu_tpetra_double_data_pack();
        ~muelu_tpetra_double_data_pack();
        int setup(const mxArray* mxa);
        int status();
        int solve(Teuchos::RCP<Teuchos::ParameterList> TPL, Teuchos::RCP<Tpetra_CrsMatrix_double> Amat, double* b, double* x, int numVecs, int &iters);
        //note: I typedef'd mm_node_t at the top of this file as the Kokkos default type
        Teuchos::RCP<Tpetra_CrsMatrix_double> GetMatrix()
        {
            return A;
        }
        Teuchos::RCP<Tpetra::Operator<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> GetPrec()
        {
            return prec;
        }
        int NumMyRows()
        {
            if(A.is_null())
                return 0;
            else            
                return A->getNodeNumRows();
		}
        int NumMyCols()
        {
            if(A.is_null())
                return 0;
            else
                return A->getNodeNumCols();
        }
        double operatorComplexity;
        Teuchos::RCP<Hierarchy_double> getHierarchy();
    private:
        Teuchos::RCP<Tpetra_CrsMatrix_double> A;
        Teuchos::RCP<Tpetra::Operator<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> prec;
};

//Only define this datapack type if complex scalars are supported
#ifdef HAVE_COMPLEX_SCALARS

class muelu_tpetra_complex_data_pack : public muelu_data_pack
{
    public:
        muelu_tpetra_complex_data_pack();
        ~muelu_tpetra_complex_data_pack();
        int setup(const mxArray* mxa);
        int status();
        int solve(Teuchos::RCP<Teuchos::ParameterList> TPL, Teuchos::RCP<Tpetra_CrsMatrix_complex> Amat, complex_t* b, complex_t* x, int numVecs, int &iters);
        //note: I typedef'd mm_node_t at the top of this file as the Kokkos default type
        Teuchos::RCP<Tpetra_CrsMatrix_complex> GetMatrix()
        {
            return A;
        }
        Teuchos::RCP<Tpetra::Operator<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> GetPrec()
        {
            return prec;
        }
        int NumMyRows()
        {
            if(A.is_null())
                return 0;
            else
                return A->getNodeNumRows();
		}
        int NumMyCols()
        {
            if(A.is_null())
                return 0;
            else
                return A->getNodeNumCols();
        }
        double operatorComplexity;
       	Teuchos::RCP<Hierarchy_complex> getHierarchy();
    private:
        Teuchos::RCP<Tpetra_CrsMatrix_complex> A;
        Teuchos::RCP<Tpetra::Operator<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> prec;
};
#endif	//complex scalars

/* nonmember utility functions */

//create an array of ints from an array of MATLAB array indices 
int* mwIndex_to_int(int N, mwIndex* mwi_array);
//Convert individual Belos verbosity setting name to its enum value
int strToMsgType(const char* str);
//attempt to get hierarchy data type from the name of the field (A, P, Nullspace, etc)
HierAttribType strToHierAttribType(const char* str);
//Parse belos output style (Brief or General)
int strToOutputStyle(const char* str);
//Parse belos verbosity settings (returns enum values | together)
int getBelosVerbosity(const char* input);
//Get an int from a MATLAB double or int input
int parseInt(const mxArray* mxa);
//create a sparse array in Matlab
template<typename Scalar = double>
mxArray* createMatlabSparse(int numRows, int numCols, int nnz);
//create an int32 dense vector in Matlab
mxArray* createMatlabLOVector(Teuchos::RCP<Xpetra_ordinal_vector> vec);
//copy a sparse Xpetra matrix (double or complex) to Matlab
template<typename Scalar = double>
mxArray* saveMatrixToMatlab(Teuchos::RCP<Xpetra::Matrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mat);
template<typename Scalar = double>
mxArray* createMatlabMultiVector(int numRows, int numCols);
template<typename Scalar = double>
mxArray* saveMultiVectorToMatlab(Teuchos::RCP<Xpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mv);
template<typename Scalar = double>
void fillMatlabArray(Scalar* array, const mxArray* mxa, int n);
//Set up an Epetra matrix from CSC arrays
Teuchos::RCP<Epetra_CrsMatrix> epetra_setup(int Nrows, int Ncols, int* rowind, int* colptr, double* vals);
//Set up an Epetra matrix from a MATLAB array
Teuchos::RCP<Epetra_CrsMatrix> epetra_setup_from_prhs(const mxArray* mxa);
//Set up an Epetra_MultiVector from MATLAB array
Teuchos::RCP<Epetra_MultiVector> epetra_setup_multivector(const mxArray* mxa);
//Set up Tpetra real matrix from MATLAB array
Teuchos::RCP<Tpetra_CrsMatrix_double> tpetra_setup_real_prhs(const mxArray* mxa);
//Set up Xpetra real matrix from MATLAB array
Teuchos::RCP<Xpetra::Matrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> xpetra_setup_real(const mxArray* mxa);
//Set up Tpetra real multivector from MATLAB array
Teuchos::RCP<Tpetra::MultiVector<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> tpetra_setup_real_multivector(const mxArray* mxa);
//Set up tpetra complex matrix from MATLAB
Teuchos::RCP<Tpetra_CrsMatrix_complex> tpetra_setup_complex_prhs(const mxArray* mxa);
//Set up Tpetra complex multivector from MATLAB
Teuchos::RCP<Tpetra::MultiVector<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> tpetra_setup_complex_multivector(const mxArray* mxa);
//Set up xpetra complex matrix from MATLAB
Teuchos::RCP<Xpetra::Matrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> xpetra_setup_complex(const mxArray* mxa);
//Solve an Epetra system
int epetra_solve(Teuchos::RCP<Teuchos::ParameterList> SetupList, Teuchos::RCP<Teuchos::ParameterList> TPL, Teuchos::RCP<Epetra_CrsMatrix> A, Teuchos::RCP<Epetra_Operator> prec, double* b, double* x, int numVecs, int &iters);
//Solve a real-valued Tpetra system
int tpetra_double_solve(Teuchos::RCP<Teuchos::ParameterList> SetupList, Teuchos::RCP<Teuchos::ParameterList> TPL, Teuchos::RCP<Tpetra_CrsMatrix_double> A, Teuchos::RCP<Tpetra::Operator<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> prec, double* b, double* x, int numVecs, int& iters);
//Solve a complex-valued Tpetra system
int tpetra_complex_solve(Teuchos::RCP<Teuchos::ParameterList> SetupList, Teuchos::RCP<Teuchos::ParameterList> TPL,
Teuchos::RCP<Tpetra_CrsMatrix_complex> A, Teuchos::RCP<Tpetra::Operator<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> prec, complex_t* b, complex_t* x, int numVecs, int& iters);
//Get a hierarchy from a muelu_data_pack
template<typename Scalar = double>
Teuchos::RCP<MueLu::Hierarchy<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> getDatapackHierarchy(muelu_data_pack* dp);

namespace muelu_data_pack_list
{
	extern std::vector<Teuchos::RCP<muelu_data_pack>> list;
	extern int nextID;
	int add(Teuchos::RCP<muelu_data_pack> D);
	Teuchos::RCP<muelu_data_pack> find(int id);
	int remove(int id);
	int size();
	int status_all();
	bool isInList(int id);
	void clearAll();
}

#endif //MUEMEX_H
