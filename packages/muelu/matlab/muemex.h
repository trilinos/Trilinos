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
#include "muemexCallbacks.h"
#include "muemexTypes.h"

#define HAVE_COMPLEX_SCALARS

typedef enum
{
    EPETRA,
    TPETRA,
    TPETRA_COMPLEX
} DataPackType;

//Mode value passed to MATLAB muemex function as 1st arg (int)
typedef enum
{
    MODE_SETUP,	//0
    MODE_SOLVE, //1
    MODE_CLEANUP, //2
    MODE_STATUS, //3
    MODE_AGGREGATE, //4
    MODE_GET, //5
    MODE_ERROR
} MODE_TYPE;

typedef enum
{
	MATRIX,
	MULTIVECTOR,
	LOVECTOR,
	SCALAR,
	UNKNOWN
} HierAttribType;

class MuemexSystem
{
    public:
        MuemexSystem(DataPackType type);	//type is one of EPETRA, TPETRA or TPETRA_COMPLEX
        virtual ~MuemexSystem() = 0;
		virtual int status() = 0;
        virtual int setup(const mxArray* mxa) = 0;
        int id;
        Teuchos::RCP<Teuchos::ParameterList> List;
		DataPackType type;
		mxArray* getHierarchyData(std::string dataName, HierAttribType dataType, int levelID); //Works for all dp types
};

class EpetraSystem : public MuemexSystem
{
    public:
        EpetraSystem();
        ~EpetraSystem();
        int setup(const mxArray* mxa);
        int status();
        mxArray* solve(Teuchos::RCP<Teuchos::ParameterList> params, Teuchos::RCP<Epetra_CrsMatrix> matrix, const mxArray* rhs, int &iters);
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

template<typename Scalar>
class TpetraSystem : public MuemexSystem
{
    public:
        TpetraSystem();
        ~TpetraSystem();
        typedef Tpetra::CrsMatrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t> TMatrix;
        typedef Tpetra::Operator<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t> TOperator;
        int setup(const mxArray* mxa);
        int status();
        mxArray* solve(Teuchos::RCP<Teuchos::ParameterList> params, Teuchos::RCP<TMatrix> matrix, const mxArray* rhs, int &iters);
        //note: I typedef'd mm_node_t at the top of this file as the Kokkos default type
        Teuchos::RCP<TMatrix> GetMatrix()
        {
            return A;
        }
        Teuchos::RCP<TOperator> GetPrec()
        {
            return prec;
        }
        int NumMyRows()
        {
            if(A.is_null())
                return 0;
            else            
                return A->getGlobalNumRows();
		}
        int NumMyCols()
        {
            if(A.is_null())
                return 0;
            else
                return A->getGlobalNumCols();
        }
        double operatorComplexity;
        Teuchos::RCP<MueLu::Hierarchy<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> getHierarchy();
    private:
        Teuchos::RCP<TMatrix> A;
        Teuchos::RCP<TOperator> prec;
};

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
//Load a multivector from a MATLAB array
template<typename Scalar>
Teuchos::RCP<Tpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> loadTpetraMV(const mxArray* mxa);
//create a sparse array in Matlab
template<typename Scalar>
mxArray* createMatlabSparse(int numRows, int numCols, int nnz);
//create an ordinal (int32) vector in Matlab
mxArray* createMatlabLOVector(Teuchos::RCP<Xpetra_ordinal_vector> vec);
//copy a sparse Xpetra matrix (double or complex) to Matlab
template<typename Scalar>
mxArray* saveMatrixToMatlab(Teuchos::RCP<Xpetra::Matrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mat);
template<typename Scalar>
mxArray* createMatlabMultiVector(int numRows, int numCols);
template<typename Scalar>
mxArray* saveTpetraMV(Teuchos::RCP<Tpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mv);
template<typename Scalar>
mxArray* saveMultiVectorToMatlab(Teuchos::RCP<Xpetra::MultiVector<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> mv);
template<typename Scalar>
void fillMatlabArray(Scalar* array, const mxArray* mxa, int n);
//set up Tpetra matrix from MATLAB array
template<typename Scalar>
Teuchos::RCP<Tpetra::CrsMatrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> tpetraLoadMatrix(const mxArray* mxa);
//same as above but for Xpetra
template<typename Scalar>
Teuchos::RCP<Xpetra::Matrix<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> xpetraLoadMatrix(const mxArray* mxa);
//Get a hierarchy from a MuemexSystem
template<typename Scalar>
Teuchos::RCP<MueLu::Hierarchy<Scalar, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> getDatapackHierarchy(MuemexSystem* dp);
//Get an Epetra_MultiVector from MATLAB array
Teuchos::RCP<Epetra_MultiVector> loadEpetraMV(mxArray* mxa);
//Save an Epetra MV to MATLAB array
mxArray* saveEpetraMV(Teuchos::RCP<Epetra_MultiVector> mv);
//Load an Epetra matrix from MATLAB array
Teuchos::RCP<Epetra_CrsMatrix> epetraLoadMatrix(const mxArray* mxa);

namespace MuemexSystemList
{
	extern std::vector<Teuchos::RCP<MuemexSystem>> list;
	extern int nextID;
	int add(Teuchos::RCP<MuemexSystem> D);
	Teuchos::RCP<MuemexSystem> find(int id);
	int remove(int id);
	int size();
	int status_all();
	bool isInList(int id);
	void clearAll();
}

#endif //MUEMEX_H
