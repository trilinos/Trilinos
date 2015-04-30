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

#ifdef HAVE_MUELU_MATLAB
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
    EPETRA_UNPREC,
    EPETRA,
    TPETRA,
    TPETRA_COMPLEX
} DataPackType;

typedef enum
{
    MODE_SETUP,
    MODE_SOLVE,
    MODE_CLEANUP,
    MODE_STATUS,
    MODE_AGGREGATE,
    MODE_SETUP_MAXWELL,
    MODE_SOLVE_NEWMATRIX,
    MODE_ERROR
} MODE_TYPE;

//Default Tpetra node type for CrsMatrix (could replace here with custom types)
typedef Tpetra::Vector<>::node_type mm_node_t;
typedef Tpetra::Vector<>::local_ordinal_type mm_LocalOrd;
typedef Tpetra::Vector<>::global_ordinal_type mm_GlobalOrd;
typedef Tpetra::Map<> muemex_map_type;
typedef Tpetra::CrsMatrix<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_CrsMatrix_double;

//Scalar could be double or std::complex, and Matrix could be the CrsMatrix from [E|T]petra
//(therefore the muelu_epetra_data_pack has the default template arguments)
class muelu_data_pack
{
    public:
        muelu_data_pack(DataPackType type);
        virtual ~muelu_data_pack();
        virtual int status() = 0;
        virtual int setup(const mxArray* mxa) = 0;
        virtual int NumMyRows() = 0;
        virtual int NumMyCols() = 0;
        int id;
        Teuchos::RCP<Teuchos::ParameterList> List;
        DataPackType type;
};

class muelu_epetra_unprec_data_pack : public muelu_data_pack
{
    public:
        muelu_epetra_unprec_data_pack();
        ~muelu_epetra_unprec_data_pack();
        int setup(const mxArray* mxa);
        int status();
        int solve(Teuchos::RCP<Teuchos::ParameterList> TPL, Teuchos::RCP<Epetra_CrsMatrix> Amat, double* b, double* x, int numVecs, int &iters);
        Teuchos::RCP<Epetra_CrsMatrix> GetMatrix()
        {
            return A;
        }
        int NumMyRows()
        {
            return A->NumMyRows();
        }
        int NumMyCols()
        {
            return A->NumMyCols();
        }
    private:
        Teuchos::RCP<Epetra_CrsMatrix> A;
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
        int NumMyRows()
        {
            return A->NumMyRows();
        }
        int NumMyCols()
        {
            return A->NumMyCols();
        }
        double operatorComplexity;
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
    private:
        Teuchos::RCP<Tpetra_CrsMatrix_double> A;
        Teuchos::RCP<Tpetra::Operator<double, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> prec;
};

//Only define this datapack type if complex scalars are supported
#ifdef HAVE_COMPLEX_SCALARS
typedef std::complex<double> complex_t;
typedef Tpetra::CrsMatrix<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t> Tpetra_CrsMatrix_complex;

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
    private:
        Teuchos::RCP<Tpetra_CrsMatrix_complex> A;
        Teuchos::RCP<Tpetra::Operator<complex_t, mm_LocalOrd, mm_GlobalOrd, mm_node_t>> prec;
};
#endif

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

#endif
#endif
