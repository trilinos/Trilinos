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

#include <stdio.h>
#include <string>
#include <vector>

#include "MueLu_config.hpp"
#include "MueLu.hpp"

/*Epetra headers*/
#include "Epetra_SerialComm.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"

/* TODO: Tpetra hdrs */

#include "Teuchos_ParameterList.hpp"
#include "BelosSolverFactory.hpp"
#include "BelosEpetraAdapter.hpp"

#ifdef HAVE_MUELU_MATLAB
#include "mex.h"

class muelu_data_pack;
class muelu_data_pack
{
    public:
        muelu_data_pack();
        virtual ~muelu_data_pack();
        virtual int setup(int N, int* rowind, int* colptr, double* vals) = 0;
        virtual int status() = 0;
        virtual int solve(Teuchos::ParameterList* TPL, Epetra_CrsMatrix* A, double* b, double* x, int &iters) = 0;
        virtual int NumMyRows() = 0;
        virtual int NumMyCols() = 0;
        virtual Epetra_CrsMatrix* GetMatrix() = 0;
        int id;
        Teuchos::ParameterList* List;
        double operator_complexity;
        muelu_data_pack* next;
};


//Temporary, pretend mueluapi_data_pack is a muelu_epetra_data_pack

#define mueluapi_data_pack muelu_epetra_data_pack

/*
class mueluapi_data_pack : public muelu_data_pack
{
    public:
        mueluapi_data_pack();
        ~mueluapi_data_pack();
        int setup(int N, int* rowind, int* colptr, double* vals);
        int status();
        int NumMyRows()
        {
            return A->NumMyRows();
        }
        int NumMyCols()
        {
            return A->NumMyCols();
        }
    private:
        Epetra_CrsMatrix* GetMatrix()
        {
            return 0;
        }
        //TODO: I assume MueLu multigrid needs some stuff here?
};
*/

class muelu_epetra_data_pack : public muelu_data_pack
{
    public:
        muelu_epetra_data_pack();
        ~muelu_epetra_data_pack();
        int setup(int N, int* rowind, int* colptr, double* vals);
        int status();
        int solve(Teuchos::ParameterList *TPL, Epetra_CrsMatrix *Amat, double*b, double*x,int &iters);
        Epetra_CrsMatrix* GetMatrix()
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
        Epetra_CrsMatrix* A;
};

/**************************************************************/
/**************************************************************/
/**************************************************************/
/* MueLu data pack list */
class muelu_data_pack_list
{
public:
  muelu_data_pack_list();
  ~muelu_data_pack_list();

  /* add - Adds an MueLu_DATA_PACK to the list.
     Parameters:
     D       - The MueLu_DATA_PACK. [I]
     Returns: problem id number of D
  */
  int add(muelu_data_pack *D);

  /* find - Finds problem by id
     Parameters:
     id      - ID number [I]
     Returns: pointer to MueLu_DATA_PACK matching 'id', if found, NULL if not
     found.
  */
  muelu_data_pack* find(int id);

  /* remove - Removes problem by id
     Parameters:
     id      - ID number [I]
     Returns: IS_TRUE if remove was succesful, IS_FALSE otherwise
  */
  int remove(int id);

  /* size - Number of stored problems
     Returns: num_probs
  */
  int size();

  /* Returns the status of all members of the list
     Returns IS_TRUE
  */
  int status_all();

protected:
  int num_probs;
  /* Note: This list is sorted */
  muelu_data_pack *L;
};

/*

TODO: Implement the Tpetra option with this datapack type
class muelu_tpetra_data_pack : public muelu_data_pack
{
    public:
        muelu_tpetra_data_pack();
        ~muelu_tpetra_data_pack();
        int setup(int N, int* rowind, int* colptr, double* vals);
        int status();
        Tpetra_CrsMatrix* GetMatrix()
        {
            return A;
        }
        int NumMyRows()
        {
            return A->
    private:
        Tpetra::CrsMatrix<>* A;
};
*/

#endif

#endif
