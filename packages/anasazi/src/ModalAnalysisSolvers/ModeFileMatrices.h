//**************************************************************************
//
//                                 NOTICE
//
// This software is a result of the research described in the report
//
// " A comparison of algorithms for modal analysis in the absence 
//   of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//  Sandia National Laboratories, Technical report SAND2003-1028J.
//
// It is based on the Epetra, AztecOO, and ML packages defined in the Trilinos
// framework ( http://software.sandia.gov/trilinos/ ).
//
// The distribution of this software follows also the rules defined in Trilinos.
// This notice shall be marked on any reproduction of this software, in whole or
// in part.
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//
// Code Authors: U. Hetmaniuk (ulhetma@sandia.gov), R. Lehoucq (rblehou@sandia.gov)
//
//**************************************************************************

#ifndef MODE_FILE_MATRICES_H
#define MODE_FILE_MATRICES_H

#include <fstream>

#include "Epetra_ConfigDefs.h"

#include "Epetra_Comm.h"
#include "Epetra_CrsGraph.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_Vector.h"

#include "ModalProblem.h"

// Chaco partition routine
#ifdef _USE_CHACO
extern "C" {
  int interface(int, int*, int*, int *, float*, float*, float*, float*, char*,
                char*, short int*, int, int, int[3], double*, int, int, int, int,
                int, double, long);
}
#endif

class ModeFileMatrices : public ModalProblem {

  private:

    const Epetra_Comm &MyComm;

    Epetra_Map *Map;

    Epetra_Operator *K;
    Epetra_Operator *M;

    const char *stiffnessFileName;
    const char *massFileName;

    double kScale;
    double mScale;

    double shift;

    // Private member functions
    void preProcess();
    void fillMatrices();

    // Don't define these functions
    ModeFileMatrices(const ModeFileMatrices &ref);
    ModeFileMatrices& operator=(const ModeFileMatrices &ref);

    public:

    ModeFileMatrices(const Epetra_Comm &_Comm, const char *_stif, const char *_mass);

    ~ModeFileMatrices();

    const Epetra_Operator* getStiffness() const { return K; }
    const Epetra_Operator* getMass()      const { return M; }

    int eigenCheck(const Epetra_MultiVector &Q, double *lambda, double *normWeight) const;

    void memoryInfo() const;
    void problemInfo() const;

};

#endif
