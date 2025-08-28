// @HEADER
// *****************************************************************************
//                 Anasazi: Block Eigensolvers Package
//
// Copyright 2004 NTESS and the Anasazi contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This software is a result of the research described in the report
//
//     "A comparison of algorithms for modal analysis in the absence
//     of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//     Sandia National Laboratories, Technical report SAND2003-1028J.
//
// It is based on the Epetra, AztecOO, and ML packages defined in the Trilinos
// framework ( http://trilinos.org/ ).

#ifndef MODE_LAPLACE_1D_Q2_H
#define MODE_LAPLACE_1D_Q2_H

#include "Epetra_ConfigDefs.h"

#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include "CheckingTools.h"
#include "ModeLaplace.h"
#include "SortingTools.h"

class ModeLaplace1DQ2 : public ModeLaplace {

  private:

    const CheckingTools myVerify;
    const Epetra_Comm &MyComm;
    const SortingTools mySort;

    Epetra_Map *Map;
    Epetra_Operator *K;
    Epetra_Operator *M;

    double Lx;
    int nX;

    double *x;

    static const int dofEle;
    static const int maxConnect;
#ifndef M_PI
    static const double M_PI;
#endif

    // Private member functions
    void preProcess();
    void makeMap();
    int countElements(bool *isTouched);
    void makeMyElementsTopology(int *elemTopo, bool *isTouched);
    void makeMyConnectivity(int *elemTopo, int numEle, int *connectivity, int *numNz);
    void makeStiffness(int *elemTopo, int numEle, int *connectivity, int *numNz);
    void makeElementaryStiffness(double *kel) const;
    void makeMass(int *elemTopo, int numEle, int *connectivity, int *numNz);
    void makeElementaryMass(double *mel) const;

    // Don't define these functions
    ModeLaplace1DQ2(const ModeLaplace1DQ2 &ref);
    ModeLaplace1DQ2& operator=(const ModeLaplace1DQ2 &ref);

  public:

    ModeLaplace1DQ2(const Epetra_Comm &_Comm, double _Lx, int _nX);

    ~ModeLaplace1DQ2();

    const Epetra_Operator* getStiffness() const { return K; }
    const Epetra_Operator* getMass()      const { return M; }

    double getFirstMassEigenValue() const;

    int eigenCheck(const Epetra_MultiVector &Q, double *lambda, double *normWeight) const;

    void memoryInfo() const;
    void problemInfo() const;

};

#endif
