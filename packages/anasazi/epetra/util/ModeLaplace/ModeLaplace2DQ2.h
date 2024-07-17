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

#ifndef ANASAZI_MODE_LAPLACE_2D_Q2_H
#define ANASAZI_MODE_LAPLACE_2D_Q2_H

#include "Epetra_ConfigDefs.h"
#include "Anasaziepetra_ModeLaplace_DLLExportMacro.h"

#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MultiVector.h"
#include "Epetra_RowMatrix.h"

#include "CheckingTools.h"
#include "ModeLaplace.h"
#include "SortingTools.h"


// Chaco partition routine
#ifdef _USE_CHACO
extern "C" {
  int interface(int, int*, int*, int *, float*, float*, float*, float*, char*,
                char*, short int*, int, int, int[3], double*, int, int, int, int,
                int, double, long);
}
#endif


class ANASAZIEPETRA_MODELAPLACE_LIB_DLL_EXPORT ModeLaplace2DQ2 : public ModeLaplace {

  private:

    const CheckingTools myVerify;
    const Epetra_Comm &MyComm;
    const SortingTools mySort;

    Epetra_Map *Map;
    Epetra_CrsMatrix *K;
    Epetra_CrsMatrix *M;

    double Lx;
    int nX;

    double Ly;
    int nY;

    double *x;
    double *y;

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
    ModeLaplace2DQ2(const ModeLaplace2DQ2 &ref);
    ModeLaplace2DQ2& operator=(const ModeLaplace2DQ2 &ref);

  public:

    ModeLaplace2DQ2(const Epetra_Comm &_Comm, double _Lx, int _nX, double _Ly, int _nY);

    ~ModeLaplace2DQ2();

    const Epetra_CrsMatrix* getStiffness() const { return K; }
    const Epetra_CrsMatrix* getMass()      const { return M; }

    double getFirstMassEigenValue() const;

    int eigenCheck(const Epetra_MultiVector &Q, double *lambda, double *normWeight, bool smallest = true) const;

    void memoryInfo() const;
    void problemInfo() const;

};

#endif
