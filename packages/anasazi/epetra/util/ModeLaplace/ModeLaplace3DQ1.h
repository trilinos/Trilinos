// @HEADER
// ***********************************************************************
//
//                 Anasazi: Block Eigensolvers Package
//                 Copyright 2004 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

// This software is a result of the research described in the report
//
//     "A comparison of algorithms for modal analysis in the absence
//     of a sparse direct method", P. Arbenz, R. Lehoucq, and U. Hetmaniuk,
//     Sandia National Laboratories, Technical report SAND2003-1028J.
//
// It is based on the Epetra, AztecOO, and ML packages defined in the Trilinos
// framework ( http://trilinos.org/ ).

#ifndef ANASAZI_MODE_LAPLACE_3D_Q1_H
#define ANASAZI_MODE_LAPLACE_3D_Q1_H

#include "Epetra_ConfigDefs.h"

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

class ModeLaplace3DQ1 : public ModeLaplace {

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

    double Lz;
    int nZ;

    double *x;
    double *y;
    double *z;

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
    void makeMass(int *elemTopo, int numEle, int *connectivity, int *numNz);

    // Don't define these functions
    ModeLaplace3DQ1(const ModeLaplace3DQ1 &ref);
    ModeLaplace3DQ1& operator=(const ModeLaplace3DQ1 &ref);

  public:

    ModeLaplace3DQ1(const Epetra_Comm &_Comm, double _Lx, int _nX, double _Ly, int _nY,
                    double _Lz, int _nZ);

    ~ModeLaplace3DQ1();

    const Epetra_CrsMatrix* getStiffness() const { return K; }
    const Epetra_CrsMatrix* getMass()      const { return M; }

    double getFirstMassEigenValue() const;

    int eigenCheck(const Epetra_MultiVector &Q, double *lambda, double *normWeight, bool smallest = true) const;

    void memoryInfo() const;
    void problemInfo() const;

};

#endif
