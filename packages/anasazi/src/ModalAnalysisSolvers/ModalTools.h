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

#ifndef MODAL_TOOLS_H
#define MODAL_TOOLS_H

#include "Epetra_ConfigDefs.h"

#include "Epetra_BLAS.h"
#include "Epetra_Comm.h"
#include "Epetra_LAPACK.h"
#include "Epetra_LocalMap.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Operator.h"
#include "Epetra_Time.h"
#include "Epetra_Vector.h"

#include "FortranRoutines.h"

class ModalTools {

  private:

    const FortranRoutines callFortran;
    const Epetra_BLAS callBLAS;
    const Epetra_LAPACK callLAPACK;

    const Epetra_Comm &MyComm;
    const Epetra_Time MyWatch;

    double eps;

    double timeQtMult;
    double timeQMult;
    double timeProj_MassMult;
    double timeNorm_MassMult;
    double timeProj;
    double timeNorm;

    int numProj_MassMult;
    int numNorm_MassMult;

  public:

    ModalTools(const Epetra_Comm &_Comm);

    int makeSimpleLumpedMass(const Epetra_Operator *M, double *weight) const;

    int massOrthonormalize(Epetra_MultiVector &X, Epetra_MultiVector &MX,
                           const Epetra_Operator *M, const Epetra_MultiVector &Q, int howMany,
                           int type = 0, double *WS = 0, double kappa = 1.5625);

    void localProjection(int numRow, int numCol, int length,
                         double *U, int ldU, double *MatV, int ldV, 
                         double *UtMatV, int ldUtMatV, double *work) const;

    int directSolver(int, double*, int, double*, int, int&, double*, int, double*, int,
                     int = 0) const;

    double getTimeProj()  const { return timeProj; }
    double getTimeProj_QtMult()    const { return timeQtMult; }
    double getTimeProj_QMult()     const { return timeQMult; }
    double getTimeProj_MassMult()  const { return timeProj_MassMult; }
    int getNumProj_MassMult()      const { return numProj_MassMult; }

    double getTimeNorm()  const { return timeNorm; }
    double getTimeNorm_MassMult() const { return timeNorm_MassMult; }
    int getNumNorm_MassMult()     const { return numNorm_MassMult; }

};


#endif

