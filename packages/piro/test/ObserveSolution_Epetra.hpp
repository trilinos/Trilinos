// @HEADER
// *****************************************************************************
//        Piro: Strategy package for embedded analysis capabilitites
//
// Copyright 2010 NTESS and the Piro contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PIROTEST_OBSERVESOLUTION_EPETRA
#define PIROTEST_OBSERVESOLUTION_EPETRA


#include "NOX_Epetra_Observer.H"

class ObserveSolution_Epetra : public NOX::Epetra::Observer
{
public:
   ObserveSolution_Epetra () {};

   ~ObserveSolution_Epetra () { };

  void observeSolution(
    const Epetra_Vector& solution)
    {
      double norm; solution.Norm2(&norm);
      if (solution.Comm().MyPID()==0)
        std::cout << "ObserveSolution: Norm = " << norm << std::endl;
    }

  void observeSolution(
    const Epetra_Vector& solution, double time_or_param_val)
    {
      double norm; solution.Norm2(&norm);
      if (solution.Comm().MyPID()==0)
        std::cout << "ObserveSolution: Norm = " << norm 
             << "  for param/time = " << time_or_param_val << std::endl;
    }

private:

};

#endif //PIROTEST_OBSERVESOLUTION_EPETRA
