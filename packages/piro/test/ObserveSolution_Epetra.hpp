#ifndef PIROTEST_OBSERVESOLUTION_EPETRA
#define PIROTEST_OBSERVESOLUTION_EPETRA


#include "Piro_Epetra_NOXObserver.hpp"

class ObserveSolution_Epetra : public Piro::Epetra::NOXObserver
{
public:
   ObserveSolution_Epetra () {};

   ~ObserveSolution_Epetra () { };

  void observeSolution(
    const Epetra_Vector& solution)
    {
      double norm; solution.Norm2(&norm);
      if (solution.Comm().MyPID()==0)
        cout << "ObserveSolution: Norm = " << norm << endl;
    }

private:

};

#endif //PIROTEST_OBSERVESOLUTION_EPETRA
