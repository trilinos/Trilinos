#ifndef PIRO_EPETRA_NOXOBSERVER
#define PIRO_EPETRA_NOXOBSERVER

#include "Epetra_Vector.h"

namespace Piro {
namespace Epetra {

class NOXObserver 
{
public:
  NOXObserver () {};
  virtual ~NOXObserver () {};

  virtual void observeSolution(const Epetra_Vector& soln) = 0;

private:

};

}
}

#endif //PIRO_EPETRA_NOXOBSERVER
