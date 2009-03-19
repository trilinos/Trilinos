#ifndef __PB_EpetraLSCStablePreconditioner_hpp__
#define __PB_EpetraLSCStablePreconditioner_hpp__

#include "PB_EpetraInverseOpWrapper.hpp"

#include "Epetra_Vector.h"

namespace PB {
namespace Epetra {

class EpetraLSCStablePreconditioner : public EpetraInverseOpWrapper {
public:
    EpetraLSCStablePreconditioner(const Epetra_Operator * A,
                                  const Epetra_Operator * invF,
                                  const Epetra_Operator * invS,
                                  const Epetra_Vector * invMass=0);

}; // end class EpetraLSCStablePreconditioner

}
} // end PB

#endif
