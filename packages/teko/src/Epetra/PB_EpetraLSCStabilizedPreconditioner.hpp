#ifndef __PB_EpetraLSCStabilizedPreconditioner_hpp__
#define __PB_EpetraLSCStabilizedPreconditioner_hpp__

#include "PB_EpetraInverseOpWrapper.hpp"

#include "Epetra_Vector.h"

namespace PB {
namespace Epetra {

class EpetraLSCStabilizedPreconditioner : public EpetraInverseOpWrapper {
public:
    EpetraLSCStabilizedPreconditioner(const Epetra_Operator * A,
                                  const Epetra_Operator * invF,
                                  const Epetra_Operator * invBQBtmC,
                                  const Epetra_Vector * invD,
                                  const Epetra_Vector * invMass=0);
}; // end class EpetraLSCStabilizedPreconditioner

}
} // end PB

#endif
