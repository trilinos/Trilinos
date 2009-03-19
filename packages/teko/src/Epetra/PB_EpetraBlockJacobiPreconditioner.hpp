#ifndef __PB_EpetraBlockJacobiPreconditioner_hpp__
#define __PB_EpetraBlockJacobiPreconditioner_hpp__

#include "PB_EpetraInverseOpWrapper.hpp"

#include "Epetra_Vector.h"

namespace PB {
namespace Epetra {

class EpetraBlockJacobiPreconditioner : public EpetraInverseOpWrapper {
public:
    EpetraBlockJacobiPreconditioner(const Epetra_Operator * A,
                                  const Epetra_Operator * invD1,
                                  const Epetra_Operator * invD2);

}; // end class EpetraBlockJacobiPreconditioner

}
} // end PB

#endif
