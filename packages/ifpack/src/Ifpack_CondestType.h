#ifndef IFPACK_CONDESTTYPE_H
#define IFPACK_CONDESTTYPE_H

//! Ifpack_CondestType: enum to define the type of condition number estimate.

enum Ifpack_CondestType {
  Ifpack_Cheap,  //!< cheap estimate
  Ifpack_CG,     //!< Uses AztecOO's CG
  Ifpack_GMRES   //!< Uses AztecOO's GMRES
};

#endif
