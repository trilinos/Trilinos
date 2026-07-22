// @HEADER
// *****************************************************************************
//                    Teuchos: Common Tools Package
//
// Copyright 2004 NTESS and the Teuchos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TEUCHOS_MPI_REDUCTION_OP_SETTER_HPP
#define TEUCHOS_MPI_REDUCTION_OP_SETTER_HPP

/// \file Teuchos_MpiReductionOpSetter.hpp
/// \brief Implementation detail of Teuchos' MPI wrapper.
///
/// \warning Everything in this file is an implementation detail of
///   Teuchos.  Do not use it and do not depend on it.
///
/// This file only contains meaningful content if building with MPI.

#include <Teuchos_ConfigDefs.hpp>

#ifdef HAVE_MPI
#include "Teuchos_ReductionOp.hpp"
#include <mpi.h>

namespace Teuchos {
namespace Details {

/// \brief Base class for an MPI-compatible reduction operator.
///
/// This class exists as an intermediate wrapper between
/// ValueTypeReductionOp and MPI_Op.  Its pure virtual method reduce()
/// takes the same arguments as an MPI_Op reduction / scan operator.
///
/// \note This class only exists if building with MPI enabled.
class TEUCHOSCOMM_LIB_DLL_EXPORT MpiReductionOpBase :
    virtual public Describable {
public:
  virtual void
  reduce (void* invec, void* inoutvec,
          int* len, MPI_Datatype* datatype) const = 0;
};

/// \brief Subclass of MpiReductionOpBase that implements reduce()
///   using a ValueTypeReductionOp instance.
///
/// \note This class only views the ValueTypeReductionOp instance; it
///   does not own it.  The ValueTypeReductionOp must be in scope;
///   otherwise this class' behavior is undefined.
template<typename OrdinalType>
class MpiReductionOp : public MpiReductionOpBase {
public:
  MpiReductionOp (const ValueTypeReductionOp<OrdinalType,char>& reductOp)
    : reductOp_ (reductOp)
  {}

  void
  reduce (void* invec, void* inoutvec, int* len, MPI_Datatype* datatype) const
  {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(!len);
    TEUCHOS_TEST_FOR_EXCEPT(!datatype);
#endif
    // mfh 23 Nov 2014: Ross made the unfortunate decision initially
    // to mash everything into MPI_CHAR.  This means that we have to
    // do a lot of type casting here.  Of course, most implementations
    // of ValueTypeReductionOp will immediately cast back from char*
    // to the actual type of interest.
    int sz;
    MPI_Type_size (*datatype, &sz);
    (void) datatype;
    reductOp_.reduce ((*len) * sz, reinterpret_cast<char*> (invec),
                      reinterpret_cast<char*> (inoutvec));
  }

private:
  const ValueTypeReductionOp<OrdinalType, char>& reductOp_;
  // Not defined and not to be called
  MpiReductionOp ();
  MpiReductionOp (const MpiReductionOp&);
  MpiReductionOp& operator= (const MpiReductionOp&);
};

/// \brief Set the current reduction or scan operation for MpiComm.
///   Get the resulting MPI_Op to pass into MPI functions.
///
/// \return Global singleton MPI_Op to pass into MPI functions.
///   This is valid only before MPI_Finalize has been called.
///
/// \warning This is an implementation detail of Teuchos.
///   Users should never call this function directly.
TEUCHOSCOMM_LIB_DLL_EXPORT MPI_Op setMpiReductionOp (const MpiReductionOpBase& reductOp);

} // namespace Details
} // namespace Teuchos

#endif // HAVE_MPI

#endif // TEUCHOS_MPI_REDUCTION_OP_SETTER_HPP
