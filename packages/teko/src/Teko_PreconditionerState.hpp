/*
// @HEADER
//
// ***********************************************************************
//
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
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
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
//
// ***********************************************************************
//
// @HEADER

*/

#ifndef __Teko_PreconditionerState_hpp__
#define __Teko_PreconditionerState_hpp__

#include "Teuchos_ParameterListAcceptor.hpp"

// Thyra includes
// #include "Thyra_SolveSupportTypes.hpp"

// Teko includes
#include "Teko_Utilities.hpp"

namespace Teko {

/** \brief An implementation of a state object
 *        preconditioners.
 *
 * An implementation of a state object for
 * preconditioners. This implementation takes a parameter list
 * and has capacity to store linear operators.
 * However, it is easy to override this class and fill it
 * with what ever is neccessary to build the preconditioner.
 * This is essentially a "bag" that can be filled with whatever
 * a user wants.
 */
class PreconditionerState : public Teuchos::ParameterListAcceptor {
 public:
  //! \name Default and copy constructors
  //@{
  PreconditionerState() : isInitialized_(false) {}
  PreconditionerState(const PreconditionerState& src)
      : paramList_(Teuchos::rcp(new Teuchos::ParameterList(*src.paramList_))),
        srcVector_(src.srcVector_),
        inverses_(src.inverses_),
        isInitialized_(src.isInitialized_) {}
  //@}

  //! \name for ParameterListAcceptor
  //@{

  //! Set parameters from a parameter list and return with default values.
  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList>& paramList);

  //! Get the parameter list that was set using setParameterList().
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();

  //! Unset the parameter list that was set using setParameterList().
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();

  //@}

  /** Has this state object been initialized for the particular operator
   * being used?
   */
  virtual bool isInitialized() const { return isInitialized_; }

  /** Set that this state object has been initialized for this operator.
   */
  virtual void setInitialized(bool init = true) { isInitialized_ = init; }

  //! Set the vector associated with this operator (think nonlinear system)
  virtual void setSourceVector(const Teko::MultiVector& srcVec) { srcVector_ = srcVec; }

  //! Set the vector associated with this operator (think nonlinear system)
  virtual const Teko::MultiVector getSourceVector() const { return srcVector_; }

  //! Add a named inverse to the state object
  virtual void addInverse(const std::string& name, const Teko::InverseLinearOp& ilo) {
    inverses_[name] = ilo;
  }

  //! Get a named inverse from the state object
  virtual Teko::InverseLinearOp getInverse(const std::string& name) const {
    std::map<std::string, Teko::InverseLinearOp>::const_iterator itr;
    itr = inverses_.find(name);
    if (itr == inverses_.end()) return Teuchos::null;
    return itr->second;
  }

  //! Add a named operator to the state object
  virtual void addLinearOp(const std::string& name, const Teko::LinearOp& lo) {
    linearOps_[name] = lo;
  }

  //! Add a named operator to the state object
  virtual Teko::LinearOp getLinearOp(const std::string& name) { return linearOps_[name]; }

  //! Add a named operator to the state object
  virtual void addModifiableOp(const std::string& name, const Teko::ModifiableLinearOp& mlo) {
    modifiableOp_[name] = mlo;
  }

  //! Add a named operator to the state object
  virtual Teko::ModifiableLinearOp& getModifiableOp(const std::string& name) {
    std::map<std::string, Teko::ModifiableLinearOp>::iterator itr;
    itr = modifiableOp_.find(name);
    if (itr == modifiableOp_.end()) return modifiableOp_[name];
    return itr->second;
  }

  //! Merge internal storage of another PreconditionerState object into this one
  virtual void merge(const PreconditionerState& ps, int position = -1);

  //! Get the tag for this operator
  unsigned int getTag() const;

  //! Set the tag for this operator
  void setTag(unsigned int tag);

 protected:
  //! for ParameterListAcceptor
  Teuchos::RCP<Teuchos::ParameterList> paramList_;

  //! Store a source vector
  Teko::MultiVector srcVector_;

  //! Store a map of inverse linear operators
  std::map<std::string, Teko::InverseLinearOp> inverses_;
  std::map<std::string, Teko::ModifiableLinearOp> modifiableOp_;
  std::map<std::string, Teko::LinearOp> linearOps_;

  //! Stores the initialization state
  bool isInitialized_;

  unsigned int tag_;
};

}  // end namespace Teko

#endif
