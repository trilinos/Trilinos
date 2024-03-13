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

#ifndef __Teko_SolveInverseFactory_hpp__
#define __Teko_SolveInverseFactory_hpp__

#include "Teko_InverseFactory.hpp"

namespace Teko {

class SolveInverseFactory : public InverseFactory {
 public:
  //! \name Constructors
  //@{

  /** \brief Constructor that takes a Thyra solve factory and
   *        makes it look like an InverseFactory
   *
   * Constructor that takes a Thyra solve factory and
   * makes it look like an InverseFactory.
   *
   * \param[in] lowsFactory Thyra LineaerOpWithSolveFactoryBase used for building
   *                        the inverse.
   */
  SolveInverseFactory(
      const Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >& lowsFactory);

  //! Copy constructor
  SolveInverseFactory(const SolveInverseFactory& siFactory);
  //@}

  virtual ~SolveInverseFactory() {}

  /** \brief Build an inverse operator
   *
   * Build the inverse operator using this factory.
   *
   * \param[in] linearOp Linear operator needing to be inverted.
   *
   * \returns New linear operator that functions as the inverse
   *          of <code>linearOp</code>.
   */
  virtual InverseLinearOp buildInverse(const LinearOp& linearOp) const;

  /** \brief Build a preconditioned inverse operator
   *
   * Build the inverse operator using this factory and a user specified
   * preconditioning operator. The default behavior is to call buildInverse
   * ignoring the preconditioner.
   *
   * \param[in] linearOp Linear operator needing to be inverted.
   * \param[in] precOp Preconditioning operator
   *
   * \returns New linear operator that functions as the inverse
   *          of <code>linearOp</code>.
   */
  virtual InverseLinearOp buildInverse(const LinearOp& linearOp, const LinearOp& precOp) const;

  /** \brief Pass in an already constructed inverse operator. Update
   *        the inverse operator based on the new source operator.
   *
   * Pass in an already constructed inverse operator. Update
   * the inverse operator based on the new source operator.
   *
   * \param[in]     source Source operator to be inverted.
   * \param[in,out] dest   Pre constructed inverse operator to be
   *                        rebuilt using the <code>source</code>
   *                        object.
   */
  virtual void rebuildInverse(const LinearOp& source, InverseLinearOp& dest) const;

  /** \brief Pass in an already constructed inverse operator. Update
   *        the inverse operator based on the new source operator.
   *
   * Pass in an already constructed inverse operator. Update
   * the inverse operator based on the new source operator.
   *
   * \param[in]     source Source operator to be inverted.
   * \param[in]     precOp Preconditioning operator
   * \param[in,out] dest   Pre constructed inverse operator to be
   *                        rebuilt using the <code>source</code>
   *                        object.
   */
  virtual void rebuildInverse(const LinearOp& source, const LinearOp& precOp,
                              InverseLinearOp& dest) const;

  /** \brief A function that permits inspection of the parameters used to create
   *        this object.
   *
   * A function that permits inspection of the parameters used to create this
   * object. Useful for determining defaults and settings used.
   *
   * \returns A list used to parameterize this object.
   */
  virtual Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;

  //! Accessor primarily for testing purposes
  Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<double> > getLowsFactory() const {
    return lowsFactory_;
  }

  /** Return a string that describes this factory */
  virtual std::string toString() const { return lowsFactory_->description(); }

 protected:
  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > lowsFactory_;

 private:
  // hide me!
  SolveInverseFactory();
};

}  // end namespace Teko

#endif
