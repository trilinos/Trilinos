#ifndef __PB_BlockInvDiagonalStrategy_hpp__
#define __PB_BlockInvDiagonalStrategy_hpp__

#include <vector>

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_LinearOpBase.hpp"

namespace PB {

// this should be paired with a BlockJacobiPreconditionerFactory
// or BlockGSPreconditionerFactory.  The idea is that this object
// provides an approximate inverse operator for each of the diagonal
// blocks.  Then, the [whatever]PreconditionerFactory can easily
// construct an approximate inverse. The system under consideration
// is
// 
//    A = [ D0  U01 U02 ...]
//        [ L10  D1 U12 ...]
//        [ L20 L21  D2 ...]
//        [  .   .   .  ...]
//
// where inverses of D0,D1,D2...DN are needed.
//
class BlockInvDiagonalStrategy {
public:
   // returns an (approximate) inverse of the diagonal blocks of A
   virtual const std::vector<Teuchos::RCP<const Thyra::LinearOpBase<double> > > &
   getInvD(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const = 0;

   // return the number of diagonal blocks in the matrix
   virtual int numDiagonalBlocks() const = 0;
};

// this is a simple strategy for a [whatever]PreconditionerFactory
// it simply returns statically set RCP pointers to the passed in
// inv(D0) and inv(D1) operators. Not this will _not_ permit 
// efficient implementations when the preconditioner has to be rebuilt
// or reused often.
//
class StaticInvDiagStrategy : public BlockInvDiagonalStrategy {
public:
   StaticInvDiagStrategy(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invD0,
                          const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invD1)
   { invDiag_.push_back(invD0); invDiag_.push_back(invD1); }

   StaticInvDiagStrategy(const std::vector<Teuchos::RCP<const Thyra::LinearOpBase<double> > > & invD)
      : invDiag_(invD)
   { }

   // returns an (approximate) inverse of the diagonal blocks of A
   // where A is closely related to the original source for invD0 and invD1
   virtual const std::vector<Teuchos::RCP<const Thyra::LinearOpBase<double> > > &
   getInvD(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const
   { return invDiag_; }

   // return the number of diagonal blocks in the matrix
   virtual int numDiagonalBlocks() const { return invDiag_.size(); }

protected:
   // stored inverse operators
   std::vector<Teuchos::RCP<const Thyra::LinearOpBase<double> > > invDiag_;
};

} // end namespace PB

#endif
