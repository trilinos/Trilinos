#ifndef __PB_StridedEpetraOperator_hpp__
#define __PB_StridedEpetraOperator_hpp__

// Epetra includes
#include "Epetra_Operator.h"

// Teuchos includes
#include "Teuchos_RCP.hpp"

#include "Thyra_LinearOpBase.hpp"

#include "Epetra/PB_EpetraOperatorWrapper.hpp"

// NOX includes
#include "NOX.H"
#include "NOX_Epetra.H"

namespace PB {
namespace Epetra {

class StridedEpetraOperator : public EpetraOperatorWrapper {
public:
   StridedEpetraOperator(int numVars,const Teuchos::RCP<Epetra_Operator> & content);
   StridedEpetraOperator(const std::vector<int> & vars,const Teuchos::RCP<Epetra_Operator> & content);

   virtual void SetContent(const std::vector<int> & vars,const Teuchos::RCP<Epetra_Operator> & content);

   virtual void RebuildOps()
   { BuildBlockedOperator(); }

   virtual const Teuchos::RCP<const Epetra_Operator> GetContent() const
   { return fullContent_; }

   virtual const Teuchos::RCP<Epetra_Operator> GetContent()
   { return fullContent_; }

   const Teuchos::RCP<const Epetra_Operator> GetBlock(int i,int j) const;

   // functions overloading Epetra_Operator
   ////////////////////////////////////////////////

   // destructor
   virtual ~StridedEpetraOperator() {}

   // attribute set methods
   
   // don't use transpose...ever!
   virtual int SetUseTranspose(bool useTranspose)
   { return -1; }

   // mathematical functions 
   virtual int Apply(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const;

   virtual int ApplyInverse(const Epetra_MultiVector &X, Epetra_MultiVector &Y) const
   { TEUCHOS_ASSERT(false); return -1; }

   virtual double NormInf() const
   { TEUCHOS_ASSERT(false); return 0.0; }

   // attribute access functions
   virtual bool UseTranspose() const { return false; }
   virtual bool HasNormInf() const { return false; }
   virtual const Epetra_Comm & Comm() const { return fullContent_->Comm(); }

protected:
   // gooey center of this shell
   Teuchos::RCP<Epetra_Operator> fullContent_;

   void BuildBlockedOperator();
};

} // end namespace Epetra
} // end namespace PB

#endif
