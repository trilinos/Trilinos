#ifndef AMESOS_SCHWARZPRECONDITIONER_H
#define AMESOS_SCHWARZPRECONDITIONER_H

#include "Amesos_ConfigDefs.h"
#include "Epetra_Operator.h" 
class Epetra_RowMatrix;
class Epetra_MultiVector;
class Amesos_Container;
class Amesos_LocalRowMatrix;
class Amesos_Partitioner;
class Amesos_Graph;

#include "Amesos_InverseFactory.h"
#include "Teuchos_ParameterList.hpp"

class Amesos_BlockJacobiPreconditioner : public Epetra_Operator {

public:

  Amesos_BlockJacobiPreconditioner(char* Type,
			       Epetra_RowMatrix* Matrix,
			       Amesos_InverseFactory& Factory,
			       int NumLocalDomains_,
			       Teuchos::ParameterList& List);

  ~Amesos_BlockJacobiPreconditioner();

  bool IsSymmetric() const 
  {
    return(IsSymmetric_);
  }

  bool IsPreconditionerComputed() const
  {
    return(IsPreconditionerComputed_);
  }

  //@{ \name Atribute set methods.

  virtual int SetUseTranspose(bool UseTranspose)
  {
    AMESOS_CHK_ERR(-1); // FIXME: can I work with the transpose?
  }
  //@}

  //@{ \name Mathematical functions.

  virtual int Apply(const Epetra_MultiVector& X, 
		    Epetra_MultiVector& Y) const;

  virtual int ApplyInverse(const Epetra_MultiVector& X, 
			   Epetra_MultiVector& Y) const;

  virtual double NormInf() const
  {
    return(-1.0);
  }

  //@}

  //@{ \name Atribute access functions

  virtual char * Label() const;
 
  //! Returns the current UseTranspose setting.
  virtual bool UseTranspose() const
  {
    return(false);
  }

  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  virtual bool HasNormInf() const
  {
    return(false);
  }

  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  virtual const Epetra_Comm & Comm() const;

  //! Returns the Epetra_Map object associated with the domain of this operator.
  virtual const Epetra_Map & OperatorDomainMap() const;

  //! Returns the Epetra_Map object associated with the range of this operator.
  virtual const Epetra_Map & OperatorRangeMap() const;
  //@}

  int NumLocalDomains() const 
  {
    return(NumLocalDomains_);
  }

private:

  int ExtractSubmatrices();

  string Label_;
  bool IsSymmetric_;
  bool IsPreconditionerComputed_;
  string Type_;
  string SubmatricesType_;

  int NumLocalDomains_;
  Teuchos::ParameterList List_;

  Epetra_RowMatrix* RowMatrix_;
  Amesos_InverseFactory& Factory_;

  Amesos_Graph* Graph_;

  Amesos_Partitioner* Partitioner_;
  // FIXME: deleteme
  Amesos_LocalRowMatrix* LocalizedMatrix_;

  mutable vector<Amesos_Container*> Containers_;
  
}; // class Amesos_OverlappingAdditiveSchwarz

#endif // AMESOS_OVERLAPPINGSCHWARZPRECONDITIONER_H

