#ifndef _ML_EPETRA_OPERATOR_H_
#define _ML_EPETRA_OPERATOR_H_

class Epetra_Map;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_RowMatrix;
class Epetra_Operator;
class Epetra_Comm;

#include "ml_config.h"

#ifdef HAVE_ML_TEUCHOS

#include "Teuchos_ParameterList.hpp"

using namespace Teuchos;

class Epetra_ML_Preconditioner: public virtual Epetra_Operator {
      
 public:

  Epetra_ML_Preconditioner( const Epetra_RowMatrix & RowMatrix,
			    ParameterList & List);

  ~Epetra_ML_Preconditioner() {
    Destroy_ML_Preconditioner(); 
 }
  
  int SetOwnership(bool ownership){ ownership_ = ownership; return(-1);};

  int SetUseTranspose(bool UseTranspose){return(-1);}

  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(-1);}

  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

  double NormInf() const {return(0.0);};

  char * Label() const{return(Label_);};
  
  //! Returns the current UseTranspose setting.
  bool UseTranspose() const {return(false);};
  
  //! Returns true if the \e this object can provide an approximate Inf-norm, false otherwise.
  bool HasNormInf() const{return(false);};
  
  //! Returns a pointer to the Epetra_Comm communicator associated with this operator.
  const Epetra_Comm & Comm() const{return(Comm_);};
  
  //! Returns the Epetra_Map object associated with the domain of this operator.
  const Epetra_Map & OperatorDomainMap() const {return(DomainMap_);};
  //! Returns the Epetra_Map object associated with the range of this operator.
  const Epetra_Map & OperatorRangeMap() const {return(RangeMap_);};
  //@}
  
  void Destroy_ML_Preconditioner();

protected:

  ML * ml_;
  ML_Aggregate *agg_;
  
  char * Label_;

 private:

  int CreateLabel();
  
  int NumLevels_;
  const Epetra_Map & DomainMap_;
  const Epetra_Map & RangeMap_;
  const Epetra_Comm & Comm_;
  bool  ownership_;
  int   ProcConfig_[AZ_PROC_SIZE];
  int   SmootherOptions_[20][AZ_OPTIONS_SIZE];
  double SmootherParams_[20][AZ_PARAMS_SIZE];
  double SmootherStatus_[AZ_STATUS_SIZE];
  ParameterList & List_;
  int MaxLevels_;
  
};

#endif

#endif /* _ML_EPETRA_OPERATOR_H_ */
