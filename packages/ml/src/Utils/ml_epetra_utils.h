#ifndef _ML_EPETRA_UTILS_H_
#define _ML_EPETRA_UTILS_H_

class Epetra_Comm;
class Epetra_BlockMap;
class Epetra_MultiVector;
class Epetra_RowMatrix;
class Epetra_Map;
class Epetra_Vector;
class Epetra_Import;
class Epetra_Object;
class Epetra_CrsMatrix;
class Epetra_RowMatrix;
class Epetra_LinearProblem;

#include "ml_common.h"

#ifdef ML_MPI
#ifndef EPETRA_MPI
#define EPETRA_MPI
#endif
#include "mpi.h"
#endif
#include "ml_include.h"

// External prototypes
int Epetra_ML_matvec(void *data, int in, double *p, int out,
                 double *ap);

#ifdef WKC
int Epetra_ML_matvec_WKC(void *data, int in, double *p, int out,
                 double *ap);
#endif 


int Epetra_ML_getrow(void *data, int N_requested_rows,
                 int requested_rows[], int allocated_space, int columns[],
                 double values[], int row_lengths[]);

int Epetra_ML_comm_wrapper(double vec[], void *data);

int EpetraMatrix2MLMatrix(ML *ml_handle, int level,
                                Epetra_RowMatrix * Amat);

int Epetra2MLMatrix(Epetra_RowMatrix * A, ML_Operator *Result);

Epetra_CrsMatrix *Epetra_MatrixMult(Epetra_RowMatrix *B, Epetra_RowMatrix *Bt);
Epetra_CrsMatrix *Epetra_MatrixAdd(Epetra_RowMatrix *B, Epetra_RowMatrix *Bt, double scalar);
int ML_Epetra_CRSinsert(ML_Operator *, int, int *, double *, int);

int ML_Operator2EpetraCrsMatrix(ML_Operator *Ke, Epetra_CrsMatrix * &
				CrsMatrix, int & MaxNumNonzeros,
				bool CheckNonzeroRow, double &);

/* This Proto-type is in ml_rap.h. This prevents ml_rap.c and ml_matmat_mult.c     */
/* from including ml_epetra_utils.h which would require the C++ compiler for these */
/* files.             
extern int  ML_back_to_epetraCrs(ML_Operator *Mat1Mat2,  ML_Operator *Result, 
				 ML_Operator *Mat1,  ML_Operator *Mat2); 
*/


#define ML_HAVE_TEUCHOS
#ifdef ML_HAVE_TEUCHOS

#include "Teuchos_ParameterList.hpp"
#include "ml_amesos_wrap.h"
#include "ml_ifpack_wrap.h"
#include "ml_agg_METIS.h"

#include "Epetra_Operator.h"

using namespace Teuchos;

class Epetra_ML_Preconditioner: public virtual Epetra_Operator {
      
 public:

  Epetra_ML_Preconditioner( const Epetra_RowMatrix & RowMatrix,
			    ParameterList & List);

  ~Epetra_ML_Preconditioner() {
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

#endif /* _ML_EPETRA_UTILS_H_ */
