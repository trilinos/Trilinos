
/* Copyright (2001) Sandia Corportation. Under the terms of Contract 
 * DE-AC04-94AL85000, there is a non-exclusive license for use of this 
 * work by or on behalf of the U.S. Government.  Export of this program
 * may require a license from the United States Government. */


/* NOTICE:  The United States Government is granted for itself and others
 * acting on its behalf a paid-up, nonexclusive, irrevocable worldwide
 * license in ths data to reproduce, prepare derivative works, and
 * perform publicly and display publicly.  Beginning five (5) years from
 * July 25, 2001, the United States Government is granted for itself and
 * others acting on its behalf a paid-up, nonexclusive, irrevocable
 * worldwide license in this data to reproduce, prepare derivative works,
 * distribute copies to the public, perform publicly and display
 * publicly, and to permit others to do so.
 * 
 * NEITHER THE UNITED STATES GOVERNMENT, NOR THE UNITED STATES DEPARTMENT
 * OF ENERGY, NOR SANDIA CORPORATION, NOR ANY OF THEIR EMPLOYEES, MAKES
 * ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LEGAL LIABILITY OR
 * RESPONSIBILITY FOR THE ACCURACY, COMPLETENESS, OR USEFULNESS OF ANY
 * INFORMATION, APPARATUS, PRODUCT, OR PROCESS DISCLOSED, OR REPRESENTS
 * THAT ITS USE WOULD NOT INFRINGE PRIVATELY OWNED RIGHTS. */

#ifndef _ML_EPETRA_OPERATOR_H_
#define _ML_EPETRA_OPERATOR_H_

class Epetra_MultiVector;
class Epetra_BlockMap;
class Epetra_Comm;

#include "ml_common.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Object.h"
#include "ml_include.h"

//! Operator: An implementation of the Epetra_Operator class.
/*! Operator class implements Epetra_Operator using a
    pre-constructed ML solver object.
    Once constructed, an Operator can be used as a preconditioner
    within an AztecOO solver object.
*/    

namespace ML_Epetra 
{
  
class MultiLevelOperator: public virtual Epetra_Operator {
      
 public:

  //@{ \name Constructor.
    //! Uses an ML instance to implement the Epetra_Operator interface.
  /*! This is designed
      for using ML as a preconditioner within an AztecOO solver instance.
    \param In - A fully-constructed ML object.
    \param In - An Epetra communicator
    \param In - An Epetra domain map
    \param In - An Epetra range map
  */
  MultiLevelOperator(ML * ml_handle, const Epetra_Comm & myComm,
                     const Epetra_Map & DomainMap,
                     const Epetra_Map & RangeMap);
  //@{ \name Destructor.
    //! Destructor
  ~MultiLevelOperator();
  //@}

  
  //@{ \name Atribute set methods.

    //! If set true, the multigrid hierarchy is destroyed when the Operator is destroyed.
    /*! This flag determines the ownership of the multigrid hierarchy. When set to true, this
        object owns the multigrid hierarchy and so it destroys it when freed. Otherwise, it is
        assumed that the multigrid hierarchy is owned by another object and so it is not freed.
        By default, the multigrid hierarchy is not owned by this object.
      
    \param In
	   ownership - If true, this object owns the corresponding multigrid hierarchy. 

  */
  int SetOwnership(bool ownership){ ownership_ = ownership; return(-1);};
  //@}

  
  //@{ \name Atribute set methods.

    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used implicitly.  Setting this flag
        affects only the Apply() and ApplyInverse() methods.  If the implementation of this interface 
	does not support transpose use, this method should return a value of -1.
      
    \param In
	   UseTranspose - If true, multiply by the transpose of operator, otherwise just use operator.

    \warning - This method has no effect and returns -1 as error code.
  */
  int SetUseTranspose(bool UseTranspose){ ML_avoid_unused_param((void *) UseTranspose); 
    return(-1);}
  //@}
  
  //@{ \name Mathematical functions.

    //! Returns the result of a Operator applied to a Epetra_MultiVector X in Y.
    /*! 
    \param In
	   X - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Out
	   Y -A Epetra_MultiVector of dimension NumVectors containing result.

    \warning - This method has no effect and returns -1 as error code.
  */
  int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const {
    return(-1);}

  //! Returns the result of a Operator inverse applied to an Epetra_MultiVector X in Y.
  /*! 
    \param In
    X - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
    Y -A Epetra_MultiVector of dimension NumVectors containing result.
    
    \return Integer error code, set to 0 if successful.
  */
#ifdef WKC

#if WKC < 1
#error WKC is not properly defined!
#endif

  int ApplyInverse_WKC(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
;
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y, int iBloc
kSize = WKC) const;
#else
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
#endif

  //int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
  
  //! Returns the infinity norm of the global matrix.
  /* Returns the quantity \f$ \| A \|_\infty\f$ such that
     \f[\| A \|_\infty = \max_{1\lei\lem} \sum_{j=1}^n |a_{ij}| \f].
     
     \warning This method must not be called unless HasNormInf() returns true.
  */ 
  double NormInf() const {return(0.0);};
  //@}
  
  //@{ \name Atribute access functions

  //! Returns a character string describing the operator
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
  
 protected:

  ML * solver_;
  char * Label_;

 private:
  const Epetra_Map & DomainMap_;
  const Epetra_Map & RangeMap_;
  const Epetra_Comm & Comm_;
  bool  ownership_;
};
 
}

#endif /* _ML_EPETRA_OPERATOR_H_ */
