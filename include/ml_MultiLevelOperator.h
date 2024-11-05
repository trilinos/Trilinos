/*!
 * \file ml_MultiLevelOperator.h
 *
 * \brief Defines an ML preconditioner as a Epetra_Operator derived class.
 *
 * ML offers two preconditioners suitable for the solution of
 * Epetra_LinearProblem objects. This file define one the two, called
 * MultiLevelOperator (in the ML_Epetra namespace). This preconditioner is
 * simple wrapper of the ML_Solve() function, so that ML can be applied to
 * Epetra_MultiVector's.
 *
 * When you should use MultiLevelOperator:
 * - when your code already defines the required ML objects, with the optimal
 *   choice of parameters, and you want to use ML for Epetra_LinearProblem or
 *   AztecOO problems;
 *
 * When you should use MultiLevelPreconditioner:
 * - when you have an Epetra_RowMatrix, and you don't want to code the
 *   conversion to ML_Operator, the creation of the hierarchy and the
 *   aggregates, and/or you want to experiment various combinations of the
 *   parameters, simply changing some parameters in a Teuchos::ParameterList.
 *
 *
 * \date Last update to Doxygen: 22-Jul-04
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#ifndef ML_MULTILEVELOPERATOR_H
#define ML_MULTILEVELOPERATOR_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

class Epetra_MultiVector;
class Epetra_BlockMap;
class Epetra_Comm;

#include "ml_common.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_Object.h"
#include "ml_include.h"

//! ML_Epetra: default namespace for all Epetra interfaces.

namespace ML_Epetra
{

//! MultiLevelOperator: An implementation of the Epetra_Operator class.
/*! MultiLevelOperator class implements Epetra_Operator using a
    pre-constructed ML solver object. This allows ML to be used as
    preconditioner within an AztecOO solver object.
*/

class MultiLevelOperator: public virtual Epetra_Operator {

 public:

  //@{ \name Constructor.
    //! Uses an ML instance to implement the Epetra_Operator interface.
  /*! This is designed
      for using ML as a preconditioner within an AztecOO solver instance.
    \param ml_handle A fully-constructed ML object (In)
    \param myComm - Epetra communicator (In)
    \param DomainMap - Epetra domain map (In)
    \param RangeMap - Epetra range map (In)
  */
  MultiLevelOperator(ML * ml_handle, const Epetra_Comm & myComm,
                     const Epetra_Map & DomainMap,
                     const Epetra_Map & RangeMap);
  //@{ \name Destructor.
    //! Destructor
  ~MultiLevelOperator();
  //@}


  //@{ \name Attribute set methods.

    //! If set true, the multigrid hierarchy is destroyed when the Operator is destroyed.
    /*! This flag determines the ownership of the multigrid
        hierarchy. When set to true, this object owns the multigrid
        hierarchy and so it destroys it when freed. Otherwise, it is
        assumed that the multigrid hierarchy is owned by another object
        and so it is not freed.  By default, the multigrid hierarchy is
        not owned by this object.

    \param ownership (In) - If true, this object owns the corresponding
    multigrid hierarchy.

  */
  int SetOwnership(bool ownership){ ownership_ = ownership; return(-1);};
  //@}


  //@{ \name Attribute set methods.

    //! If set true, transpose of this operator will be applied.
    /*! This flag allows the transpose of the given operator to be used
        implicitly.  Setting this flag affects only the Apply() and
        ApplyInverse() methods.  If the implementation of this interface
        does not support transpose use, this method should return a
        value of -1.

	\param UseTranspose (In) - If true, multiply by the transpose of
	operator, otherwise just use operator.

    \warning - This method has no effect and returns -1 as error code.
  */
  int SetUseTranspose(bool UseTrans){ ML_avoid_unused_param((void *) UseTrans);
    return(-1);}
  //@}

  //@{ \name Mathematical functions.

  //! Returns the result of a Operator applied to a Epetra_MultiVector X in Y.
  /*!
    \param X (In) - A Epetra_MultiVector of dimension NumVectors to multiply with matrix.
    \param Y (Out) -A Epetra_MultiVector of dimension NumVectors containing result.
    \warning - This method has no effect and returns -1 as error code.
  */
  int Apply(const Epetra_MultiVector& /* X */, Epetra_MultiVector& /* Y */) const {
    return(-1);}

  //! Returns the result of a Operator inverse applied to an Epetra_MultiVector X in Y.
  /*!
    \param X (In) - A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Y (Out) -A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
  */
#ifdef WKC

#if WKC < 1
#error Blocking parameter (WKC) is not properly defined!
#endif


  int ApplyInverse_WKC(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
;
//  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y, int iBlockSize = WKC) const;
  int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;
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

  //@{ \name Attribute access functions

  //! Returns a character string describing the operator
  const char * Label() const{return(Label_);};

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

  //@{ Private data
  //! Pointer to the ML_Structure.
  ML * solver_;
  //! Label for this object.
  const char * Label_;

 private:
  //! Copy constructor (NOT DEFINED)
  MultiLevelOperator(const MultiLevelOperator& RHS) :
    Comm_(RHS.Comm()),
    DomainMap_(RHS.OperatorDomainMap()),
    RangeMap_(RHS.OperatorRangeMap())
  { }

  //! Operator= (NOT DEFINED)
  MultiLevelOperator& operator=(const MultiLevelOperator& /* RHS */)
  {
    return(*this);
  }

  //! Reference to Epetra communicator.
  const Epetra_Comm& Comm_;
  //! Reference to Domain Map.
  const Epetra_Map& DomainMap_;
  //! Reference to Range Map.
  const Epetra_Map& RangeMap_;
  bool  ownership_;
};

}

#endif /* ML_MULTILEVELOPERATOR_H */
