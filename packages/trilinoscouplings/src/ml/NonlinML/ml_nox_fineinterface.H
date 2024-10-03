// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
/*!
 * \file ml_nox_fineinterface.H
 *
 * \class Ml_Nox_Fineinterface
 *
 * \brief virtual class Ml_Nox_Fineinterface used by ML_Nox_Preconditioner
 *
 * \date Last update do Doxygen: 31-Mar-05
 *
 */
#include "TrilinosCouplings_config.h"

#if defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA) && defined(HAVE_ML_AZTECOO) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_IFPACK) && defined(HAVE_ML_AMESOS) && defined(HAVE_ML_EPETRAEXT)

#ifndef ML_NOX_FINEINTERFACE_H
#define ML_NOX_FINEINTERFACE_H

#include <iostream>
// ML headers
#include "ml_common.h"
#include "ml_include.h"
// Epetra headers
#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"
#include "Epetra_CrsGraph.h"
// NOX interfaces
#include "NOX_EpetraNew_Interface_Required.H"
#include "NOX_EpetraNew_Interface_Jacobian.H"
#include "NOX_EpetraNew_Interface_Preconditioner.H"
#include "NOX_Parameter_PrePostOperator.H"

//! ML_NOX: default namespace for all ML/NOX nonlinear methods
namespace ML_NOX
{
/*!
 
   \brief Ml_Nox_Fineinterface: a virtual class used by ML_Nox_Preconditioner

   ML_NOX::Ml_Nox_Fineinterface is a pure virtual class that serves as
   an interface between the ML_Nox_Preconditioner and any application.
   It implements all three NOX interfaces NOX::EpetraNew::Interface::Required,
   NOX::EpetraNew::Interface::Jacobian and NOX::EpetraNew::Interface::Preconditioner.
   Besides serving as an interface for the nonlinear iteration, it also serves as an 
   interface to do finite differencing of Jacobians and as the fine level interface all
   coarser level interfaces refer to. It therefore plays a major role in using
   the ML_Nox_Preconditioner.
   
    To make use of the full functionality of this class it 
    requires ML to be configured with the following options:
    - \c --with-ml_nox
    - \c --enable-epetra
    - \c --enable-epetraext
    - \c --enable-nox
    - \c --enable-nox-epetra
    - \c --enable-prerelease
    - \c --enable-aztecoo
    - \c --enable-amesos
    
   \note
    Currently, using this class makes it necessary to configure and make install
    of Trilinos WITHOUT the --with-ml_nox, and then repeating the configure
    and make install process WITH the --with-ml_nox option to be able to resolve the
    dependencies of this class on NOX header files

    \author Michael Gee, SNL 9214
*/  
class  Ml_Nox_Fineinterface : public virtual NOX::EpetraNew::Interface::Required,
                              public virtual NOX::EpetraNew::Interface::Jacobian,
                              public virtual NOX::EpetraNew::Interface::Preconditioner,
                              public virtual NOX::Parameter::PrePostOperator
                              
{
public:

  //! Constructor
  /*! Constructor
   */
  Ml_Nox_Fineinterface() 
  { isnewJacobian_=false; numJacobian_=0; printlevel_=10; t_ = 0.; ncalls_computeF_=0; return;}

  //! virtual Constructor
  /*! virtual Constructor
   */
  virtual ~Ml_Nox_Fineinterface() {};

  //! Compute and return F (derived from NOX::EpetraNew::Interface::Required)
  /*! Compute and return F (derived from NOX::EpetraNew::Interface::Required). This method is called
  by the nonlinear preconditioner. It supplies a solution vector x and needs a residual vector F computed
  by this funtion on output.
   \param x        (In)  : current solution vector
   \param F        (Out) : residual vector matching x
   \param fillFlag (In)  : flag indicating for what purpose this routine was called, see NOX documentation of
                           NOX::EpetraNew::Interface::Required
   */
  virtual bool computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag) = 0;

  //! Compute the current Jacobian (derived from NOX::EpetraNew::Interface::Jacobian)
  /*! Compute the current Jacobian (derived from NOX::EpetraNew::Interface::Jacobian). 
      This method is called by the nonlinear preconditioner. It supplies a solution vector x and wants
      a Jacobian computed matching this solution vector on output. 
      getJacobian will be used to obtain the Jacobian
      \sa getJacobian
   \param x        (In)  : current solution vector
   */
  virtual bool computeJacobian(const Epetra_Vector& x) = 0;
  
  // Compute Preconditioner (derived from NOX::EpetraNew::Interface::Preconditioner)
  //! Compute the current Preconditioner (derived from NOX::EpetraNew::Interface::Preconditioner)
  /*! FIXME: just found out that this interface is never used by ML_Nox_Preconditioner, 
     as ML_Nox_Preconditioner implements it's own NOX::EpetraNew::Interface::Preconditioner.
     delete this?
   */
  virtual bool computePreconditioner(const Epetra_Vector& x,
				     NOX::Parameter::List* precParams = 0) = 0;

  //! Return the Jacobian from the application.
  /*! Return the Jacobian from the application. The ML_Nox_Preconditioner uses this
      method to get the jacobian from the application. It will never be used in the matrixfree case.
      \sa ML_Nox_Preconditioner
   */
  virtual Epetra_CrsMatrix* getJacobian()  = 0;

  // return the graph from the problem
  // Return the Graph of the fine level problem from the application.
  //! Return the Graph of the fine level problem from the application.
  /*! Return the Graph of the fine level problem from the application. 
      The ML_Nox_Preconditioner uses this method to get the graph of
      the fine level problem from the application.
   */
  virtual const Epetra_CrsGraph* getGraph() = 0;
  virtual const Epetra_CrsGraph* getModifiedGraph() = 0;
  
  //! Return the current solution from the application.
  /*! Return the current solution from the application.
      The ML_Nox_Preconditioner uses this method to get the current solution vector of
      the fine level problem from the application.
   */
  virtual const Epetra_Vector* getSolution() = 0;
                               
  //! Return the Map of the current solution from the application.
  /*! Return the Map of the current solution from the application.
      The ML_Nox_Preconditioner uses this method to get the map of the current solution vector of
      the fine level problem from the application.
   */
  virtual const Epetra_Map& getMap() = 0;
                               
  //! Query whether the current Jacobian matches the current solution vector.
  /*! Query whether the current Jacobian matches the current solution vector.
   */
  virtual bool isnewJacobian() { return isnewJacobian_; }
  
  //! Query how often the Jacobian has been recomputed.
  /*! Query how often the Jacobian has been recomputed.
   */
  virtual int getnumJacobian() { return numJacobian_; }

  //! Compute and return the Nullspace of the problem from the application.
  /*! Compute and return the Nullspace of the problem from the application. This method is called
  by the nonlinear preconditioner to obtain the Nullspace of the problem from the application.
  If the application returns NULL on exit, the ML_Nox_Preconditioner will use a default Nullspace of
  correct size (given by dim_nullsp) and will assume constant number of dofs per node (constant block size).
  This will lead to errors in the aggregation process if the problem is NOT of constant block size!
   \param nummyrows     (In)  : dimension of the nullspace vectors
   \param numpde        (In) : number of PDEs (normally dofs per node)
   \param dim_nullsp    (In)  : dimension of the nullspace
   */
  virtual double* Get_Nullspace(const int nummyrows, const int numpde, 
                                const int dim_nullsp) = 0;
  
  
  // Apply constraints
  virtual void ApplyAllConstraints(Epetra_Vector& gradient, int level) 
  { return;}
  
  // get block information
  //! Compute and return the block node information from the application.
  /*! Compute and return the block node information from the application. This method is called
  by the nonlinear preconditioner to obtain information about non-constant block sizes. It is used
  only in the case of the VBMETIS aggregation scheme. If on exit blocks and block_pde are NULL, the
  ML_Nox_Preconditioner will assume a constant vblock size and will use METIS instead of the VBMETIS
  aggregation scheme.
   \param nblocks     (Out)  : on exit, number of variable sized nodal blocks
   \param blocks      (Out)  : on exit, allocated vector of length NumMyRows containing the block indize
                               of each row in C-style numbering (starting with 0)
   \param block_pde   (Out)  : on exit, allocated vector of length NumMyRows containing the number of
                               PDE a certain row belongs to.
   */
  virtual bool getBlockInfo(int *nblocks, int** blocks, int **block_pde) = 0;

  //! Query the summed time spent in this interface.
  /*! Query the summed time spent in this interface. This is for performance analysis.
   */
  double getsumtime() {return t_;}
  
  //! Reset the summed time spent in this interface to a certain value.
  /*! Reset the summed time spent in this interface to a certain value.
      This is for performance analysis.
   */
  void resetsumtime() { t_ = 0.; return; }
  
  //! Query the number of calls to the computeF method of this class.
  /*! Query the number of calls to the computeF method of this class. This is for performance analysis.
   */
  int getnumcallscomputeF() { return ncalls_computeF_; }
  
  //! Reset the number of calls to the computeF method of this class.
  /*! Reset the number of calls to the computeF method of this class. This is for performance analysis.
   */
  bool setnumcallscomputeF(int ncalls) { ncalls_computeF_=ncalls; return true; }
  
  // get the output level 
  inline int OutLevel() { return printlevel_; }

  // output method to look at the iterate
#if 0
  virtual bool PrintSol(int id) = 0;
#endif  

public:

  //! Flag to indicate, that the Jacobian has been evaluated by the application.
  /*! Flag to indicate, that the Jacobian has been evaluated by the application.
      the ML_NOX::Ml_Nox_Fineinterface should set this to true, whenever
      the evaluated Jacobian matches the F in the computeF.
      It should be set to false if there is either no Jacobian evaluated or
      it does not match the F from the computeF function
   */
  bool isnewJacobian_;
  
  //! count number of evaluated Jacobians
  /*! The user's ML_NOX::Ml_Nox_Fineinterface derived class should use this to
      count the number of Jacobian evaluations for statistical and debug-purpose
  */
  int numJacobian_;
  
  //! Printlevel
  /*! The user's ML_NOX::Ml_Nox_Fineinterface derived class should set a
  printlevel value between 0 and 10
  */  
  int printlevel_;

  //! t_ should be used for time measurement of time spent in this class
  /*! t_ should be used for time measurement. The methods getsumtime and resetsumtime will be used
      to query this value for performance measurement
  */  
  double t_;
  
  //! ncalls_computeF_ is used to count the number of calls to computeF
  /*! ncalls_computeF_ is used to count the number of calls to computeF
  */  
  int ncalls_computeF_;
  
private:  
  // don't want copy-ctor and = operator
  Ml_Nox_Fineinterface(ML_NOX::Ml_Nox_Fineinterface& old);
  //Ml_Nox_Fineinterface operator = (const ML_NOX::Ml_Nox_Fineinterface& old);
  
};
} // namespace ML_NOX


#endif // ML_NOX_FINEINTERFACE_H
#endif // defined(HAVE_ML_NOX) && defined(HAVE_ML_EPETRA)
