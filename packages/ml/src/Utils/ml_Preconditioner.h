/*!
 * \file ml_Preconditioner.h
 *
 * \class ML_Preconditioner
 *
 * \brief (Mostly) abstract base class wrapper for Epetra_Operator-based ML
 * preconditioners.  Implements a few Teuchos-related query functions.
 *
 * \date Last update to Doxygen: 13-Nov-06
 *
 */
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */

#ifndef ML_PRECONDITIONER_H
#define ML_PRECONDITIONER_H

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS)
#include "Epetra_Comm.h"
#include "Epetra_Map.h"
#include "Epetra_Operator.h"
#include "Teuchos_ParameterList.hpp"
#include "ml_utils.h"

namespace ML_Epetra
{

  class ML_Preconditioner: public virtual Epetra_Operator
  {
  public:
    //! @name Constructor
    //@{
    //! Constructor
    ML_Preconditioner():Label_(0),IsComputePreconditionerOK_(0){};
    //@}

    //! @name Destructor
    //@{
    //! Destructor
    virtual ~ML_Preconditioner() {if(Label_) delete [] Label_;};
    //@}

    //@{ \name Query & Set functions
    //! Prints label associated to this object.
    virtual const char* Label() const{return(Label_);};

    //! Prints unused parameters in the input ParameterList on standard output.
    virtual void PrintUnused() const
    {
      List_.unused(std::cout);
    }

    //! Prints unused parameters in the input ParameterList on the specified stream.
    virtual void PrintUnused(std::ostream & os) const
    {
      List_.unused(os);
    }

    //! Prints unused parameters in the input ParameterList to std::cout on proc \c MyPID.
    /*! Mispelled parameters are simply ignored. Therefore, it is often the best
     * choice to print out the parameters that have not been used in the
     * construction phase.
     * - \param MyPID (In) : ID of process that should print the unused parameters.
     */
    virtual void PrintUnused(const int MyPID) const
    {
      if( Comm().MyPID() == MyPID ) {
        ML_print_line("-",78);
        std::cout << "Unused parameters:" << std::endl;
        PrintUnused();
        ML_print_line("-",78);
      }
    }

    //! Gets a reference to the internally stored parameters' list.
    virtual Teuchos::ParameterList& GetList()
    {
      return List_;
    }

    //! Prints on \c std::cout the values of the internally stored parameter list for processor \c MyPID
    virtual void PrintList(int MyPID)
    {
      if( Comm().MyPID() == MyPID ) {
        ML_print_line("-",78);
        std::cout << List_;
        ML_print_line("-",78);
      }
    }

    //! Copies \c List into the internally stored parameter list object.
    virtual int SetParameterList(const Teuchos::ParameterList& List)
    {
      if( IsComputePreconditionerOK_ == true ) DestroyPreconditioner();
      List_ = List;
      return 0;
    }

    //@}


    //@{ \name Attribute access functions

    //! Computes the multilevel hierarchy.
    /*! Computes the multilevel hierarchy. This function retrives the user's defines parameters (as
      specified in the input ParameterList), or takes default values otherwise, and creates the ML
      objects for aggregation and hierarchy. Allocated data can be freed used DestroyPreconditioner(),
      or by the destructor,

      In a Newton-type procedure, several linear systems have to be solved, Often, these systems
      are not too different. In this case, it might be convenient to keep the already
      computed preconditioner (with hierarchy, coarse solver, smoothers), and use it to
      precondition the next linear system. ML offers a way to determine whether the
      already available preconditioner is "good enough" for the next linear system.
      The user should proceed as follows:
      - define \c "reuse: enable" == \c true
      - solve the first linear system. ML tries to estimate the rate of convergence, and record it;
      - change the values of the linear system matrix (but NOT its structure)
      - compute the new preconditioner as \c ComputePreconditioner(true)
      It is supposed that the pointer to the Epetra_RowMatrix remains constant. Currently,
      it is not possible to modify this pointer (other than creating a new preconditioner)  */

    //! Computes the preconditioner
    virtual int ComputePreconditioner(const bool CheckFiltering)=0;

    //! Recomputes the preconditioner
    virtual int ReComputePreconditioner()=0;

    //! Print the individual operators in the multigrid hierarchy.
    virtual void Print(int whichHierarchy=-2)=0;

    //! Queries whether multilevel hierarchy has been computed or not.
    virtual int IsPreconditionerComputed() const
    {
      return(IsComputePreconditionerOK_);
    }

    //! Destroys all structures allocated in \c ComputePreconditioner() if the preconditioner has been computed.
    virtual int DestroyPreconditioner()=0;

    //! Return operator complexity and #nonzeros in fine grid matrix.
    virtual void Complexities(double &complexity, double &fineNnz)=0;

    //@}


  protected:
    //@{ \name Internal data

    //! Label for this object
    char* Label_;

    //! Specifies whether a hierarchy already exists or not.
    bool IsComputePreconditionerOK_;

    //! List containing all input parameters.
    Teuchos::ParameterList List_;
    //@}
  };//ML_Preconditioner

}//end namespace ML_Epetra


#endif

#endif
