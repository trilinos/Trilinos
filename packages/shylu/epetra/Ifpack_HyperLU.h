/** \file Ifpack_HyperLU.h

    \brief Use HyperLU as a preconditioner within IFPACK.

    \author Siva Rajamanickam

*/

#ifndef IFPACK_HYPERLU_H
#define IFPACK_HYPERLU_H

#include <assert.h>
#include <iostream>
#include <sstream>

/* This define will make S block diagonal, To make this happen even some
   zero columns/rows of C and R will be stored.
   */
#define BLOCK_DIAGONAL_Si

#include "Isorropia_config.h" // Just for HAVE_MPI

// Epetra includes
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_Time.h" 
#include "Epetra_CrsMatrix.h" 
#include "Epetra_Map.h" 
#include "Epetra_MultiVector.h" 
#include "Epetra_LinearProblem.h"
#include "Epetra_Import.h" 
#include "Epetra_Export.h" 

// Teuchos includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_LAPACK.hpp"

// AztecOO includes
#include "AztecOO.h"


#include "hyperlu_util.h"

// Ifpack includes
#include "Ifpack_Preconditioner.h"

//Isorropia includes
#include "Isorropia_Epetra.hpp"
#include "Isorropia_EpetraRedistributor.hpp"
#include "Isorropia_EpetraPartitioner.hpp"

// EpetraExt includes
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_RowMatrixOut.h"

// Amesos includes
#include "Amesos.h"
#include "Amesos_BaseSolver.h"

using namespace std;

class Ifpack_HyperLU: public Ifpack_Preconditioner
{
    public:
    // @{ Constructors and destructors.
    //! Constructor
    Ifpack_HyperLU(Epetra_CrsMatrix* A);

    //! Destructor
    ~Ifpack_HyperLU()
    {
        Destroy();
    }

    // @}
    // @{ Construction methods

    //! Initialize the preconditioner, does not touch matrix values.
    int Initialize();

    //! Returns \c true if the preconditioner has been successfully initialized.
    bool IsInitialized() const
    {
        return(IsInitialized_);
    }

    //! Compute ILU factors L and U using the specified parameters.
    int Compute();

    //! If factor is completed, this query returns true, otherwise it returns false.
    bool IsComputed() const 
    {
        return(IsComputed_);
    }

    //! Set parameters using a Teuchos::ParameterList object.
    /* This method is only available if the Teuchos package is enabled.
     This method recognizes four parameter names: relax_value,
     absolute_threshold, relative_threshold and overlap_mode. These names are
     case insensitive, and in each case except overlap_mode, the ParameterEntry
     must have type double. For overlap_mode, the ParameterEntry must have
     type Epetra_CombineMode.
    */
    int SetParameters(Teuchos::ParameterList& parameterlist);

    int SetUseTranspose(bool UseTranspose_in) {
        UseTranspose_ = UseTranspose_in;
     return(0);
     };

    //! Returns the current UseTranspose setting.
    bool UseTranspose() const {return(UseTranspose_);};




    // @}

    // @{ Mathematical functions.
    // Applies the matrix to X, returns the result in Y.
    int Apply(const Epetra_MultiVector& X, 
           Epetra_MultiVector& Y) const
    {
    return(Multiply(false,X,Y));
    }

    int Multiply(bool Trans, const Epetra_MultiVector& X, 
           Epetra_MultiVector& Y) const{return A_->Multiply(Trans,X,Y);}

    //! Returns the result of a Epetra_Operator inverse applied to an Epetra_MultiVector X in Y.
    /*!
    \param 
       X - (In) A Epetra_MultiVector of dimension NumVectors to solve for.
    \param Out
       Y - (Out) A Epetra_MultiVector of dimension NumVectors containing result.

    \return Integer error code, set to 0 if successful.
    */
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Computes the estimated condition number and returns the value.
    double Condest(const Ifpack_CondestType CT = Ifpack_Cheap, 
                 const int MaxIters = 1550,
                 const double Tol = 1e-9,
         Epetra_RowMatrix* Matrix_in = 0);

    //! Returns the computed estimated condition number, or -1.0 if not computed.
int JustTryIt() ;

    double Condest() const
    {
        return(Condest_);
    }

    // @}
    // @{ Query methods

    //! Returns a character string describing the operator
    const char* Label() const {return(Label_);}

    //! Sets label for \c this object.
    int SetLabel(const char* Label_in)
    {
        strcpy(Label_,Label_in);
        return(0);
    }


    //! Returns 0.0 because this class cannot compute Inf-norm.
    double NormInf() const {return(0.0);};

    //! Returns false because this class cannot compute an Inf-norm.
    bool HasNormInf() const {return(false);};

    //! Returns the Epetra_Map object associated with the domain of this operator.
    const Epetra_Map & OperatorDomainMap() const {return(A_->OperatorDomainMap());};

    //! Returns the Epetra_Map object associated with the range of this operator.
    const Epetra_Map & OperatorRangeMap() const{return(A_->OperatorRangeMap());};

    //! Returns the Epetra_BlockMap object associated with the range of this matrix operator.
    const Epetra_Comm & Comm() const{return(A_->Comm());};

    //! Returns a reference to the matrix to be preconditioned.
    const Epetra_RowMatrix& Matrix() const
    { 
        return(*A_);
    }

    //! Prints on stream basic information about \c this object.
    virtual ostream& Print(ostream& os) const;

    //! Returns the number of calls to Initialize().
    virtual int NumInitialize() const
    {
        return(NumInitialize_);
    }

    //! Returns the number of calls to Compute().
    virtual int NumCompute() const
    {
        return(NumCompute_);
    }

    //! Returns the number of calls to ApplyInverse().
    virtual int NumApplyInverse() const
    {
        return(NumApplyInverse_);
    }

    //! Returns the time spent in Initialize().
    virtual double InitializeTime() const
    {
        return(InitializeTime_);
    }

    //! Returns the time spent in Compute().
    virtual double ComputeTime() const
    {
        return(ComputeTime_);
    }

    //! Returns the time spent in ApplyInverse().
    virtual double ApplyInverseTime() const
    {
        return(ApplyInverseTime_);
    }

    //! Returns the number of flops in the initialization phase.
    virtual double InitializeFlops() const
    {
        return(0.0);
    }

    virtual double ComputeFlops() const
    {
        return(ComputeFlops_);
    }

    virtual double ApplyInverseFlops() const
    {
        return(ApplyInverseFlops_);
    }

    private:

    // @}
    // @{ Private methods

    //! Copy constructor (should never be used)
    Ifpack_HyperLU(const Ifpack_HyperLU& RHS) :
    Time_(RHS.Comm())
    {}

    //! operator= (should never be used)
    Ifpack_HyperLU& operator=(const Ifpack_HyperLU& RHS)
    {
        return(*this);
    }

    //! Destroys all internal data
    void Destroy();

    //! Returns the number of global matrix rows.
    int NumGlobalRows() const {return(A_->NumGlobalRows());};

    //! Returns the number of global matrix columns.
    int NumGlobalCols() const {return(A_->NumGlobalCols());};

    //! Returns the number of local matrix rows.
    int NumMyRows() const {return(A_->NumMyRows());};

    //! Returns the number of local matrix columns.
    int NumMyCols() const {return(A_->NumMyCols());};

    // @}
    // @{ Internal data

    //! Pointer to the Epetra_RowMatrix to factorize
    Epetra_CrsMatrix *A_;
    Epetra_CrsMatrix *C_;
    Teuchos::ParameterList List_;
    Isorropia::Epetra::Partitioner *partitioner_;
    Isorropia::Epetra::Redistributor *rd_;
    Epetra_LinearProblem *LP_;
    Amesos_BaseSolver *Solver_;
    int Dnr_;
    int Snr_;
    int *DRowElems_;
    int *SRowElems_;
    Teuchos::RCP<Epetra_CrsMatrix> Sbar_;
    mutable AztecOO *solver_; // Ugh !!! Mutable ! Ugh! Workaround AztecOO bug
    Epetra_MultiVector *Xs_;
    Epetra_MultiVector *Bs_;
    Epetra_Operator *precop_ ;
    string libName_;

    bool UseTranspose_;

    double Condest_;

    bool IsParallel_;
    //! If \c true, the preconditioner has been successfully initialized.
    bool IsInitialized_;
    //! If \c true, the preconditioner has been successfully computed.
    bool IsComputed_;
    //! Label of \c this object.
    char Label_[160];
    //! Contains the number of successful calls to Initialize().
    int NumInitialize_;
    //! Contains the number of successful call to Compute().
    int NumCompute_;

    //! Contains the number of successful call to ApplyInverse().
    mutable int NumApplyInverse_;
    //! Contains the time for all successful calls to Initialize().
    double InitializeTime_;
    //! Contains the time for all successful calls to Compute().
    double ComputeTime_;
    //! Contains the time for all successful calls to ApplyInverse().
    mutable double ApplyInverseTime_;
    //! Contains the number of flops for Compute().
    double ComputeFlops_;
    //! Contain sthe number of flops for ApplyInverse().
    mutable double ApplyInverseFlops_;
    //! Used for timing issues
    mutable Epetra_Time Time_;
    // @}

};

#endif /* IFPACK_HYPERLU_H */
