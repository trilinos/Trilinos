/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifndef ML_MATRIX_FREE_PRECONDITIONER
#define ML_MATRIX_FREE_PRECONDITIONER

#if defined(ML_SHOW_DEPRECATED_WARNINGS)
#ifdef __GNUC__
#warning "The ML package is deprecated"
#endif
#endif

/*!
 * \file ml_MatrixFreePreconditioner.h
 *
 * \class MatrixFreePreconditioner
 *
 * \brief ML preconditioner for Matrix-Free operators
 *
 * \author Marzio Sala, ETHZ/D-INFK
 *
 * \date Last update to Doxygen: 01-Apr-06
 *
 */

#include "ml_include.h"
#if defined(HAVE_ML_EPETRA) && defined(HAVE_ML_TEUCHOS) && defined(HAVE_ML_EPETRAEXT) && defined(HAVE_ML_IFPACK)
#include <string>
#include "ml_epetra.h"
#include "Epetra_Time.h"
#include "Epetra_Operator.h"
#include "Epetra_Comm.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include <vector>
#include <map>
#include "ml_MultiLevelPreconditioner.h"

class Epetra_Map;
class Epetra_BlockMap;
class Epetra_CrsGraph;
class Epetra_Vector;
class Epetra_MultiVector;
class Epetra_RowMatrix;
class Epetra_FECrsMatrix;
class Ifpack_Chebyshev;

namespace ML_Epetra {

class MultiLevelPreconditioner;

/*!
\brief MatrixFreePreconditioner: a class to define preconditioners for Epetra_Operator's

This file requires ML to be configured with the following options:
- \c --enable-epetra
- \c --enable-epetraext
- \c --enable-teuchos

The following options are suggested:
- \c --enable-amesos
- \c --enable-ifpack

This class does not support Maxwell problems. It has been tested on symmetric problems; however it can in principle be used with non-symmetric problems as well.

\author Marzio Sala, ETHZ/D-INFK
*/
class MatrixFreePreconditioner : public Epetra_Operator
{
  public:
    //@{ \name Constructors and destructors.

    //! Constructor
    MatrixFreePreconditioner(const Epetra_Operator& Operator,
                             const Epetra_CrsGraph& Graph,
                             Epetra_MultiVector& NullSpace,
                             const Epetra_Vector& PointDiagonal,
                             Teuchos::ParameterList& List);

    //! destructor
    virtual ~MatrixFreePreconditioner();

    // @}
    // @{ \name Query methods.

    //! Sets the use of the transpose of the operator (NOT SUPPORTED).
    int SetUseTranspose(bool UseTranspose);

    //! Applies the operator to a std::vector (NOT SUPPORTED).
    int Apply(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Applies the preconditioner to std::vector \c X, returns the result in \c Y.
    int ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const;

    //! Returns the infinite norm of the operator (NOT SUPPORTED).
    double NormInf() const
    {
      return(-1.0);
    }

    //! Returns the label of \c this operator.
    const char * Label() const
    {
      return(Label_.c_str());
    }

    //! Returns \c true if the tranpose of the operator is considerd (NOT SUPPORTED).
    bool UseTranspose() const
    {
      return(false);
    }

    //! Returns \c false.
    bool HasNormInf() const
    {
      return(false);
    }

    //! Returns a reference to the communicator object.
    const Epetra_Comm& Comm() const
    {
      return(Comm_);
    }

    //! Returns the domain map of the operator.
    const Epetra_Map & OperatorDomainMap() const
    {
      return(Operator_.OperatorDomainMap());
    }

    //! Returns the range map of the operator.
    const Epetra_Map & OperatorRangeMap() const
    {
      return(Operator_.OperatorRangeMap());
    }

    //! Returns the coarser-level operator as an Epetra_RowMatrix.
    const Epetra_RowMatrix& C() const
    {
      if (!IsComputed())
        throw(-1); // prec not computed yet

      return(*C_);
    }

    const MultiLevelPreconditioner& MLP() const
    {
      if (!IsComputed())
        throw(-1); // prec not computed yet

      return(*MLP_);
    }

    //! Returns the restriction operator as an Epetra_CrsMatrix.
    const Epetra_CrsMatrix& R() const
    {
      return(*R_);
    }

    //! Returns the ML communicator of \c this object.
    ML_Comm* Comm_ML()
    {
      return(Comm_ML_);
    }

    //! Returns the PID of the calling processor.
    inline int MyPID() const
    {
      return(Comm().MyPID());
    }

    //! Returns the number of processors in the communicator.
    inline int NumProc() const
    {
      return(Comm().NumProc());
    }

    //! Returns \c true if the preconditioner has been successfully computed.
    bool IsComputed() const
    {
      return(IsComputed_);
    }

    //! Returns the total CPU time spent in \c this class.
    double TotalCPUTime() const;

    bool CheckSPD(const Epetra_Operator& Op,
                  const bool UseApply = true,
                  const int NumChecks = 1,
                  const int NumVectors = 1) const;

    // @}
    // @{ \name Construction methods.

    //! Performs coarsening for a given operator \c A.
    int Coarsen(ML_Operator* A, ML_Aggregate** aggr, ML_Operator** P,
                ML_Operator** R, ML_Operator** C, int NumPDEEqns = 1,
                int NullSpaceDim = 1, double* NullSpace = NULL);

    //! Probes for the block diagonal of the given operator.
    int GetBlockDiagonal(const Epetra_CrsGraph& Graph, std::string DiagonalColoringType);

  private:

    //! Computes the preconditioner.
    int Compute(const Epetra_CrsGraph& Graph, Epetra_MultiVector& NullSpace);

    // @}
    // @{ \name Basic smoothers

    //! Applies the pre-smoother (using zero starting solution).
    int ApplyPreSmoother(Epetra_MultiVector& X) const;

    //! Applies the post-smoother (using non-zero starting solution).
    int ApplyPostSmoother(Epetra_MultiVector& X, const Epetra_MultiVector& Y,
                          Epetra_MultiVector& tmp) const;

    //! Applies one sweep of Jacobi to std::vector \c X.
    int ApplyJacobi(Epetra_MultiVector& X, const double omega) const;

    //! Applies one sweep of Jacobi to std::vector \c X, using \c X as starting solution.
    int ApplyJacobi(Epetra_MultiVector& X, const Epetra_MultiVector& B,
                    const double omega, Epetra_MultiVector& tmp) const;

    //! Applies one sweep of block Jacobi to std::vector \c X.
    int ApplyBlockJacobi(Epetra_MultiVector& X, const double omega) const;

    //! Applies one sweep of block Jacobi to std::vector \c X, using \c X as starting solution.
    int ApplyBlockJacobi(Epetra_MultiVector& X, const Epetra_MultiVector& B,
                    const double omega, Epetra_MultiVector& tmp) const;

    int ApplyInvBlockDiag(const double alpha, Epetra_MultiVector& X,
                          const double gamma, const Epetra_MultiVector& B) const;

    // @}
    // @{ \name Timing

    inline void ResetStartTime() const
    {
      Time_->ResetStartTime();
    }

    inline void AddAndResetStartTime(const std::string& Label, const int print = false) const
    {
      TimeTable[Label] += Time_->ElapsedTime();
      Time_->ResetStartTime();
      if (print)
      {
        if (MyPID() == 0 && ML_Get_PrintLevel() > 5)
          std::cout << "Time for " << Label << " = " << TimeTable[Label] << " (s)" << std::endl;
      }
    }

    void PrintTimings() const
    {
      if (MyPID() == 0)
      {
        double Total = 0.0;
        std::cout << "Cumulative timing so far:" << std::endl;
        std::cout << "- for coarsening = " << TimeTable["coarsening"] << std::endl;
        std::cout << "- total time     = " << Total << std::endl;
      }
    }

    // @}
    // @{ \name Private data

    //! Toggles output level.
    bool verbose_;
    //! Communicator for ML.
    ML_Comm* Comm_ML_;
    //! Communicator object for Epetra.
    const Epetra_Comm& Comm_;

    //! Label of this object
    std::string Label_;
    //! Set to \c true if the preconditioner has been successfully computed.
    bool IsComputed_;
    //! Type of preconditioner (additive or hybrid)
    int PrecType_;
    //! Type of smoother (Jacobi, block Jacobi, or Chebyshev)
    int SmootherType_;
    //! Damping parameter for Jacobi.
    double omega_;
    //! List containing all the parameters
    Teuchos::ParameterList List_;

    //! Fine-level operator
    const Epetra_Operator& Operator_;
    //! Inverse of the point diagonal of the operator.
    Teuchos::RefCountPtr<Epetra_Vector> InvPointDiagonal_;
    //! Inverse of the diagonal of \c Operator_ as provided by the user.
    std::vector<double> InvBlockDiag_;

    //! Presmoother
    Teuchos::RefCountPtr<Ifpack_Chebyshev> PreSmoother_;
    //! Presmoother
    Teuchos::RefCountPtr<Ifpack_Chebyshev> PostSmoother_;
    //! Restriction from fine to coarse.
    Teuchos::RefCountPtr<Epetra_CrsMatrix> R_;
    //! Coarser-level operator as an Epetra_RowMatrix (wrapper for C_ML_).
    Teuchos::RefCountPtr<Epetra_RowMatrix> C_;
    //! Coarser-level operator as an ML_Operator.
    ML_Operator* C_ML_;
    //! Preconditioner that approximates the inverse of \c C_.
    Teuchos::RefCountPtr<MultiLevelPreconditioner> MLP_;

    //! Number of PDE equations
    int NumPDEEqns_;
    int NumMyBlockRows_;

    //! Time object.
    mutable Teuchos::RefCountPtr<Epetra_Time> Time_;
    mutable std::map<std::string, double> TimeTable;
    // @}

}; // class MatrixFreePreconditioner

} // namespace ML_Epetra

#endif // HAVE_ML_EPETRA
#endif // ML_MATRIX_FREE_PRECONDITIONER
