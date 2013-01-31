/*@HEADER
// ***********************************************************************
//
//       Ifpack2: Tempated Object-Oriented Algebraic Preconditioner Package
//                 Copyright (2009) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
//@HEADER
*/

//-----------------------------------------------------
// Ifpack2::KRYLOV is based on the Krylov iterative
// solvers in Belos.
// written by Paul Tsuji.
//-----------------------------------------------------

#ifndef IFPACK2_KRYLOV_DECL_HPP
#define IFPACK2_KRYLOV_DECL_HPP

#include "Ifpack2_ConfigDefs.hpp"
#include "Ifpack2_Preconditioner.hpp"
#include "Ifpack2_Condest.hpp"
#include "Ifpack2_Heap.hpp"
#include "Ifpack2_Parameters.hpp"

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>

#include <BelosConfigDefs.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosSolverManager.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosTpetraAdapter.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosBlockCGSolMgr.hpp>

#include <Teuchos_Assert.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_TypeNameTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>

#include <iostream>
#include <string>
#include <sstream>


#include <string>
#include <sstream>
#include <iostream>
#include <cmath>

namespace Teuchos {
  // forward declaration
  class ParameterList;
}

namespace Ifpack2 {

  //! A traits class for determining the scalar type to use for Belos
  /*
   * This exists to allow a ScalarType to override what scalar type is used
   * for Belos, if those are not equal, by specializing this class.
   */
  template <typename ScalarType>
  struct BelosScalarType {
    typedef ScalarType type;
  };
  
  //! A class for constructing and using a CG/GMRES smoother
  // for a given Tpetra::RowMatrix.
  
  /*! Ifpack2::Krylov computes a few iterations of CG/GMRES with zero
    initial guess as a smoother for a given Tpetra::RowMatrix.
    
    For all valid parameters, see the method Krylov::setParameters.
  */
  template<class MatrixType,class PrecType>
  class Krylov: virtual public Ifpack2::Preconditioner<typename MatrixType::scalar_type,typename MatrixType::local_ordinal_type,typename MatrixType::global_ordinal_type,typename MatrixType::node_type> {
    
  public:
    // \name Typedefs
    //@{
    
    //! The type of the entries of the input MatrixType.
    typedef typename MatrixType::scalar_type scalar_type;

    //! The scalar type used by Belos (which may not be scalar_type)
    typedef typename BelosScalarType<scalar_type>::type belos_scalar_type;
    
    //! Preserved only for backwards compatibility.  Please use "scalar_type".
    TEUCHOS_DEPRECATED typedef typename MatrixType::scalar_type Scalar;
    
    
    //! The type of local indices in the input MatrixType.
    typedef typename MatrixType::local_ordinal_type local_ordinal_type;
    
    //! Preserved only for backwards compatibility.  Please use "local_ordinal_type".
    TEUCHOS_DEPRECATED typedef typename MatrixType::local_ordinal_type LocalOrdinal;
    

    //! The type of global indices in the input MatrixType.
    typedef typename MatrixType::global_ordinal_type global_ordinal_type;
    
    //! Preserved only for backwards compatibility.  Please use "global_ordinal_type".
    TEUCHOS_DEPRECATED typedef typename MatrixType::global_ordinal_type GlobalOrdinal;
    
    
    //! The type of the Kokkos Node used by the input MatrixType.
    typedef typename MatrixType::node_type node_type;
    
    //! Preserved only for backwards compatibility.  Please use "node_type".
    TEUCHOS_DEPRECATED typedef typename MatrixType::node_type Node;


    //! The type of the magnitude (absolute value) of a matrix entry.
    typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;
    
    //! Preserved only for backwards compatibility.  Please use "magnitude_type".
    TEUCHOS_DEPRECATED typedef typename Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitudeType;

    //@}
    // \name Constructors and Destructors
    //@{
    
    //! Krylov explicit constuctor with Tpetra::RowMatrix input.
    explicit Krylov(const Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > &A);
    
    //! Krylov Destructor
    virtual ~Krylov();
    
    //@}
    //@{ Construction methods
    //! Set parameters for the preconditioner.
    /**
       <ul>
       <li> "krylov: number of iterations" (int)<br>
       <li> "krylov: residual tolerance" (double)<br>
       <li> "krylov: inner preconditioner" (string)<br>
       <li> "krylov: number of sweeps for inner preconditioning" (int)<br>
       <li> "krylov: damping parameter for inner preconditioning" (double)<br>
       </ul>
    */
    void setParameters(const Teuchos::ParameterList& params);
    
    //! Initialize Krylov preconditioner object.
    /*! Clear away any previously-allocated objects.
     */
    void initialize();
    
    //! Returns \c true if the preconditioner has been successfully initialized.
    inline bool isInitialized() const {
      return(IsInitialized_);
    }
    
    //! Setup iteration
    /*! Setups problem for Belos and computes any necessary preconditioners
      <ol>
      <li>
      <li>
      <li>
      </ol>
    */
    void compute();
    
    //! If compute() is completed, this query returns true, otherwise it returns false.
    inline bool isComputed() const {
      return(IsComputed_);
    }

    //@}

    //! @name Methods implementing Tpetra::Operator.
    //@{
    
    //! Returns the result of a few iterations of CG/GMRES on a Tpetra::MultiVector X in Y.
    /*!
      \param
      X - (In) A Tpetra::MultiVector of dimension NumVectors to solve for.
      \param
      Y - (Out) A Tpetra::MultiVector of dimension NumVectors containing result.
    */
    void apply(
	       const Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& X,
	       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>& Y,
	       Teuchos::ETransp mode = Teuchos::NO_TRANS,
	       scalar_type alpha = Teuchos::ScalarTraits<scalar_type>::one(),
	       scalar_type beta = Teuchos::ScalarTraits<scalar_type>::zero()) const;
    
    //! Tpetra::Map representing the domain of this operator.
    const Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> >& getDomainMap() const;
    
    //! Tpetra::Map representing the range of this operator.
    const Teuchos::RCP<const Tpetra::Map<local_ordinal_type,global_ordinal_type,node_type> >& getRangeMap() const;
    
    //! Whether this object's apply() method can apply the transpose (or conjugate transpose, if applicable).
    bool hasTransposeApply() const;
    
    //! Returns the computed estimated condition number, or -1.0 if no computed.
    magnitude_type getCondEst() const { return Condest_; }

    //@}

    //@{
    //! \name Mathematical functions.
    
    //! Returns the Tpetra::BlockMap object associated with the range of this matrix operator.
    const Teuchos::RCP<const Teuchos::Comm<int> > & getComm() const;

    //! Returns a reference to the matrix to be preconditioned.
    Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > getMatrix() const;
    
    //! Returns the number of calls to Initialize().
    int getNumInitialize() const;
    
    //! Returns the number of calls to Compute().
    int getNumCompute() const;

    //! Returns the number of calls to apply().
    int getNumApply() const;
    
    //! Returns the time spent in Initialize().
    double getInitializeTime() const;

    //! Returns the time spent in Compute().
    double getComputeTime() const;
    
    //! Returns the time spent in apply().
    double getApplyTime() const;

    //! Computes the estimated condition number and returns the value.
    magnitude_type computeCondEst(CondestType CT = Cheap,
				  local_ordinal_type MaxIters = 1550,
				  magnitude_type Tol = 1e-9,
				  const Teuchos::Ptr<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > &Matrix_in = Teuchos::null);
    
    //! @name Overridden from Teuchos::Describable
    //@{
    
    /** \brief Return a simple one-line description of this object. */
    std::string description() const;
    
    /** \brief Print the object with some verbosity level to an FancyOStream object. */
    void describe(Teuchos::FancyOStream &out, const Teuchos::EVerbosityLevel verbLevel=Teuchos::Describable::verbLevel_default) const;
    
    //@}
    
  private:
    
    // @{ Internal methods
    
    //! Copy constructor (should never be used)
    Krylov(const Krylov<MatrixType,PrecType>& RHS);
    
    //! operator= (should never be used)
    Krylov<MatrixType,PrecType>& operator=(const Krylov<MatrixType,PrecType>& RHS);
    
    //@}
    
    // @{ Internal data and parameters
    
    //! reference to the matrix to be preconditioned.
    const Teuchos::RCP<const Tpetra::RowMatrix<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > A_;
    //! Reference to the communicator object.
    const Teuchos::RCP<const Teuchos::Comm<int> > Comm_;

    //! General CG/GMRES parameters
    //!  Type of iteration
    // 1 for GMRES (default)
    // 2 for CG
    int IterationType_;
    //! Number of Iterations
    int Iterations_;
    //! Residual Tolerance
    magnitude_type ResidualTolerance_;
    //! Block size
    int BlockSize_;
    //! If true, the starting solution is always the zero vector.
    bool ZeroStartingSolution_;
    //! Preconditioner Type
    // 1 for relaxation (default)
    // 2 for ILUT
    // 3 for Chebyshev
    int PreconditionerType_;
    //! Teuchos parameter list
    Teuchos::ParameterList params_;

    //! Condition number estimate
    magnitude_type Condest_;
    //! \c true if \c this object has been initialized
    bool IsInitialized_;
    //! \c true if \c this object has been computed
    bool IsComputed_;
    //! Contains the number of successful calls to Initialize().
    int NumInitialize_;
    //! Contains the number of successful call to Compute().
    int NumCompute_;
    //! Contains the number of successful call to apply().
    mutable int NumApply_;
    //! Contains the time for all successful calls to Initialize().
    double InitializeTime_;
    //! Contains the time for all successful calls to Compute().
    double ComputeTime_;
    //! Contains the time for all successful calls to apply().
    mutable double ApplyTime_;
    //! Used for timing purposes
    mutable Teuchos::Time Time_;
    //! Number of local rows.
    local_ordinal_type NumMyRows_;
    //! Global number of nonzeros in L and U factors
    global_size_t NumGlobalNonzeros_;
    
    //! Belos Objects
    Teuchos::RCP< Belos::LinearProblem<belos_scalar_type,
				       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>,
				       Tpetra::Operator<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > > belosProblem_;
    Teuchos::RCP< Belos::SolverManager<belos_scalar_type,
				       Tpetra::MultiVector<scalar_type,local_ordinal_type,global_ordinal_type,node_type>,
				       Tpetra::Operator<scalar_type,local_ordinal_type,global_ordinal_type,node_type> > > belosSolver_;
    Teuchos::RCP<Teuchos::ParameterList> belosList_;
    Teuchos::RCP<PrecType> ifpack2_prec_;

  //@}

}; // class Krylov

}//namespace Ifpack2

#endif /* IFPACK2_KRYLOV_HPP */
