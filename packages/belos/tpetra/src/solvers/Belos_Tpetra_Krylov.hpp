#ifndef BELOS_TPETRA_KRYLOV_HPP
#define BELOS_TPETRA_KRYLOV_HPP

#include "Belos_Tpetra_Krylov_parameters.hpp"
#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_oblackholestream.hpp"
#include "Teuchos_ParameterList.hpp"

namespace BelosTpetra {
namespace Impl {

/// \brief Base class for iterative solvers that only work with Tpetra
///   objects, and that only solve for one right-hand-side vector at a
///   time.
///
/// \warning This is an implementation detail of Belos.
///
/// \tparam SC The type of of a dot product result.
/// \tparam MV Tpetra::MultiVector specialization.
/// \tparam OP Tpetra::MultiVector specialization.
///
/// This is a base class for Tpetra-based "non-block" solvers.  Here,
/// "block" means that the solver can solve for multiple right-hand
/// sides at a time.  Subclasses need only implement \c solveOneVec.
///
/// Users do not interact directly with this class.  Instead, when
/// users select the appropriate solver type, they get a
/// BelosTpetra::Impl::SolverManager instance from
/// Belos::SolverFactory.  The BelosTpetra::Impl::SolverManager
/// instance wraps an instance of a subclass of Krylov.
template<class SC = Tpetra::MultiVector<>::scalar_type,
	 class MV = Tpetra::MultiVector<SC>,
	 class OP = Tpetra::Operator<SC>>
class Krylov {
private:
  static_assert (std::is_same<MV, Tpetra::MultiVector<typename MV::scalar_type,
		 typename MV::local_ordinal_type,
		 typename MV::global_ordinal_type,
		 typename MV::node_type>>::value,
		 "MV must be a Tpetra::MultiVector specialization.");
  static_assert (std::is_same<OP, Tpetra::Operator<typename OP::scalar_type,
		 typename OP::local_ordinal_type,
		 typename OP::global_ordinal_type,
		 typename OP::node_type>>::value,
		 "OP must be a Tpetra::Operator specialization.");
public:
  using vec_type = Tpetra::Vector<typename MV::scalar_type,
				  typename MV::local_ordinal_type,
				  typename MV::global_ordinal_type,
				  typename MV::node_type>;
  
private:
  using STS = Teuchos::ScalarTraits<SC>;
  using mag_type = typename STS::magnitudeType;
  using STM = Teuchos::ScalarTraits<mag_type>;
  using device_type = typename MV::device_type;

public:
  Krylov () :
    verbosity_ (0)
  {}

  Krylov (const Teuchos::RCP<const OP>& A) :
    Krylov ()
  {
    A_ = A;
  }

  virtual ~Krylov () {}

  //! Set the matrix A in the linear system to solve.
  void setMatrix (const Teuchos::RCP<const OP>& A) {
    if (A_.get () != A.get ()) {
      A_ = A;
    }
  }

  void setLeftPrec (Teuchos::RCP<const OP> M) {
    input_.precoSide = "left";
    M_ = M;
  }

  void setRightPrec (Teuchos::RCP<const OP> M) {
    input_.precoSide = "right";
    M_ = M;
  }

  Teuchos::RCP<const OP>
  getMatrix () const {
    return A_;
  }

  Teuchos::RCP<const OP>
  getPreconditioner () const {
    return M_;
  }

  /// \brief Fill \c params with all parameters this solver
  ///   understands, and either their current values, or their default
  ///   values.
  ///
  /// \param params [out] To be filled with this solver's parameters.
  /// \param defaultValues [in] Whether to use default values (true)
  ///   or current values (false).
  virtual void
  getParameters (Teuchos::ParameterList& params,
                 const bool defaultValues) const
  {
    if (defaultValues) {
      getDefaultParameters (params);
    }
    else {
      const std::string implResScal = input_.needToScale ?
	"Norm of Preconditioned Initial Residual" : "None"; // ???
    
      params.set ("Convergence Tolerance", input_.tol);
      params.set ("Implicit Residual Scaling", implResScal);
      params.set ("Maximum Iterations", input_.maxNumIters);
      params.set ("Verbosity", verbosity_);
    }
  }

  virtual void
  getDefaultParameters (Teuchos::ParameterList& params) const
  {
    const SolverInput<SC> input;
    const int verbosity = 0;
    const std::string implResScal = input.needToScale ?
      "Norm of Preconditioned Initial Residual" : "None"; // ???
    
    params.set ("Convergence Tolerance", input.tol);
    params.set ("Implicit Residual Scaling", implResScal);
    params.set ("Maximum Iterations", input.maxNumIters);
    params.set ("Verbosity", verbosity);
  }

  /// \brief Set the solver's parameters.
  ///
  /// This solver takes a subset of the parameters that
  /// Belos::PseudoBlockCGSolMgr (Belos' CG implementation) takes, and
  /// ignores the rest.  The point is minimal sufficient compatibility
  /// with Belos' "generic" solvers.
  virtual void
  setParameters (Teuchos::ParameterList& params)
  {
    if (params.isParameter ("Convergence Tolerance")) {
      const mag_type tol = params.get<mag_type> ("Convergence Tolerance");
      TEUCHOS_TEST_FOR_EXCEPTION
        (tol < STM::zero (), std::invalid_argument,
         "\"Convergence tolerance\" = " << tol << " < 0.");
      input_.tol = tol;
    }

    if (params.isParameter ("Implicit Residual Scaling")) {
      const std::string implScal =
        params.get<std::string> ("Implicit Residual Scaling");
      if (implScal == "Norm of Initial Residual") {
        // FIXME (mfh 26 Oct 2016) Once we implement left
        // preconditioning, we'll have to keep separate preconditioned
        // and unpreconditioned absolute residuals.
        input_.needToScale = true;
      }
      else if (implScal == "Norm of Preconditioned Initial Residual") {
        input_.needToScale = true;
      }
      else if (implScal == "Norm of RHS") {
        // FIXME (mfh 26 Oct 2016) If we want to implement this, it
        // would make sense to combine that all-reduce with the
        // all-reduce for computing the initial residual norms.  We
        // could modify computeResiduals to have an option to do this.
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::logic_error,
           "\"Norm of RHS\" scaling option not implemented");
      }
      else if (implScal == "None") {
        input_.needToScale = false;
      }
      else {
        TEUCHOS_TEST_FOR_EXCEPTION
          (true, std::invalid_argument, "\"Implicit Residual Scaling\""
           " has an invalid value \"" << implScal << "\".");
      }
    }

    if (params.isParameter ("Maximum Iterations")) {
      const int maxNumIters = params.get<int> ("Maximum Iterations");
      TEUCHOS_TEST_FOR_EXCEPTION
        (maxNumIters < 0, std::invalid_argument,
         "\"Maximum Iterations\" = " << maxNumIters << " < 0.");
      input_.maxNumIters = maxNumIters;
    }

    int verbosity = verbosity_;
    if (params.isType<int> ("Verbosity")) {
      verbosity = params.get<int> ("Verbosity");
    }
    else if (params.isType<bool> ("Verbosity")) {
      const bool verbBool = params.get<bool> ("Verbosity");
      verbosity = verbBool ? 1 : 0;
    }
    verbosity_ = verbosity;
  }

private:
  static void
  computeResiduals (Kokkos::DualView<mag_type*, device_type>& norms,
                    MV& R,
                    const OP& A,
                    const MV& X,
                    const MV& B)
  {
    typedef typename device_type::memory_space dev_mem_space;

    const SC ONE = STS::one ();
    A.apply (X, R);
    R.update (ONE, B, -ONE); // R := B - A*X

    norms.template modify<dev_mem_space> ();
    Kokkos::View<mag_type*, device_type> norms_d =
      norms.template view<dev_mem_space> ();
    R.norm2 (norms_d);
  }

  static mag_type
  computeResidual (vec_type& R,
                   const OP& A,
                   const vec_type& X,
                   const vec_type& B)
  {
    const SC ONE = STS::one ();
    A.apply (X, R);
    R.update (ONE, B, -ONE); // R := B - A*X
    const mag_type r_norm = R.norm2 ();
    return r_norm;
  }

public:
  //! Solve the linear system(s) AX=B.
  SolverOutput<SC>
  solve (MV& X, const MV& B)
  {
    using Teuchos::FancyOStream;
    using Teuchos::getFancyOStream;
    using Teuchos::oblackholestream;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcpFromRef;
    using std::endl;

    TEUCHOS_TEST_FOR_EXCEPTION
      (A_.is_null (), std::runtime_error, "Matrix A is null.  Please call "
       "setMatrix() with a nonnull argument before calling solve().");
    TEUCHOS_TEST_FOR_EXCEPTION
      (X.getNumVectors () != B.getNumVectors (), std::runtime_error,
       "X.getNumVectors() = " << X.getNumVectors () <<
       " != B.getNumVectors() = " << B.getNumVectors () << ".");

    RCP<FancyOStream> outPtr;
    if (verbosity_) {
      const int myRank = A_->getDomainMap ()->getComm ()->getRank ();
      if (myRank == 0) {
        outPtr = getFancyOStream (rcpFromRef (std::cout));
      }
      else {
        outPtr = getFancyOStream (rcp (new oblackholestream ()));
      }
    }
    return solveImpl (outPtr.get (), X, B);
  }

  //! Return solver inputs
  SolverInput<SC> getSolverInput() {
    return input_;
  }

protected:
  static mag_type
  getConvergenceMetric (const mag_type r_norm_new,
                        const SolverInput<SC>& input)
  {
    if (input.needToScale) {
      return input.r_norm_orig == STM::zero () ?
        r_norm_new :
        (r_norm_new / input.r_norm_orig);
    }
    else {
      return r_norm_new;
    }
  }

  static mag_type
  getConvergenceMetric (const mag_type r_norm_new,
                        const mag_type r_norm_orig,
                        const SolverInput<SC>& input)
  {
    if (input.needToScale) {
      return r_norm_orig == STM::zero () ?
        r_norm_new :
        (r_norm_new / r_norm_orig);
    }
    else {
      return r_norm_new;
    }
  }

private:
  //! Solve the linear system(s) AX=B.
  SolverOutput<SC>
  solveImpl (Teuchos::FancyOStream* outPtr, MV& X, const MV& B)
  {
    using Teuchos::RCP;
    using std::endl;

    TEUCHOS_TEST_FOR_EXCEPTION
      (M_.get () == nullptr && input_.precoSide != "none",
       std::logic_error, "BelosTpetra::Impl::Krylov::solveImpl: M_ is null "
       "but input_.precoSide=\"" << input_.precoSide << "\" != \"none\".  "
       "This should never happen.  "
       "Please report this bug to the Belos developers.");
    
    const size_t numVecs = B.getNumVectors ();
    Kokkos::DualView<mag_type*, device_type> norms ("norms", numVecs);
    MV R (B.getMap (), numVecs);

    computeResiduals (norms, R, *A_, X, B);
    norms.template sync<Kokkos::HostSpace> ();
    auto norms_h = norms.template view<Kokkos::HostSpace> ();
    SolverOutput<SC> allOutput {};
    for (size_t j = 0; j < numVecs; ++j) {
      if (outPtr != nullptr) {
        *outPtr << "Solve for column " << (j+1) << " of " << numVecs << ":"
		<< endl;
        outPtr->pushTab ();
      }
      RCP<vec_type> R_j = R.getVectorNonConst (j);
      RCP<vec_type> X_j = X.getVectorNonConst (j);
      input_.r_norm_orig = norms_h(j);
      const SolverOutput<SC> curOutput =
        solveOneVec (outPtr, *X_j, *R_j, *A_,
		     (input_.precoSide == "none" ? *A_ : *M_), input_);
      combineSolverOutput (allOutput, curOutput);
    }
    return allOutput;
  }
  
protected:
  virtual SolverOutput<SC>
  solveOneVec (Teuchos::FancyOStream* outPtr,
               vec_type& X, // in X/out X
               vec_type& R, // in B/out R
               const OP& A,
               const OP& M,
               const SolverInput<SC>& input) = 0;

protected:
  SolverInput<SC> input_;

private:
  Teuchos::RCP<const OP> A_;
  std::string precoType_;
  Teuchos::RCP<const OP> M_;
  int verbosity_;
};

} // namespace Impl
} // namespace BelosTpetra

#endif // BELOS_TPETRA_KRYLOV_HPP
