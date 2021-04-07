#ifndef BELOS_TPETRA_KRYLOV_PARAMETERS_HPP
#define BELOS_TPETRA_KRYLOV_PARAMETERS_HPP

#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

namespace BelosTpetra {
namespace Impl {

class Indent {
public:
  Indent (Teuchos::FancyOStream* out) : out_ (out) {
    if (out_ != nullptr) {
      out_->pushTab ();
    }
  }

  Indent () = delete;
  Indent (const Indent&) = delete;
  Indent& operator= (const Indent&) = delete;

  ~Indent () {
    if (out_ != nullptr) {
      out_->popTab ();
    }
  }
private:
  Teuchos::FancyOStream* out_;
};

template<class ScalarType>
struct SolverInput {
public:
  using mag_type = typename Teuchos::ScalarTraits<ScalarType>::magnitudeType;

private:
  using STS = Teuchos::ScalarTraits<ScalarType>;
  using STM = Teuchos::ScalarTraits<mag_type>;

public:
  SolverInput () = default;
  SolverInput (const SolverInput<ScalarType>&) = default;
  SolverInput& operator= (const SolverInput<ScalarType>&) = default;

  mag_type r_norm_orig = STM::zero ();
  mag_type tol = STM::squareroot (STS::eps ());
  int maxNumIters = 1000;
  int resCycle = 30;
  int stepSize = 5;
  bool needToScale = true;
  bool needToReortho = false;
  int maxOrthoSteps = 0;
  std::string orthoType {"ICGS"};
  std::string precoSide {"none"};
  bool computeRitzValues = true;
  bool computeRitzValuesOnFly = false;
};

// The default constructor creates output corresponding to "solving
// a linear system with zero right-hand sides."  This means that
// the solve succeeded trivially (converged == true), in zero
// iterations, with zero residual norm.
template<class SC>
struct SolverOutput {
  using val_type = SC;
  using mag_type = typename Teuchos::ScalarTraits<SC>::magnitudeType;
  using complex_type = std::complex<mag_type>;

  SolverOutput () = default;
  SolverOutput (const SolverOutput<SC>&) = default;
  SolverOutput& operator= (const SolverOutput<SC>&) = default;

  //! Absolute residual norm.
  mag_type absResid = Teuchos::ScalarTraits<mag_type>::zero ();
  //! Relative residual norm (if applicable).
  mag_type relResid = Teuchos::ScalarTraits<mag_type>::zero ();
  //! Number of iterations executed.
  int numIters = 0;
  //! Number of restarts.
  int numRests = -1;
  //! Whether the solve converged.
  bool converged = true;
  //! Ritz values if requested
  std::vector<complex_type> ritzValues;
};

/// \brief Combine two solver outputs.
///
/// "Combining" two solver outputs means aggregating the results of
/// solving A x_1 = b_1 and A x_2 = b_2, that is, two solves with the
/// same matrix, but different right-hand sides.  Combining is
/// associative and commutative.
template<class SC>
void
combineSolverOutput (SolverOutput<SC>& dst, const SolverOutput<SC>& src)
{
  // max of the residuals and iteration counts
  dst.relResid = dst.relResid > src.relResid ? dst.relResid : src.relResid;
  dst.absResid = dst.absResid > src.absResid ? dst.absResid : src.absResid;
  dst.numIters = dst.numIters > src.numIters ? dst.numIters : src.numIters;
  dst.numRests = dst.numRests > src.numRests ? dst.numRests : src.numRests;
  // "converged" if all converged
  dst.converged = dst.converged && src.converged;
  // copy ritz values
  if (src.ritzValues.size () > dst.ritzValues.size ()) {
    dst.ritzValues.resize (src.ritzValues.size ());
    std::copy (std::begin (src.ritzValues), std::end (src.ritzValues),
               std::begin (dst.ritzValues));
  }
}

template<class SC>
std::ostream&
operator<< (std::ostream& out,
            const SolverOutput<SC>& so)
{
  using std::endl;

  out << "Solver output:" << endl
      << " Absolute residual norm: " << so.absResid << endl
      << " Relative residual norm: " << so.relResid << endl
      << " Number of iterations: " << so.numIters << endl;
  if (so.numRests >= 0) {
    out << " Number of restarts: " << so.numRests << endl;
  }
  out << " Converged: " << (so.converged ? "true" : "false") << endl;
  if (so.ritzValues.size () != 0) {
    out << " Ritz values: [";
    for (std::size_t k = 0; k < so.ritzValues.size (); ++k) {
      out << so.ritzValues[k];
      if (k + std::size_t (1) < so.ritzValues.size ()) {
        out << ", ";
      }
    }
    out << "]" << endl;
  }
  return out;
}

template<class SC>
std::ostream&
operator<< (std::ostream& out,
            const SolverInput<SC>& si)
{
  using std::endl;

  out << "Solver input:" << endl;
  out << " Original residual norm: " << si.r_norm_orig << endl
      << " Residual norm tolerance: " << si.tol << endl
      << " Max. number of iterations: " << si.maxNumIters << endl;
  if (si.resCycle > 0) {
    out << " Restart cycle: " << si.resCycle << endl;
    out << " Orthogonalization: " << si.orthoType << endl;
  }
  out << " Step size: " << si.stepSize << endl
      << " Preconditioner: " << si.precoSide << endl;
  return out;
}

} // namespace Impl
} // namespace BelosTpetra

#endif // BELOS_TPETRA_KRYLOV_PARAMETERS_HPP
