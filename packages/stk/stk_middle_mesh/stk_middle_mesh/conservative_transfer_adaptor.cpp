#include "conservative_transfer_adaptor.hpp"

namespace stk {
namespace middle_mesh {
namespace disc1d {
namespace impl {

void ConservativeTransferAdaptor::correct_coeffs()
{
  m_interpolator->get_interpolation_matrix(m_rowCoeffs, m_rowIndices);
  m_srcMass  = compute_summed_mass_matrix(m_discSrc);
  m_destMass = compute_summed_mass_matrix(m_discDest);

  get_column_indices();

  std::vector<double> coeffsJ, optCoeffsJ, destMassValues;
  double maxDiffNorm = 0, maxConstraintViolation = 0;
  for (int j = 0; j < m_discSrc->get_num_dofs(); ++j)
  {
    auto indices = m_columnIndices[j];
    int ncoeffs  = indices.size();

    coeffsJ.resize(ncoeffs);
    optCoeffsJ.resize(ncoeffs);
    destMassValues.resize(ncoeffs);

    for (int i = 0; i < ncoeffs; ++i)
    {
      CSRIndex idx      = indices[i];
      coeffsJ[i]        = m_rowCoeffs[idx.row][idx.colIdx];
      destMassValues[i] = m_destMass[idx.row];
    }
    double srcMassValue = m_srcMass[j];

    double constraintViolation = std::abs(dot(destMassValues, coeffsJ) - srcMassValue);
    maxConstraintViolation     = std::max(constraintViolation, maxConstraintViolation);

    solve_quadratic_program(coeffsJ, destMassValues, srcMassValue, optCoeffsJ);

    // compute norm of change in coeffs;
    double diffNorm = 0;
    for (int i = 0; i < ncoeffs; ++i)
    {
      double diff = std::abs(optCoeffsJ[i] - coeffsJ[i]);
      diffNorm += diff * diff;
    }
    diffNorm    = std::sqrt(diffNorm);
    maxDiffNorm = std::max(diffNorm, maxDiffNorm);

    for (int i = 0; i < ncoeffs; ++i)
    {
      CSRIndex idx                     = indices[i];
      m_rowCoeffs[idx.row][idx.colIdx] = optCoeffsJ[i];
    }
  }

  // std::cout << "max constraint violation    = " << max_constraint_violation << std::endl;
  // std::cout << "max coefficient change norm = " << max_diff_norm << std::endl;
}

void ConservativeTransferAdaptor::get_column_indices()
{
  m_columnIndices.resize(m_discSrc->get_num_dofs());
  for (int row = 0; row < int(m_rowIndices.size()); ++row)
  {
    auto& columns = m_rowIndices[row];
    for (int columnIdx = 0; columnIdx < int(columns.size()); ++columnIdx)
    {
      int column = columns[columnIdx];
      m_columnIndices[column].push_back(CSRIndex{row, columnIdx});
    }
  }

  // TODO: sort them?
}

// solves a quadratic program for the special case where G = I and A is a single row
void ConservativeTransferAdaptor::solve_quadratic_program(const std::vector<double>& x0, const std::vector<double>& a,
                                                          double b, std::vector<double>& x)
{
  // Note: take x = x0 initially (the initial value of x is arbitrary)
  //       Maybe x=0 would be better
  // TODO: dynamic memory allocation in std::vector
  assert(x0.size() == x.size());
  assert(a.size() == x0.size());
  /*
    std::cout << "\nconstraint matrix = ";
    for (auto& a : A)
      std::cout << a << ", ";
    std::cout << std::endl;

    std::cout << "b = " << b << std::endl;
    std::cout << "initial constraint violation = " << (dot(A, x0)  - b) << std::endl;
  */
  int nvars = x0.size();
  std::vector<double> g(nvars);
  for (int i = 0; i < nvars; ++i)
    g[i] = -2 * x0[i] + x0[i];
  double h = dot(a, x0) - b;

  // solve for lambda
  double lambdaRhs = dot(a, g) - h;
  double lambdaLhs = dot(a, a);
  double lambda    = lambdaRhs / lambdaLhs;

  // Note: when G = I, we can compute the step p directly
  std::vector<double> p(nvars);
  for (int i = 0; i < nvars; ++i)
    p[i] = a[i] * lambda - g[i];

  for (int i = 0; i < nvars; ++i)
    x[i] = p[i] + x0[i];

  // verify constraints are satisfied
  //      double residual = dot(A, x) - b;
  /*
    std::cout << "final constraint violation = " << residual << std::endl;

    std::cout << "solution:" << std::endl;
    for (int i=0; i < nvars; ++i)
      std::cout << "coeff " << i << ": initial " << x0[i] << ", final " << x[i] << std::endl;
  */
  //     if (std::abs(residual) > 1e-13)
  //       throw std::runtime_error("constraints are not satisfied");

  // verify overall problem is solved
  std::vector<double> block1(nvars);
  double block2;

  for (int i = 0; i < nvars; ++i)
    block1[i] = -p[i] + a[i] * lambda + 2 * x0[i] - x0[i];
  block2 = -dot(a, p) - dot(a, x0) + b;

  double block1Norm = std::sqrt(dot(block1, block1));
  double block2Norm = std::abs(block2);

  // std::cout << "block1_norm = " << block1_norm << ", block2_norm = " << block2_norm << std::endl;
  if (block1Norm > 1e-13 || block2Norm > 1e-13)
    throw std::runtime_error("quadratic program not solved correctly");
}

double ConservativeTransferAdaptor::dot(const std::vector<double>& a, const std::vector<double>& b)
{
  assert(a.size() == b.size());

  double c = 0;
  for (size_t i = 0; i < a.size(); ++i)
    c += a[i] * b[i];

  return c;
}
} // namespace impl
} // namespace disc1d
} // namespace middle_mesh
} // namespace stk
