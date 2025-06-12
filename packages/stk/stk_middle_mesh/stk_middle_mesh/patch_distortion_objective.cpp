#include "patch_distortion_objective.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

double PatchDistortionObjective::compute_quality(const utils::Point& ptIn)
{
  set_min_denominator(ptIn);
  return compute_qualityT(m_metric, ptIn);
}

utils::Point PatchDistortionObjective::compute_quality_rev(const utils::Point& ptIn, double q3Bar)
{
  set_min_denominator(ptIn);
  utils::Point grad = compute_quality_revT(m_metric, ptIn, q3Bar);

  return grad;
}

utils::impl::Mat2x2<double> PatchDistortionObjective::compute_hessian(const utils::Point& ptIn)
{
  set_min_denominator(ptIn);
  utils::impl::Mat2x2<double> hess = compute_hessianT(m_metric, ptIn);

  return hess;
}


template <typename T>
T PatchDistortionObjective::compute_qualityT(std::shared_ptr<DistortionMetric<T>> /*metric*/, const utils::PointT<T>& ptIn)
{
  utils::impl::Matrix<utils::PointT<T>> triPts;
  std::vector<utils::impl::Mat2x2<double>> wMats;
  compute_parameterization(ptIn, triPts, wMats);

  T q = 0;
  for (int i = 0; i < m_data->get_num_elements(); ++i)
  {
    TriPts<T> ptsI{triPts(i, 0), triPts(i, 1), triPts(i, 2)};

    T qI = m_metric->get_value(ptsI, wMats[i]);

    // apply constraint aggregation function
    q += qI * qI;
  }

  q = std::sqrt(q / m_data->get_num_elements());

  return q;
}


// computes derivative of computeQuality wrt the active vert only
// Computes derivative wrt the parameterized coords
// pt is the current position of the active vert
template <typename T>
utils::PointT<T> PatchDistortionObjective::compute_quality_revT(std::shared_ptr<DistortionMetric<T>> metric, const utils::PointT<T>& ptIn, T q3Bar)
{
  utils::impl::Matrix<utils::PointT<T>> triPts;
  std::vector<utils::impl::Mat2x2<double>> wMats;
  compute_parameterization(ptIn, triPts, wMats);

  T q = 0;
  for (int i = 0; i < m_data->get_num_elements(); ++i)
  {
    TriPts<T> ptsI{triPts(i, 0), triPts(i, 1), triPts(i, 2)};

    T qI = metric->get_value(ptsI, wMats[i]);

    // apply constraint aggregation function
    q += qI * qI;
  }

  T q2 = q / m_data->get_num_elements();

  //-------------------------------
  // reverse
  T q2Bar = 0.5 * q3Bar / std::sqrt(q2);
  T qBar  = q2Bar / m_data->get_num_elements();

  utils::PointT<T> activePtBar;
  for (int i = 0; i < m_data->get_num_elements(); ++i)
  {    
    TriPts<T> ptsI{triPts(i, 0), triPts(i, 1), triPts(i, 2)};

    T qI = metric->get_value(ptsI, wMats[i]);

    // apply constraint aggregation function
    // q += q_i*q_i;

    T qIBar = 2 * qI * qBar;
    TriPts<T> ptsIBar;
    m_metric->get_deriv(ptsI, wMats[i], ptsIBar, qIBar);

    // accumulate only active vert
    int idx = m_data->get_current_vert_idx(i);

    activePtBar += ptsIBar[idx];
  }

  return activePtBar;
}


template <typename T>
utils::impl::Mat2x2<T> PatchDistortionObjective::compute_hessianT(std::shared_ptr<DistortionMetric<T>> metric, const utils::PointT<T>& ptIn)
{
  utils::impl::Matrix<utils::PointT<T>> triPts;
  std::vector<utils::impl::Mat2x2<double>> wMats;
  compute_parameterization(ptIn, triPts, wMats);

  using DotVec = DistortionMetric<double>::DotVec<T>;
  using PtsDot = DistortionMetric<double>::PtsDot<T>;

  double q = 0;
  DotVec qDot{};
  for (int i = 0; i < m_data->get_num_elements(); ++i)
  {
    TriPts<T> ptsI{triPts(i, 0), triPts(i, 1), triPts(i, 2)}, ptsIDot;

    double qI = metric->get_value(ptsI, wMats[i]);
    metric->get_deriv(ptsI, wMats[i], ptsIDot, 1);

    // apply constraint aggregation function
    q += qI * qI;

    int idx = m_data->get_current_vert_idx(i);
    qDot[0] += 2 * qI * ptsIDot[idx].x;
    qDot[1] += 2 * qI * ptsIDot[idx].y;
  }

  DotVec q2Dot, q3Dot;
  T q2 = q / m_data->get_num_elements();
  q2Dot[0]  = qDot[0] / m_data->get_num_elements();
  q2Dot[1]  = qDot[1] / m_data->get_num_elements();
  // double q3 = std::sqrt(q2);
  q3Dot[0] = 0.5 * q2Dot[0] / std::sqrt(q2);
  q3Dot[1] = 0.5 * q2Dot[1] / std::sqrt(q2);

  //-------------------------------
  // reverse
  T q3Bar = 1;
  T q2Bar = 0.5 * q3Bar / std::sqrt(q2);
  T qBar  = q2Bar / m_data->get_num_elements();

  DotVec qBarDot;
  for (unsigned int i = 0; i < qBarDot.size(); ++i)
  {
    auto q2BarDot = -0.25 * q3Bar * q2Dot[i] / std::pow(q2, 1.5);
    qBarDot[i]    = q2BarDot / m_data->get_num_elements();
  }

  utils::Point activePtBar;
  utils::impl::Mat2x2<T> activePtBarDot; // Hessian

  for (int i = 0; i < m_data->get_num_elements(); ++i)
  {
    TriPts<T> ptsI{triPts(i, 0), triPts(i, 1), triPts(i, 2)}, ptsIBar;

    T qI = metric->get_value(ptsI, wMats[i]);
    metric->get_deriv(ptsI, wMats[i], ptsIBar, 1);

    int idx = m_data->get_current_vert_idx(i);
    DotVec qIDot;
    qIDot[0] = ptsIBar[idx].x;
    qIDot[1] = ptsIBar[idx].y;

    // apply constraint aggregation function
    // q += q_i*q_i;

    T qIBar = 2 * qI * qBar;
    DotVec qIBarDot;
    qIBarDot[0] = 2 * qBar * qIDot[0] + 2 * qI * qBarDot[0];
    qIBarDot[1] = 2 * qBar * qIDot[1] + 2 * qI * qBarDot[1];


    TriPts<T> ptsIBar2;
    metric->get_deriv(ptsI, wMats[i], ptsIBar2, qIBar);

    PtsDot ptsBarDot;
    PtsDot ptsIDot;
    for (int j = 0; j < 6; ++j)
      for (int k = 0; k < 2; ++k)
        ptsIDot[j][k] = 0;

    ptsIDot[2 * idx][0]     = 1;
    ptsIDot[2 * idx + 1][1] = 1;

    metric->get_value_rev_dot(ptsI, wMats[i], ptsIDot, qIBar, qIBarDot, ptsBarDot);

    // accumulate only active vert
    activePtBar += ptsIBar2[idx];

    activePtBarDot(0, 0) += ptsBarDot[2 * idx][0];
    activePtBarDot(0, 1) += ptsBarDot[2 * idx][1];
    activePtBarDot(1, 0) += ptsBarDot[2 * idx + 1][0];
    activePtBarDot(1, 1) += ptsBarDot[2 * idx + 1][1];
  }

  return activePtBarDot;
}


void PatchDistortionObjective::set_min_denominator(const utils::Point& ptIn)
{
  // calling compute_parameterization is duplicitive, but we don't
  // want the gradient calculation to include the effect of
  // this parameter, so it has to occur outside the compute_quality function
  utils::impl::Matrix<utils::PointT<double>> triPts;
  std::vector<utils::impl::Mat2x2<double>> wMats;
  compute_parameterization(ptIn, triPts, wMats);

  double minDen = std::numeric_limits<double>::max();
  for (int i=0; i < triPts.extent0(); ++i)
  {
    TriPts<double> ptsI{triPts(i, 0), triPts(i, 1), triPts(i, 2)};
    double den = m_metric->compute_denominator(ptsI, wMats[i]);
    minDen = std::min(den, minDen);
  }

  m_metric->set_min_denominator(minDen);
}


} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk
