#include "patch_distortion_objective.hpp"

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

double PatchDistortionObjective::compute_quality(const utils::Point& ptIn)
{
  utils::impl::Matrix<utils::Point> triPts;
  std::vector<utils::impl::Mat2x2<double>> wMats;
  compute_parameterization(ptIn, triPts, wMats);

  double q = 0;
  for (int i = 0; i < m_data->get_num_elements(); ++i)
  {
    TriPts ptsI{triPts(i, 0), triPts(i, 1), triPts(i, 2)};

    double qI = m_metric->get_value(ptsI, wMats[i]);

    // apply constraint aggregation function
    q += qI * qI;
  }

  q = std::sqrt(q / m_data->get_num_elements());

  return q;
}

// computes derivative of computeQuality wrt the active vert only
// Computes derivative wrt the parameterized coords
// pt is the current position of the active vert
utils::Point PatchDistortionObjective::compute_quality_rev(const utils::Point& ptIn, double q3Bar)
{
  // std::cout << "\nEntered computeQuality_rev" << std::endl;
  utils::impl::Matrix<utils::Point> triPts;
  std::vector<utils::impl::Mat2x2<double>> wMats;
  compute_parameterization(ptIn, triPts, wMats);

  double q = 0;
  for (int i = 0; i < m_data->get_num_elements(); ++i)
  // for (int i=0; i < 1; ++i)
  {
    TriPts ptsI{triPts(i, 0), triPts(i, 1), triPts(i, 2)};

    //    auto W = DistortionMetric::computeW(pts_i);
    //    if (det2x2(W) > 0)
    //      continue;

    double qI = m_metric->get_value(ptsI, wMats[i]);

    // apply constraint aggregation function
    q += qI * qI;
  }

  double q2 = q / m_data->get_num_elements();
  // double q3 = std::sqrt(q2);

  //-------------------------------
  // reverse
  double q2Bar = 0.5 * q3Bar / std::sqrt(q2);
  double qBar  = q2Bar / m_data->get_num_elements();

  utils::Point activePtBar;
  for (int i = 0; i < m_data->get_num_elements(); ++i)
  // for (int i=0; i < 1; ++i)
  {
    TriPts ptsI{triPts(i, 0), triPts(i, 1), triPts(i, 2)};
    //    auto W = DistortionMetric::computeW(pts_i);
    //    if (det2x2(W) > 0)
    //      continue;

    double qI = m_metric->get_value(ptsI, wMats[i]);

    // apply constraint aggregation function
    // q += q_i*q_i;

    double qIBar = 2 * qI * qBar;
    TriPts ptsIBar;
    m_metric->get_deriv(ptsI, wMats[i], ptsIBar, qIBar);

    // accumulate only active vert
    int idx = m_data->get_current_vert_idx(i);
    activePtBar += ptsIBar[idx];
  }

  return activePtBar;
}

utils::impl::Mat2x2<double> PatchDistortionObjective::compute_hessian(const utils::Point& ptIn)
{
  utils::impl::Matrix<utils::Point> triPts;
  std::vector<utils::impl::Mat2x2<double>> wMats;
  compute_parameterization(ptIn, triPts, wMats);

  using DotVec = DistortionMetric::DotVec;
  using PtsDot = DistortionMetric::PtsDot;

  double q = 0;
  DotVec qDot{};
  for (int i = 0; i < m_data->get_num_elements(); ++i)
  // for (int i=0; i < 1; ++i)
  {
    TriPts ptsI{triPts(i, 0), triPts(i, 1), triPts(i, 2)}, ptsIDot;

    //    auto W = DistortionMetric::computeW(pts_i);
    //    if (det2x2(W) > 0)
    //      continue;

    double qI = m_metric->get_value(ptsI, wMats[i]);
    m_metric->get_deriv(ptsI, wMats[i], ptsIDot, 1);

    // apply constraint aggregation function
    q += qI * qI;

    int idx = m_data->get_current_vert_idx(i);
    qDot[0] += 2 * qI * ptsIDot[idx].x;
    qDot[1] += 2 * qI * ptsIDot[idx].y;
  }

  DotVec q2Dot, q3Dot;
  double q2 = q / m_data->get_num_elements();
  q2Dot[0]  = qDot[0] / m_data->get_num_elements();
  q2Dot[1]  = qDot[1] / m_data->get_num_elements();
  // double q3 = std::sqrt(q2);
  q3Dot[0] = 0.5 * q2Dot[0] / std::sqrt(q2);
  q3Dot[1] = 0.5 * q2Dot[1] / std::sqrt(q2);

  //-------------------------------
  // reverse
  double q3Bar = 1;
  double q2Bar = 0.5 * q3Bar / std::sqrt(q2);
  double qBar  = q2Bar / m_data->get_num_elements();

  DotVec qBarDot;
  for (unsigned int i = 0; i < qBarDot.size(); ++i)
  {
    auto q2BarDot = -0.25 * q3Bar * q2Dot[i] / std::pow(q2, 1.5);
    qBarDot[i]    = q2BarDot / m_data->get_num_elements();
  }

  utils::Point activePtBar;
  utils::impl::Mat2x2<double> activePtBarDot; // Hessian
                                              //  active_pt_bar_dot(0, 0) = 0;
                                              //  active_pt_bar_dot(0, 1) = 0;
                                              //  active_pt_bar_dot(1, 0) = 0;
                                              //  active_pt_bar_dot(1, 1) = 0;

  for (int i = 0; i < m_data->get_num_elements(); ++i)
  // for (int i=0; i < 1; ++i)
  {
    TriPts ptsI{triPts(i, 0), triPts(i, 1), triPts(i, 2)}, ptsIBar;
    //    auto W = DistortionMetric::computeW(pts_i);
    //    if (det2x2(W) > 0)
    //      continue;

    double qI = m_metric->get_value(ptsI, wMats[i]);
    m_metric->get_deriv(ptsI, wMats[i], ptsIBar, 1);

    int idx = m_data->get_current_vert_idx(i);
    DotVec qIDot;
    qIDot[0] = ptsIBar[idx].x;
    qIDot[1] = ptsIBar[idx].y;

    // apply constraint aggregation function
    // q += q_i*q_i;

    double qIBar = 2 * qI * qBar;
    DotVec qIBarDot;
    qIBarDot[0] = 2 * qBar * qIDot[0] + 2 * qI * qBarDot[0];
    qIBarDot[1] = 2 * qBar * qIDot[1] + 2 * qI * qBarDot[1];

    TriPts ptsIBar2;
    m_metric->get_deriv(ptsI, wMats[i], ptsIBar2, qIBar);

    PtsDot ptsBarDot;
    PtsDot ptsIDot;
    for (int j = 0; j < 6; ++j)
      for (int k = 0; k < 2; ++k)
        ptsIDot[j][k] = 0;

    ptsIDot[2 * idx][0]     = 1;
    ptsIDot[2 * idx + 1][1] = 1;

    m_metric->get_value_rev_dot(ptsI, wMats[i], ptsIDot, qIBar, qIBarDot, ptsBarDot);

    // accumulate only active vert
    activePtBar += ptsIBar2[idx];

    activePtBarDot(0, 0) += ptsBarDot[2 * idx][0];
    activePtBarDot(0, 1) += ptsBarDot[2 * idx][1];
    activePtBarDot(1, 0) += ptsBarDot[2 * idx + 1][0];
    activePtBarDot(1, 1) += ptsBarDot[2 * idx + 1][1];
  }

  return activePtBarDot;
}

} // namespace impl
} // namespace mesh
} // namespace middle_mesh
} // namespace stk
