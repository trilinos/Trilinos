#ifndef EDGE_INTERSECTION_PRIMITIVE_H
#define EDGE_INTERSECTION_PRIMITIVE_H

#include <stk_middle_mesh/projection.hpp>
#include <cassert>
#include <memory>

namespace stk {
namespace middle_mesh {
namespace predicates {
namespace impl {

class EdgeIntersectionResult
{
  public:
    EdgeIntersectionResult(double alpha, double beta);
    EdgeIntersectionResult(double alpha1, double beta1, double alpha2, double beta2);

    EdgeIntersectionResult();

    bool intersection1_found() const { return m_intersection1Found; }
    double get_alpha1() const
    {
      assert(m_intersection1Found);
      return m_alpha1;
    }

    double get_beta1() const
    {
      assert(m_intersection1Found);
      return m_beta1;
    }

    bool intersection2_found() const { return m_intersection2Found; }
    double get_alpha2() const
    {
      assert(m_intersection2Found);
      return m_alpha2;
    }

    double get_beta2() const
    {
      assert(m_intersection2Found);
      return m_beta2;
    }

  private:
    bool m_intersection1Found = false;
    double m_alpha1           = -1;
    double m_beta1            = -1;

    bool m_intersection2Found = false;
    double m_alpha2           = -1;
    double m_beta2            = -1;
};

class EdgeIntersectionData
{
  public:
    // d0 and d1 are the direction vectors at g0 and g1, respectively
    EdgeIntersectionData(const utils::Point& b0, const utils::Point& b1, const utils::Point& g0, const utils::Point& g1,
                         const utils::Point& d0, const utils::Point& d1, double tol = 1e-12, double errorTol = 1e-12);

    void compute_beta();
    void compute_alpha();
    double solve_alpha_linear_system(double beta, bool& alphaValid);

    bool is_beta_1_valid() { return m_beta1Valid; }
    bool is_beta_2_valid() { return m_beta2Valid; }
    bool is_alpha_1_valid() { return m_alpha1Valid; }
    bool is_alpha_2_valid() { return m_alpha2Valid; }

    double get_beta_1_value() { return m_beta1; }
    double get_beta_2_value() { return m_beta2; }
    double get_alpha_1_value() { return m_alpha1; }
    double get_alpha_2_value() { return m_alpha2; }

    // returns the cosine of the angle between the direction vector
    // and the displacement vector
    double compute_error_angle(double alpha, double beta);

    EdgeIntersectionResult compute_edge_intersection_result();

    void compute_solution_in_weird_way();

  private:
    utils::Point m_b0;
    utils::Point m_b1;
    utils::Point m_g0;
    utils::Point m_g1;
    utils::Point m_d0;
    utils::Point m_d1;
    bool m_beta1Valid;
    bool m_beta2Valid;
    bool m_alpha1Valid;
    bool m_alpha2Valid;
    double m_beta1;
    double m_beta2;
    double m_alpha1;
    double m_alpha2;
    double m_tol;
    double m_errorTol;
};

EdgeIntersectionResult compute_edge_intersection(const utils::Point& b0, const utils::Point& b1, const utils::Point& g0,
                                                 const utils::Point& g1, const utils::Point& d0, const utils::Point& d1,
                                                 double tol = 1e-12, double errorTol = 1e-12);

} // namespace impl

} // namespace predicates
} // namespace middle_mesh
} // namespace stk
#endif
