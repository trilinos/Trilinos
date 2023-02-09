#ifndef BASIS_1D_H
#define BASIS_1D_H

#include <vector>

namespace stk {
namespace middle_mesh {
namespace disc1d {
namespace impl {

// TODO: move to own file
class BasisFunctions
{
  public:
    explicit BasisFunctions(int degree);

    double get_value(double xi, int node) const;

    double get_derivative(double xi, int node) const;

    int get_num_nodes() const { return m_nodeXi.size(); }

    int get_degree() const { return m_degree; }

    const std::vector<double>& get_node_xi() { return m_nodeXi; }

    double get_xi_min() const { return -1; }

    double get_xi_max() const { return 1; }

  private:
    double get_value_degree1(double xi, int node) const;

    double get_deriv_degree1(double /*xi*/, int node) const;

    double get_value_degree2(double xi, int node) const;

    double get_deriv_degree2(double xi, int node) const;

    void check_xi(double xi) const;

    int m_degree;
    std::vector<double> m_nodeXi;
};

class BasisFunctionValues
{
  public:
    BasisFunctionValues(BasisFunctions& basis, const std::vector<double>& xiPoints);

    double get_value(int node, int point) const { return m_vals[get_idx(node, point)]; }

    double get_deriv(int node, int point) const { return m_derivs[get_idx(node, point)]; }

    int get_num_nodes() const { return m_numNodes; }

    int get_num_points() const { return m_numPoints; }

  private:
    int get_idx(int node, int point) const { return node * m_numPoints + point; }

    std::vector<double> m_vals;
    std::vector<double> m_derivs;
    int m_numNodes;
    int m_numPoints;
};

} // namespace impl

} // namespace disc1d
} // namespace middle_mesh
} // namespace stk
#endif