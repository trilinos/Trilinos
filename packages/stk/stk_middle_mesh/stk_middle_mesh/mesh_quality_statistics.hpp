#ifndef MESH_QUALITY_STATISTICS_H
#define MESH_QUALITY_STATISTICS_H

#include <memory>

#include "mesh.hpp"
#include "patch_objective.hpp"
#include <fstream>

namespace stk {
namespace middle_mesh {
namespace mesh {
namespace impl {

class MeshQualityStatistics
{
  public:
    MeshQualityStatistics(std::shared_ptr<Mesh> mesh, std::shared_ptr<opt::impl::PatchObjective> obj)
      : m_mesh(mesh)
      , m_obj(obj)
    {
      int nvals = 0;
      for (auto vert : mesh->get_vertices())
        if (vert)
        {
          m_activeVerts.emplace_back(mesh, vert);
          nvals++;
        }

      m_qualityValues.resize(nvals);
    }

    struct Statistics
    {
        double mean;
        double stddev;
        double min;
        double max;
    };

    Statistics compute_statistics()
    {
      int idx = 0;
      for (auto& activeVert : m_activeVerts)
      {
        m_obj->set_active_patch(activeVert);
        auto pt      = activeVert.get_current_vert()->get_point_orig(0);
        auto ptPrime = m_obj->compute_parameterization(pt);
        // Note: the objective computes (roughly) 1/quality, because it needs
        //       good quality to be the *minimum* value
        m_qualityValues[idx++] = 1.0 / m_obj->compute_quality(ptPrime);
      }

      compute_mean();
      compute_standard_deviation();
      compute_min_max();

      return {m_mean, m_stddev, m_min, m_max};
    }

    void write_to_file(const std::string& fname)
    {
      std::ofstream of;
      of.open(fname);

      for (auto val : m_qualityValues)
        of << val << std::endl;

      of.close();
    }

  private:
    void compute_mean()
    {
      double mean = 0;
      for (auto val : m_qualityValues)
        mean += val;

      m_mean = mean / m_qualityValues.size();
    }

    void compute_standard_deviation()
    {
      double stddev = 0;
      for (auto val : m_qualityValues)
        stddev = (val - m_mean) * (val - m_mean);

      m_stddev = std::sqrt(stddev / m_qualityValues.size());
    }

    void compute_min_max()
    {
      m_max = m_qualityValues[0], m_min = m_qualityValues[0];
      for (auto valI : m_qualityValues)
      {
        m_min = std::min(m_min, valI);
        m_max = std::max(m_max, valI);
      }
    }

    std::shared_ptr<Mesh> m_mesh;
    std::shared_ptr<opt::impl::PatchObjective> m_obj;
    std::vector<opt::impl::ActiveVertData> m_activeVerts;
    std::vector<double> m_qualityValues;
    double m_mean;
    double m_max;
    double m_min;
    double m_stddev;
};

} // namespace impl

} // namespace mesh
} // namespace middle_mesh
} // namespace stk
#endif