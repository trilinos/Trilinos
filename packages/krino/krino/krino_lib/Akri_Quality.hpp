#ifndef KRINO_KRINO_KRINO_LIB_AKRI_QUALITY_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_QUALITY_HPP_

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Selector; } }

namespace krino
{
class QualityMetric;

double compute_mesh_quality(const stk::mesh::BulkData & mesh,
    const stk::mesh::Selector & elementSelector,
    const QualityMetric &qualityMetric);

}

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_QUALITY_HPP_ */
