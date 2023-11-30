#ifndef KRINO_KRINO_KRINO_LIB_AKRI_REFINEMENTREBALANCE_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_REFINEMENTREBALANCE_HPP_

namespace stk { namespace mesh { class BulkData; } }

namespace krino {

class Refinement;

bool rebalance_refined_mesh(const Refinement & refinement, stk::mesh::BulkData & mesh);

} // namespace krino
#endif /* KRINO_KRINO_KRINO_LIB_AKRI_REFINEMENTREBALANCE_HPP_ */
