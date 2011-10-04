#ifndef _ZOLTAN2_TEMPLATE_MACROS_HPP
#define _ZOLTAN2_TEMPLATE_MACROS_HPP

#include <Kokkos_DefaultNode.hpp>

////////////////////////////////////////////////////////////////////////
#define HELLO std::cout << "Hello from " << __func__ << std::endl

#define Z2CLASS_TEMPLATE \
        typename Scalar=float, \
        typename LNO=int, typename GNO=int, \
        typename LID=LNO, typename GID=GNO, \
        typename Node=Kokkos::DefaultNode::DefaultNodeType

#define Z2FN_TEMPLATE \
        typename Scalar, \
        typename LNO, typename GNO, \
        typename LID, typename GID, \
        typename Node

#define Z2PARAM_TEMPLATE \
        Scalar, LNO, GNO, LID, GID, Node

// TODO : Need to add sparse ops ??
#define CONSISTENT_TRILINOS_CLASS_TEMPLATE_LINE \
        template <typename Scalar=float, \
                  typename LNO=int, typename GNO=int, \
                  typename Node=Kokkos::DefaultNode::DefaultNodeType>

#define CONSISTENT_TRILINOS_TEMPLATE_PARAMS \
        Scalar, LNO, GNO, Node
#endif
