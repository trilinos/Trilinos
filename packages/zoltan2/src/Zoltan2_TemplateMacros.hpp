#ifndef _ZOLTAN2_TEMPLATE_MACROS_HPP
#define _ZOLTAN2_TEMPLATE_MACROS_HPP

#include <Kokkos_DefaultNode.hpp>

////////////////////////////////////////////////////////////////////////
#define HELLO cout << "Hello from " << __func__ << endl

#define CONSISTENT_CLASS_TEMPLATE_LINE \
        template <typename Scalar=float, \
      typename LNO=int, typename GNO=int, typename LID=LNO, typename GID=GNO, \
      typename Node=Kokkos::DefaultNode::DefaultNodeType>

#define CONSISTENT_FN_TEMPLATE_LINE \
        template <typename Scalar, \
                  typename LNO, typename GNO, typename LID, typename GID, \
                  typename Node>

#define CONSISTENT_TEMPLATE_PARAMS \
        Scalar, LNO, GNO, LID, GID, Node

#endif
