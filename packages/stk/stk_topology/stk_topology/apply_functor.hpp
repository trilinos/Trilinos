// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
#ifndef STKTOPOLOGY_APPLY_FUNCTOR_TCC
#define STKTOPOLOGY_APPLY_FUNCTOR_TCC

#include "stk_topology/topology_decl.hpp"
#include "stk_topology/types.hpp"

namespace stk {

//*****************************************************************************
// Converts a runtime topology to a compile-time topology_type<Topology>
// and calls the given functor on the compile-time topology
//*****************************************************************************
template <typename Functor>
struct topology::apply_host_functor
{
  typedef typename Functor::result_type result_type;

  apply_host_functor()
    : m_functor()
  {}

  apply_host_functor(Functor f)
    : m_functor(f)
  {}

  result_type operator()(topology_t t) const
  {
    switch(t)
    {
    case INVALID_TOPOLOGY:             return m_functor( topology_type< INVALID_TOPOLOGY            >() );
    case NODE:                         return m_functor( topology_type< NODE                        >() );
    case LINE_2:                       return m_functor( topology_type< LINE_2                      >() );
    case LINE_3:                       return m_functor( topology_type< LINE_3                      >() );
    case TRI_3:                        return m_functor( topology_type< TRI_3                       >() );
    case TRI_4:                        return m_functor( topology_type< TRI_4                       >() );
    case TRI_6:                        return m_functor( topology_type< TRI_6                       >() );
    case QUAD_4:                       return m_functor( topology_type< QUAD_4                      >() );
    case QUAD_6:                       return m_functor( topology_type< QUAD_6                      >() );
    case QUAD_8:                       return m_functor( topology_type< QUAD_8                      >() );
    case QUAD_9:                       return m_functor( topology_type< QUAD_9                      >() );
    case PARTICLE:                     return m_functor( topology_type< PARTICLE                    >() );
    case LINE_2_1D:                    return m_functor( topology_type< LINE_2_1D                   >() );
    case LINE_3_1D:                    return m_functor( topology_type< LINE_3_1D                   >() );
    case BEAM_2:                       return m_functor( topology_type< BEAM_2                      >() );
    case BEAM_3:                       return m_functor( topology_type< BEAM_3                      >() );
    case SHELL_LINE_2:                 return m_functor( topology_type< SHELL_LINE_2                >() );
    case SHELL_LINE_3:                 return m_functor( topology_type< SHELL_LINE_3                >() );
    case SHELL_SIDE_BEAM_2:            return m_functor( topology_type< SHELL_SIDE_BEAM_2           >() );
    case SHELL_SIDE_BEAM_3:            return m_functor( topology_type< SHELL_SIDE_BEAM_3           >() );
    case SPRING_2:                     return m_functor( topology_type< SPRING_2                    >() );
    case SPRING_3:                     return m_functor( topology_type< SPRING_3                    >() );
    case TRI_3_2D:                     return m_functor( topology_type< TRI_3_2D                    >() );
    case TRI_4_2D:                     return m_functor( topology_type< TRI_4_2D                    >() );
    case TRI_6_2D:                     return m_functor( topology_type< TRI_6_2D                    >() );
    case QUAD_4_2D:                    return m_functor( topology_type< QUAD_4_2D                   >() );
    case QUAD_8_2D:                    return m_functor( topology_type< QUAD_8_2D                   >() );
    case QUAD_9_2D:                    return m_functor( topology_type< QUAD_9_2D                   >() );
    case SHELL_TRI_3:                  return m_functor( topology_type< SHELL_TRI_3                 >() );
    case SHELL_TRI_4:                  return m_functor( topology_type< SHELL_TRI_4                 >() );
    case SHELL_TRI_6:                  return m_functor( topology_type< SHELL_TRI_6                 >() );
    case SHELL_TRI_3_ALL_FACE_SIDES:   return m_functor( topology_type< SHELL_TRI_3_ALL_FACE_SIDES  >() );
    case SHELL_TRI_4_ALL_FACE_SIDES:   return m_functor( topology_type< SHELL_TRI_4_ALL_FACE_SIDES  >() );
    case SHELL_TRI_6_ALL_FACE_SIDES:   return m_functor( topology_type< SHELL_TRI_6_ALL_FACE_SIDES  >() );
    case SHELL_QUAD_4:                 return m_functor( topology_type< SHELL_QUAD_4                >() );
    case SHELL_QUAD_8:                 return m_functor( topology_type< SHELL_QUAD_8                >() );
    case SHELL_QUAD_9:                 return m_functor( topology_type< SHELL_QUAD_9                >() );
    case SHELL_QUAD_4_ALL_FACE_SIDES:  return m_functor( topology_type< SHELL_QUAD_4_ALL_FACE_SIDES >() );
    case SHELL_QUAD_8_ALL_FACE_SIDES:  return m_functor( topology_type< SHELL_QUAD_8_ALL_FACE_SIDES >() );
    case SHELL_QUAD_9_ALL_FACE_SIDES:  return m_functor( topology_type< SHELL_QUAD_9_ALL_FACE_SIDES >() );
    case TET_4:                        return m_functor( topology_type< TET_4                       >() );
    case TET_8:                        return m_functor( topology_type< TET_8                       >() );
    case TET_10:                       return m_functor( topology_type< TET_10                      >() );
    case TET_11:                       return m_functor( topology_type< TET_11                      >() );
    case PYRAMID_5:                    return m_functor( topology_type< PYRAMID_5                   >() );
    case PYRAMID_13:                   return m_functor( topology_type< PYRAMID_13                  >() );
    case PYRAMID_14:                   return m_functor( topology_type< PYRAMID_14                  >() );
    case WEDGE_6:                      return m_functor( topology_type< WEDGE_6                     >() );
    case WEDGE_12:                     return m_functor( topology_type< WEDGE_12                    >() );
    case WEDGE_15:                     return m_functor( topology_type< WEDGE_15                    >() );
    case WEDGE_18:                     return m_functor( topology_type< WEDGE_18                    >() );
    case HEX_8:                        return m_functor( topology_type< HEX_8                       >() );
    case HEX_20:                       return m_functor( topology_type< HEX_20                      >() );
    case HEX_27:                       return m_functor( topology_type< HEX_27                      >() );
    default: break;
    }
    return m_functor( topology_type<INVALID_TOPOLOGY>() );
  }
 
  result_type operator()(topology_t t)
  {
    switch(t)
    {
    case INVALID_TOPOLOGY:             return m_functor( topology_type< INVALID_TOPOLOGY            >() );
    case NODE:                         return m_functor( topology_type< NODE                        >() );
    case LINE_2:                       return m_functor( topology_type< LINE_2                      >() );
    case LINE_3:                       return m_functor( topology_type< LINE_3                      >() );
    case TRI_3:                        return m_functor( topology_type< TRI_3                       >() );
    case TRI_4:                        return m_functor( topology_type< TRI_4                       >() );
    case TRI_6:                        return m_functor( topology_type< TRI_6                       >() );
    case QUAD_4:                       return m_functor( topology_type< QUAD_4                      >() );
    case QUAD_6:                       return m_functor( topology_type< QUAD_6                      >() );
    case QUAD_8:                       return m_functor( topology_type< QUAD_8                      >() );
    case QUAD_9:                       return m_functor( topology_type< QUAD_9                      >() );
    case PARTICLE:                     return m_functor( topology_type< PARTICLE                    >() );
    case LINE_2_1D:                    return m_functor( topology_type< LINE_2_1D                   >() );
    case LINE_3_1D:                    return m_functor( topology_type< LINE_3_1D                   >() );
    case BEAM_2:                       return m_functor( topology_type< BEAM_2                      >() );
    case BEAM_3:                       return m_functor( topology_type< BEAM_3                      >() );
    case SHELL_LINE_2:                 return m_functor( topology_type< SHELL_LINE_2                >() );
    case SHELL_LINE_3:                 return m_functor( topology_type< SHELL_LINE_3                >() );
    case SHELL_SIDE_BEAM_2:            return m_functor( topology_type< SHELL_SIDE_BEAM_2           >() );
    case SHELL_SIDE_BEAM_3:            return m_functor( topology_type< SHELL_SIDE_BEAM_3           >() );
    case SPRING_2:                     return m_functor( topology_type< SPRING_2                    >() );
    case SPRING_3:                     return m_functor( topology_type< SPRING_3                    >() );
    case TRI_3_2D:                     return m_functor( topology_type< TRI_3_2D                    >() );
    case TRI_4_2D:                     return m_functor( topology_type< TRI_4_2D                    >() );
    case TRI_6_2D:                     return m_functor( topology_type< TRI_6_2D                    >() );
    case QUAD_4_2D:                    return m_functor( topology_type< QUAD_4_2D                   >() );
    case QUAD_8_2D:                    return m_functor( topology_type< QUAD_8_2D                   >() );
    case QUAD_9_2D:                    return m_functor( topology_type< QUAD_9_2D                   >() );
    case SHELL_TRI_3:                  return m_functor( topology_type< SHELL_TRI_3                 >() );
    case SHELL_TRI_4:                  return m_functor( topology_type< SHELL_TRI_4                 >() );
    case SHELL_TRI_6:                  return m_functor( topology_type< SHELL_TRI_6                 >() );
    case SHELL_TRI_3_ALL_FACE_SIDES:   return m_functor( topology_type< SHELL_TRI_3_ALL_FACE_SIDES  >() );
    case SHELL_TRI_4_ALL_FACE_SIDES:   return m_functor( topology_type< SHELL_TRI_4_ALL_FACE_SIDES  >() );
    case SHELL_TRI_6_ALL_FACE_SIDES:   return m_functor( topology_type< SHELL_TRI_6_ALL_FACE_SIDES  >() );
    case SHELL_QUAD_4:                 return m_functor( topology_type< SHELL_QUAD_4                >() );
    case SHELL_QUAD_8:                 return m_functor( topology_type< SHELL_QUAD_8                >() );
    case SHELL_QUAD_9:                 return m_functor( topology_type< SHELL_QUAD_9                >() );
    case SHELL_QUAD_4_ALL_FACE_SIDES:  return m_functor( topology_type< SHELL_QUAD_4_ALL_FACE_SIDES >() );
    case SHELL_QUAD_8_ALL_FACE_SIDES:  return m_functor( topology_type< SHELL_QUAD_8_ALL_FACE_SIDES >() );
    case SHELL_QUAD_9_ALL_FACE_SIDES:  return m_functor( topology_type< SHELL_QUAD_9_ALL_FACE_SIDES >() );
    case TET_4:                        return m_functor( topology_type< TET_4                       >() );
    case TET_8:                        return m_functor( topology_type< TET_8                       >() );
    case TET_10:                       return m_functor( topology_type< TET_10                      >() );
    case TET_11:                       return m_functor( topology_type< TET_11                      >() );
    case PYRAMID_5:                    return m_functor( topology_type< PYRAMID_5                   >() );
    case PYRAMID_13:                   return m_functor( topology_type< PYRAMID_13                  >() );
    case PYRAMID_14:                   return m_functor( topology_type< PYRAMID_14                  >() );
    case WEDGE_6:                      return m_functor( topology_type< WEDGE_6                     >() );
    case WEDGE_12:                     return m_functor( topology_type< WEDGE_12                    >() );
    case WEDGE_15:                     return m_functor( topology_type< WEDGE_15                    >() );
    case WEDGE_18:                     return m_functor( topology_type< WEDGE_18                    >() );
    case HEX_8:                        return m_functor( topology_type< HEX_8                       >() );
    case HEX_20:                       return m_functor( topology_type< HEX_20                      >() );
    case HEX_27:                       return m_functor( topology_type< HEX_27                      >() );
    default: break;
    }
    return m_functor( topology_type<INVALID_TOPOLOGY>() );
  }

  Functor m_functor;
};

template <typename Functor>
struct topology::apply_functor
{
  typedef typename Functor::result_type result_type;

  STK_FUNCTION
  apply_functor()
    : m_functor()
  {}

  STK_FUNCTION
  apply_functor(Functor f)
    : m_functor(f)
  {}

  STK_FUNCTION
  result_type operator()(topology_t t) const
  {
    switch(t)
    {
    case INVALID_TOPOLOGY:             return m_functor( topology_type< INVALID_TOPOLOGY             >() );
    case NODE:                         return m_functor( topology_type< NODE                         >() );
    case LINE_2:                       return m_functor( topology_type< LINE_2                       >() );
    case LINE_3:                       return m_functor( topology_type< LINE_3                       >() );
    case TRI_3:                        return m_functor( topology_type< TRI_3                        >() );
    case TRI_4:                        return m_functor( topology_type< TRI_4                        >() );
    case TRI_6:                        return m_functor( topology_type< TRI_6                        >() );
    case QUAD_4:                       return m_functor( topology_type< QUAD_4                       >() );
    case QUAD_6:                       return m_functor( topology_type< QUAD_6                       >() );
    case QUAD_8:                       return m_functor( topology_type< QUAD_8                       >() );
    case QUAD_9:                       return m_functor( topology_type< QUAD_9                       >() );
    case PARTICLE:                     return m_functor( topology_type< PARTICLE                     >() );
    case LINE_2_1D:                    return m_functor( topology_type< LINE_2_1D                    >() );
    case LINE_3_1D:                    return m_functor( topology_type< LINE_3_1D                    >() );
    case BEAM_2:                       return m_functor( topology_type< BEAM_2                       >() );
    case BEAM_3:                       return m_functor( topology_type< BEAM_3                       >() );
    case SHELL_LINE_2:                 return m_functor( topology_type< SHELL_LINE_2                 >() );
    case SHELL_LINE_3:                 return m_functor( topology_type< SHELL_LINE_3                 >() );
    case SHELL_SIDE_BEAM_2:            return m_functor( topology_type< SHELL_SIDE_BEAM_2            >() );
    case SHELL_SIDE_BEAM_3:            return m_functor( topology_type< SHELL_SIDE_BEAM_3            >() );
    case SPRING_2:                     return m_functor( topology_type< SPRING_2                     >() );
    case SPRING_3:                     return m_functor( topology_type< SPRING_3                     >() );
    case TRI_3_2D:                     return m_functor( topology_type< TRI_3_2D                     >() );
    case TRI_4_2D:                     return m_functor( topology_type< TRI_4_2D                     >() );
    case TRI_6_2D:                     return m_functor( topology_type< TRI_6_2D                     >() );
    case QUAD_4_2D:                    return m_functor( topology_type< QUAD_4_2D                    >() );
    case QUAD_8_2D:                    return m_functor( topology_type< QUAD_8_2D                    >() );
    case QUAD_9_2D:                    return m_functor( topology_type< QUAD_9_2D                    >() );
    case SHELL_TRI_3:                  return m_functor( topology_type< SHELL_TRI_3                  >() );
    case SHELL_TRI_4:                  return m_functor( topology_type< SHELL_TRI_4                  >() );
    case SHELL_TRI_6:                  return m_functor( topology_type< SHELL_TRI_6                  >() );
    case SHELL_TRI_3_ALL_FACE_SIDES:   return m_functor( topology_type< SHELL_TRI_3_ALL_FACE_SIDES   >() );
    case SHELL_TRI_4_ALL_FACE_SIDES:   return m_functor( topology_type< SHELL_TRI_4_ALL_FACE_SIDES   >() );
    case SHELL_TRI_6_ALL_FACE_SIDES:   return m_functor( topology_type< SHELL_TRI_6_ALL_FACE_SIDES   >() );
    case SHELL_QUAD_4:                 return m_functor( topology_type< SHELL_QUAD_4                 >() );
    case SHELL_QUAD_8:                 return m_functor( topology_type< SHELL_QUAD_8                 >() );
    case SHELL_QUAD_9:                 return m_functor( topology_type< SHELL_QUAD_9                 >() );
    case SHELL_QUAD_4_ALL_FACE_SIDES:  return m_functor( topology_type< SHELL_QUAD_4_ALL_FACE_SIDES  >() );
    case SHELL_QUAD_8_ALL_FACE_SIDES:  return m_functor( topology_type< SHELL_QUAD_8_ALL_FACE_SIDES  >() );
    case SHELL_QUAD_9_ALL_FACE_SIDES:  return m_functor( topology_type< SHELL_QUAD_9_ALL_FACE_SIDES  >() );
    case TET_4:                        return m_functor( topology_type< TET_4                        >() );
    case TET_8:                        return m_functor( topology_type< TET_8                        >() );
    case TET_10:                       return m_functor( topology_type< TET_10                       >() );
    case TET_11:                       return m_functor( topology_type< TET_11                       >() );
    case PYRAMID_5:                    return m_functor( topology_type< PYRAMID_5                    >() );
    case PYRAMID_13:                   return m_functor( topology_type< PYRAMID_13                   >() );
    case PYRAMID_14:                   return m_functor( topology_type< PYRAMID_14                   >() );
    case WEDGE_6:                      return m_functor( topology_type< WEDGE_6                      >() );
    case WEDGE_12:                     return m_functor( topology_type< WEDGE_12                     >() );
    case WEDGE_15:                     return m_functor( topology_type< WEDGE_15                     >() );
    case WEDGE_18:                     return m_functor( topology_type< WEDGE_18                     >() );
    case HEX_8:                        return m_functor( topology_type< HEX_8                        >() );
    case HEX_20:                       return m_functor( topology_type< HEX_20                       >() );
    case HEX_27:                       return m_functor( topology_type< HEX_27                       >() );
    default: break;
    }
    return m_functor( topology_type<INVALID_TOPOLOGY>() );
  }
 
  STK_FUNCTION
  result_type operator()(topology_t t)
  {
    switch(t)
    {
    case INVALID_TOPOLOGY:             return m_functor( topology_type< INVALID_TOPOLOGY             >() );
    case NODE:                         return m_functor( topology_type< NODE                         >() );
    case LINE_2:                       return m_functor( topology_type< LINE_2                       >() );
    case LINE_3:                       return m_functor( topology_type< LINE_3                       >() );
    case TRI_3:                        return m_functor( topology_type< TRI_3                        >() );
    case TRI_4:                        return m_functor( topology_type< TRI_4                        >() );
    case TRI_6:                        return m_functor( topology_type< TRI_6                        >() );
    case QUAD_4:                       return m_functor( topology_type< QUAD_4                       >() );
    case QUAD_6:                       return m_functor( topology_type< QUAD_6                       >() );
    case QUAD_8:                       return m_functor( topology_type< QUAD_8                       >() );
    case QUAD_9:                       return m_functor( topology_type< QUAD_9                       >() );
    case PARTICLE:                     return m_functor( topology_type< PARTICLE                     >() );
    case LINE_2_1D:                    return m_functor( topology_type< LINE_2_1D                    >() );
    case LINE_3_1D:                    return m_functor( topology_type< LINE_3_1D                    >() );
    case BEAM_2:                       return m_functor( topology_type< BEAM_2                       >() );
    case BEAM_3:                       return m_functor( topology_type< BEAM_3                       >() );
    case SHELL_LINE_2:                 return m_functor( topology_type< SHELL_LINE_2                 >() );
    case SHELL_LINE_3:                 return m_functor( topology_type< SHELL_LINE_3                 >() );
    case SHELL_SIDE_BEAM_2:            return m_functor( topology_type< SHELL_SIDE_BEAM_2            >() );
    case SHELL_SIDE_BEAM_3:            return m_functor( topology_type< SHELL_SIDE_BEAM_3            >() );
    case SPRING_2:                     return m_functor( topology_type< SPRING_2                     >() );
    case SPRING_3:                     return m_functor( topology_type< SPRING_3                     >() );
    case TRI_3_2D:                     return m_functor( topology_type< TRI_3_2D                     >() );
    case TRI_4_2D:                     return m_functor( topology_type< TRI_4_2D                     >() );
    case TRI_6_2D:                     return m_functor( topology_type< TRI_6_2D                     >() );
    case QUAD_4_2D:                    return m_functor( topology_type< QUAD_4_2D                    >() );
    case QUAD_8_2D:                    return m_functor( topology_type< QUAD_8_2D                    >() );
    case QUAD_9_2D:                    return m_functor( topology_type< QUAD_9_2D                    >() );
    case SHELL_TRI_3:                  return m_functor( topology_type< SHELL_TRI_3                  >() );
    case SHELL_TRI_4:                  return m_functor( topology_type< SHELL_TRI_4                  >() );
    case SHELL_TRI_6:                  return m_functor( topology_type< SHELL_TRI_6                  >() );
    case SHELL_TRI_3_ALL_FACE_SIDES:   return m_functor( topology_type< SHELL_TRI_3_ALL_FACE_SIDES   >() );
    case SHELL_TRI_4_ALL_FACE_SIDES:   return m_functor( topology_type< SHELL_TRI_4_ALL_FACE_SIDES   >() );
    case SHELL_TRI_6_ALL_FACE_SIDES:   return m_functor( topology_type< SHELL_TRI_6_ALL_FACE_SIDES   >() );
    case SHELL_QUAD_4:                 return m_functor( topology_type< SHELL_QUAD_4                 >() );
    case SHELL_QUAD_8:                 return m_functor( topology_type< SHELL_QUAD_8                 >() );
    case SHELL_QUAD_9:                 return m_functor( topology_type< SHELL_QUAD_9                 >() );
    case SHELL_QUAD_4_ALL_FACE_SIDES:  return m_functor( topology_type< SHELL_QUAD_4_ALL_FACE_SIDES  >() );
    case SHELL_QUAD_8_ALL_FACE_SIDES:  return m_functor( topology_type< SHELL_QUAD_8_ALL_FACE_SIDES  >() );
    case SHELL_QUAD_9_ALL_FACE_SIDES:  return m_functor( topology_type< SHELL_QUAD_9_ALL_FACE_SIDES  >() );
    case TET_4:                        return m_functor( topology_type< TET_4                        >() );
    case TET_8:                        return m_functor( topology_type< TET_8                        >() );
    case TET_10:                       return m_functor( topology_type< TET_10                       >() );
    case TET_11:                       return m_functor( topology_type< TET_11                       >() );
    case PYRAMID_5:                    return m_functor( topology_type< PYRAMID_5                    >() );
    case PYRAMID_13:                   return m_functor( topology_type< PYRAMID_13                   >() );
    case PYRAMID_14:                   return m_functor( topology_type< PYRAMID_14                   >() );
    case WEDGE_6:                      return m_functor( topology_type< WEDGE_6                      >() );
    case WEDGE_12:                     return m_functor( topology_type< WEDGE_12                     >() );
    case WEDGE_15:                     return m_functor( topology_type< WEDGE_15                     >() );
    case WEDGE_18:                     return m_functor( topology_type< WEDGE_18                     >() );
    case HEX_8:                        return m_functor( topology_type< HEX_8                        >() );
    case HEX_20:                       return m_functor( topology_type< HEX_20                       >() );
    case HEX_27:                       return m_functor( topology_type< HEX_27                       >() );
    default: break;
    }
    return m_functor( topology_type<INVALID_TOPOLOGY>() );
  }

  Functor m_functor;
};

} //namespace stk

#endif //STKTOPOLOGY_APPLY_FUNCTOR_TCC
