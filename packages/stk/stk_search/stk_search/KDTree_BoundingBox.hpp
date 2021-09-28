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

#ifndef STKSEARCH_BOUNDING_BOX_H_
#define STKSEARCH_BOUNDING_BOX_H_

// #######################  Start Clang Header Tool Managed Headers ########################
#include <stk_search/BoundingBox.hpp>
#include <Kokkos_Core.hpp>
#include <iostream>                     // for operator<<, basic_ostream, etc
// #######################   End Clang Header Tool Managed Headers  ########################

namespace stk {
  namespace search {

    /**
     *  Bounding box of an object plus an index to that object.  This object is used to create bounding box arrays for
     *  things like processors, faces, or nodes.
     */
    template <typename BoxType>
      class ObjectBoundingBox_T  {
    public:

      typedef typename BoxType::value_type coordinate_t;


      KOKKOS_DEFAULTED_FUNCTION ObjectBoundingBox_T() = default;                ///< Default empty box and invalid object id
      KOKKOS_FORCEINLINE_FUNCTION ObjectBoundingBox_T(const BoxType& box,         ///< Explicitly defined constructor
                                                      const int obj_num_) :
      m_box(box), obj_num(obj_num_) {}             
      KOKKOS_FORCEINLINE_FUNCTION ObjectBoundingBox_T(const ObjectBoundingBox_T &box);       ///< Explicit copy constructor
      KOKKOS_DEFAULTED_FUNCTION ObjectBoundingBox_T(ObjectBoundingBox_T&&) = default;        ///< Default Move constructor
      KOKKOS_DEFAULTED_FUNCTION ~ObjectBoundingBox_T() = default;                          ///< Destructor

      KOKKOS_FORCEINLINE_FUNCTION ObjectBoundingBox_T& operator = (const ObjectBoundingBox_T&);    ///< Standard assignment (allows a = b = c)
      KOKKOS_DEFAULTED_FUNCTION ObjectBoundingBox_T& operator = (ObjectBoundingBox_T&&) = default; ///< Default move assignment
      ///
      ///  Explicity set or extract the object index for this bounding box
      ///
      KOKKOS_FORCEINLINE_FUNCTION void set_object_number(const int &obj_num_)       {obj_num = obj_num_;}
      KOKKOS_FORCEINLINE_FUNCTION int  get_object_number()                    const {return obj_num;    }
      KOKKOS_FORCEINLINE_FUNCTION const BoxType& GetBox() const        {return m_box;}
      KOKKOS_FORCEINLINE_FUNCTION void SetBox(const BoxType& box) {m_box.set_box(box.get_x_min(), box.get_y_min(), box.get_z_min(), box.get_x_max(), box.get_y_max(), box.get_z_max());}
      KOKKOS_FORCEINLINE_FUNCTION void AddBox(const BoxType& box) {stk::search::add_to_box(m_box, box);}
      KOKKOS_FORCEINLINE_FUNCTION void Reset ()                              {m_box = BoxType();}
      KOKKOS_FORCEINLINE_FUNCTION BoxType& GetBox() {return m_box;}

    private:
      BoxType m_box;
      int obj_num = -1;
    };

    template <typename BoxType>
      KOKKOS_FORCEINLINE_FUNCTION ObjectBoundingBox_T<BoxType>::ObjectBoundingBox_T(const ObjectBoundingBox_T<BoxType> &box) :
        m_box(box.m_box),
        obj_num(box.obj_num)
        {}

    template <typename BoxType>
      KOKKOS_FORCEINLINE_FUNCTION ObjectBoundingBox_T<BoxType>& ObjectBoundingBox_T<BoxType>::operator = (const ObjectBoundingBox_T<BoxType>& P) {
      obj_num = P.get_object_number();
      m_box = P.GetBox();
      return *this;
    }

    template <typename BoxType>
      KOKKOS_FORCEINLINE_FUNCTION std::ostream& operator<<(std::ostream &output, const ObjectBoundingBox_T<BoxType> &box) {
      output<<"Min corner "<<box.GetBox().get_x_min()<<" "<<box.GetBox().get_y_min()<<" "<<box.GetBox().get_z_min()<<std::endl;
      output<<"Max corner "<<box.GetBox().get_x_max()<<" "<<box.GetBox().get_y_max()<<" "<<box.GetBox().get_z_max()<<std::endl;
      output<<"object number "<<box.get_object_number()<<std::endl;
      return output;
    }

  }
}

#endif // STKSEARCH_BOUNDING_BOX_H_
