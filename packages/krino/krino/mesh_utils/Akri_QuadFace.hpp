/*
 * Akri_QuadFace.hpp
 *
 *  Created on: Apr 25, 2024
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_MESH_UTILS_AKRI_QUADFACE_HPP_
#define KRINO_KRINO_MESH_UTILS_AKRI_QUADFACE_HPP_

#include <array>
#include <bitset>
#include <unordered_map>
#include <vector>

#include <stk_mesh/base/Entity.hpp>

namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { struct topology; }

namespace krino {
struct QuadFace
{
    typedef std::bitset<128> quad_face_value_type;

    static constexpr quad_face_value_type InvalidQuadFace = quad_face_value_type();

    quad_face_value_type mValue;

    QuadFace() : mValue(InvalidQuadFace) {}

    explicit QuadFace(quad_face_value_type value) : mValue(value) {}

    QuadFace operator=(quad_face_value_type val) { mValue = val; return *this;}

    quad_face_value_type value() const { return mValue; }

    bool is_valid() const { return mValue.any(); }

    bool operator==(QuadFace entity) const { return mValue == entity.mValue; }
    bool operator==(quad_face_value_type val) const { return mValue == val; }
    bool operator!=(QuadFace entity) const { return mValue != entity.mValue; }
    bool operator!=(quad_face_value_type val) const { return mValue != val; }
    bool operator<(QuadFace entity) const
    {
      for (int i = 128-1; i >= 0; i--)
          if (mValue[i] ^ entity.mValue[i]) return entity.mValue[i];
      return false;
    }
};

QuadFace quad_face_from_unordered_nodes(const stk::mesh::BulkData & mesh, const std::array<stk::mesh::Entity,4> & faceNodes);
QuadFace quad_face_from_unordered_nodes(const stk::mesh::BulkData & mesh, stk::mesh::Entity faceNode0, stk::mesh::Entity faceNode1, stk::mesh::Entity faceNode2, stk::mesh::Entity faceNode3);

std::array<stk::mesh::Entity,4> get_quad_face_nodes_sorted_by_id(const QuadFace quadFace);

void append_entity_quad_faces(const stk::mesh::BulkData & mesh, const stk::topology entityTopology, const stk::mesh::Entity entity, std::vector<QuadFace> & entityQuadFaces);

void append_entity_quad_faces(const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity, std::vector<QuadFace> & entityQuadFaces);
void fill_entity_quad_faces(const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity, std::vector<QuadFace> & entityQuadFaces);

}

template<>
struct std::hash<krino::QuadFace>
{
    std::size_t operator()(const krino::QuadFace & quadFace) const noexcept
    {
        return std::hash<krino::QuadFace::quad_face_value_type>{}(quadFace.value());
    }
};



#endif /* KRINO_KRINO_MESH_UTILS_AKRI_QUADFACE_HPP_ */
