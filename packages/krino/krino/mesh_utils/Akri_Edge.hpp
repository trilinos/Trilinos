/*
 * Akri_Edge.hpp
 *
 *  Created on: Sep 23, 2022
 *      Author: drnoble
 */

#ifndef KRINO_KRINO_KRINO_LIB_AKRI_EDGE_HPP_
#define KRINO_KRINO_KRINO_LIB_AKRI_EDGE_HPP_
#include <array>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>
namespace stk { namespace mesh { struct Entity; } }
namespace stk { namespace mesh { class BulkData; } }
namespace stk { namespace mesh { class Selector; } }
namespace stk { struct topology; }

namespace krino {
struct Edge
{
    enum Edge_t {
        InvalidEdge = 0ull,
    };

    typedef uint64_t edge_value_type;

    edge_value_type mValue;

    Edge() : mValue(InvalidEdge) {}

    explicit Edge(edge_value_type value) : mValue(value) {}

    Edge operator=(edge_value_type val) { mValue = val; return *this;}

    edge_value_type value() const { return mValue; }

    bool is_valid() const { return mValue > 0; }

    bool operator==(Edge entity) const { return mValue == entity.mValue; }
    bool operator==(edge_value_type val) const { return mValue == val; }
    bool operator!=(Edge entity) const { return mValue != entity.mValue; }
    bool operator!=(edge_value_type val) const { return mValue != val; }
    bool operator<(Edge entity) const { return mValue < entity.mValue; }
};

Edge edge_from_edge_nodes(const stk::mesh::BulkData & mesh, stk::mesh::Entity edgeNode0, stk::mesh::Entity edgeNode1);

std::array<stk::mesh::Entity,2> get_edge_nodes(const Edge edge);

void fill_edge_nodes(const Edge edge, std::vector<stk::mesh::Entity> & edgeNodes);

void fill_entity_edges(const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity, std::vector<Edge> & entityEdges);
void append_entity_edges(const stk::mesh::BulkData & mesh, const stk::topology entityTopology, const stk::mesh::Entity entity, std::vector<Edge> & entityEdges);
void append_entity_edges(const stk::mesh::BulkData & mesh, const stk::mesh::Entity entity, std::vector<Edge> & entityEdges);

int get_edge_parallel_owner_rank(const stk::mesh::BulkData & mesh, const Edge edge);

std::string debug_edge(const stk::mesh::BulkData & mesh, const Edge edge);

std::vector<Edge> get_edges_of_selected_elements(const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector);
std::vector<Edge> get_edges_of_elements(const stk::mesh::BulkData & mesh, const std::vector<stk::mesh::Entity> & elements);

}

template<>
struct std::hash<krino::Edge>
{
    std::size_t operator()(const krino::Edge & edge) const noexcept
    {
        return edge.value();
    }
};

#endif /* KRINO_KRINO_KRINO_LIB_AKRI_EDGE_HPP_ */
