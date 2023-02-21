#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_MESHSPECS_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_MESHSPECS_HPP_
#include <array>
#include <vector>

#include <stk_math/StkVector.hpp>
#include <stk_topology/topology_decl.hpp>

namespace krino {

struct RegularTri
{
    RegularTri() = default;
    static constexpr stk::topology::topology_t TOPOLOGY = stk::topology::TRIANGLE_3_2D;
    std::vector<stk::math::Vector2d> nodeLocs
    {{
        {-0.500,  0.000 },
        { 0.500,  0.000 },
        { 0.000,  std::sqrt(3.)/2. },
    }};

    std::array<unsigned,3> TriConn{{0, 1, 2}};
    std::vector<std::array<unsigned, 3>> allElementConn{TriConn};
};

struct UMRRegularTri
{
    UMRRegularTri() = default;
    static constexpr stk::topology::topology_t TOPOLOGY = stk::topology::TRIANGLE_3_2D;
    std::vector<stk::math::Vector2d> nodeLocs
    {{
        {-0.500,  0.000 },
        { 0.500,  0.000 },
        { 0.000,  std::sqrt(3.)/2. },

        { 0.000,  0.000 },
        { 0.250,  std::sqrt(3.)/4. },
        {-0.250,  std::sqrt(3.)/4. }
    }};

    std::array<unsigned,3> TriConn0{{3,4,5}};
    std::array<unsigned,3> TriConn1{{0,3,5}};
    std::array<unsigned,3> TriConn2{{1,4,3}};
    std::array<unsigned,3> TriConn3{{2,5,4}};
    std::vector<std::array<unsigned, 3>> allElementConn{TriConn0,TriConn1,TriConn2,TriConn3};
};

struct RightTriSurroundedByEdgeTris
{
    RightTriSurroundedByEdgeTris() = default;
    static constexpr stk::topology::topology_t TOPOLOGY = stk::topology::TRIANGLE_3_2D;
    std::vector<stk::math::Vector2d> nodeLocs
    {{
        { -std::sqrt(3.)/2.,  0.500 },
        { 0.500,  -std::sqrt(3.)/2. },
        { 0.5+std::sqrt(3.)/2., 0.5+std::sqrt(3.)/2. },
        { 0.000,  0.000 },
        { 1.000,  0.000 },
        { 0.000,  1.000 }
    }};

    std::array<unsigned,3> TriConn1{{0,3,5}};
    std::array<unsigned,3> TriConn2{{1,4,3}};
    std::array<unsigned,3> TriConn3{{2,5,4}};
    std::array<unsigned,3> TriConn4{{3,4,5}};
    std::vector<std::array<unsigned, 3>> allElementConn{TriConn1,TriConn2,TriConn3,TriConn4};
};

struct Tri306090
{
    Tri306090() = default;
    static constexpr stk::topology::topology_t TOPOLOGY = stk::topology::TRIANGLE_3_2D;
    std::vector<stk::math::Vector2d> nodeLocs
    {{
        { 0.000,  0.000 },
        { 0.500,  0.000 },
        { 0.000,  std::sqrt(3.)/2. },
    }};

    std::array<unsigned,3> TriConn{{0, 1, 2}};
    std::vector<std::array<unsigned, 3>> allElementConn{TriConn};
};

struct TwoTri306090
{
    TwoTri306090() = default;
    static constexpr stk::topology::topology_t TOPOLOGY = stk::topology::TRIANGLE_3_2D;
    std::vector<stk::math::Vector2d> nodeLocs
    {{
        { 0.000,  0.000 },
        { 0.500,  0.000 },
        { 0.000,  std::sqrt(3.)/2. },
        {-0.500,  0.000 }
    }};

    std::array<unsigned,3> Tri1Conn{{0, 1, 2}};
    std::array<unsigned,3> Tri2Conn{{0, 2, 3}};
    std::vector<std::array<unsigned, 3>> allElementConn{Tri1Conn, Tri2Conn};
};

struct QuadSplit4Tri
{
    QuadSplit4Tri() = default;
    static constexpr stk::topology::topology_t TOPOLOGY = stk::topology::TRIANGLE_3_2D;
    std::vector<stk::math::Vector2d> nodeLocs
    {{
        { -0.500,  -0.500 },
        { 0.500,  -0.500 },
        { 0.500,  0.500 },
        { -0.500,  0.500 },
        { 0,  0 },
    }};

    std::array<unsigned,3> Tri1Conn{{0, 1, 4}};
    std::array<unsigned,3> Tri2Conn{{1, 2, 4}};
    std::array<unsigned,3> Tri3Conn{{2, 3, 4}};
    std::array<unsigned,3> Tri4Conn{{3, 0, 4}};
    std::vector<std::array<unsigned, 3>> allElementConn{Tri1Conn, Tri2Conn, Tri3Conn, Tri4Conn };
};


struct RegularTet
{
    RegularTet() = default;
    static constexpr stk::topology::topology_t TOPOLOGY = stk::topology::TETRAHEDRON_4;
    std::vector<stk::math::Vector3d> nodeLocs
    {{
        { 0.5,  0.0, -0.5/std::sqrt(2.) },
        {-0.5,  0.0, -0.5/std::sqrt(2.) },
        { 0.0, -0.5,  0.5/std::sqrt(2.) },
        { 0.0,  0.5,  0.5/std::sqrt(2.) },
    }};

    std::array<unsigned,4> TetConn{{0, 1, 2, 3}};
    std::vector<std::array<unsigned, 4>> allElementConn{TetConn};
};

struct UMRRegularTet
{
    UMRRegularTet() = default;
    static constexpr stk::topology::topology_t TOPOLOGY = stk::topology::TETRAHEDRON_4;
    std::vector<stk::math::Vector3d> nodeLocs
    {{
        { 0.5,  0.0, -0.5/std::sqrt(2.) },
        {-0.5,  0.0, -0.5/std::sqrt(2.) },
        { 0.0, -0.5,  0.5/std::sqrt(2.) },
        { 0.0,  0.5,  0.5/std::sqrt(2.) },

        { 0.0,  0.0,  0.5/std::sqrt(2.) },
        {-0.25, 0.25, 0.0 },
        { 0.25, 0.25, 0.0 },
        { 0.25,-0.25, 0.0 },
        {-0.25,-0.25, 0.0 },
        { 0.0,  0.0, -0.5/std::sqrt(2.) },
    }};

    std::array<unsigned,4> TetConn0{{0,9,7,6}};
    std::array<unsigned,4> TetConn1{{9,1,8,5}};
    std::array<unsigned,4> TetConn2{{7,8,2,4}};
    std::array<unsigned,4> TetConn3{{6,5,4,3}};
    std::array<unsigned,4> TetConn4{{7,6,9,5}};
    std::array<unsigned,4> TetConn5{{7,4,6,5}};
    std::array<unsigned,4> TetConn6{{7,8,4,5}};
    std::array<unsigned,4> TetConn7{{7,9,8,5}};
    std::vector<std::array<unsigned, 4>> allElementConn{TetConn0,TetConn1,TetConn2,TetConn3,TetConn4,TetConn5,TetConn6,TetConn7};
};

struct RightTet
{
    RightTet() = default;
    static constexpr stk::topology::topology_t TOPOLOGY = stk::topology::TETRAHEDRON_4;
    std::vector<stk::math::Vector3d> nodeLocs
    {{
        { 0.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 },
    }};

    std::array<unsigned,4> TetConn{{0, 1, 2, 3}};
    std::vector<std::array<unsigned, 4>> allElementConn{TetConn};
};

struct RightTetSurroundedByFaceTets
{
    RightTetSurroundedByFaceTets() = default;
    static constexpr stk::topology::topology_t TOPOLOGY = stk::topology::TETRAHEDRON_4;
    std::vector<stk::math::Vector3d> nodeLocs
    {{
        { 0.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 },
        { -0.50, 0.333, 0.333 },
        { 0.333, -0.50, 0.333 },
        { 0.333, 0.333, -0.50 },
        { 0.666, 0.666, 0.666 },
    }};

    std::array<unsigned,4> TetConn0{{0, 1, 2, 3}};
    std::array<unsigned,4> TetConn1{{0, 3, 2, 4}};
    std::array<unsigned,4> TetConn2{{0, 1, 3, 5}};
    std::array<unsigned,4> TetConn3{{0, 2, 1, 6}};
    std::array<unsigned,4> TetConn4{{1, 2, 3, 7}};

    std::vector<std::array<unsigned, 4>> allElementConn{TetConn0,TetConn1,TetConn2,TetConn3,TetConn4};
};

struct RightTetSurroundedByEdgeTets
{
    RightTetSurroundedByEdgeTets() = default;
    static constexpr stk::topology::topology_t TOPOLOGY = stk::topology::TETRAHEDRON_4;
    std::vector<stk::math::Vector3d> nodeLocs
    {{
        { 0.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 },

        { 1.0,-1.0, 0.0},
        { 1.0, 0.0,-1.0},
        { 1.0, 1.0, 0.0},
        { 1.0, 1.0,-1.0},
        {-1.0, 1.0, 0.0},
        { 0.0, 1.0,-1.0},
        {-1.0, 0.0, 1.0},
        { 0.0,-1.0, 1.0},
        { 1.0,-1.0, 1.0},
        { 1.0, 0.0, 1.0},
        { 0.0, 1.0, 1.0},
        {-1.0, 1.0, 1.0},
    }};

    std::array<unsigned,4> TetConn0{{0, 1, 2, 3}};
    std::array<unsigned,4> TetConn1{{0, 1, 4, 5}};
    std::array<unsigned,4> TetConn2{{1, 2, 6, 7}};
    std::array<unsigned,4> TetConn3{{2, 0, 8, 9}};
    std::array<unsigned,4> TetConn4{{0, 3,10,11}};
    std::array<unsigned,4> TetConn5{{1, 3,12,13}};
    std::array<unsigned,4> TetConn6{{2, 3,14,15}};

    std::vector<std::array<unsigned, 4>> allElementConn{TetConn0,TetConn1,TetConn2,TetConn3,TetConn4,TetConn5,TetConn6};
};

struct FourRightTets
{
    FourRightTets() = default;
    static constexpr stk::topology::topology_t TOPOLOGY = stk::topology::TETRAHEDRON_4;
    std::vector<stk::math::Vector3d> nodeLocs
    {{
        { 0.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        {-1.0, 0.0, 0.0 },
        { 0.0,-1.0, 0.0 },
        { 0.0, 0.0, 1.0 },
    }};

    std::array<unsigned,4> Tet1Conn{{0, 1, 2, 5}};
    std::array<unsigned,4> Tet2Conn{{0, 2, 3, 5}};
    std::array<unsigned,4> Tet3Conn{{0, 3, 4, 5}};
    std::array<unsigned,4> Tet4Conn{{0, 4, 1, 5}};
    std::vector<std::array<unsigned, 4>> allElementConn{Tet1Conn, Tet2Conn, Tet3Conn, Tet4Conn};
};

struct TwoRightTets
{
    TwoRightTets() = default;
    static constexpr stk::topology::topology_t TOPOLOGY = stk::topology::TETRAHEDRON_4;
    std::vector<stk::math::Vector3d> nodeLocs
    {{
        { 0.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        {-1.0, 0.0, 0.0 },
        { 0.0, 0.0, 1.0 },
    }};

    std::array<unsigned,4> Tet1Conn{{0, 1, 2, 4}};
    std::array<unsigned,4> Tet2Conn{{0, 2, 3, 4}};
    std::vector<std::array<unsigned, 4>> allElementConn{Tet1Conn, Tet2Conn};
};

struct TwoRightTris
{
    TwoRightTris() = default;
    static constexpr stk::topology::topology_t TOPOLOGY = stk::topology::TRIANGLE_3_2D;
    std::vector<stk::math::Vector2d> nodeLocs
    {{
        { 0.0, 0.0 },
        { 1.0, 0.0 },
        { 0.0, 1.0 },
        {-1.0, 0.0 },
    }};

    std::array<unsigned,3> Tri1Conn{{0, 1, 2}};
    std::array<unsigned,3> Tri2Conn{{0, 2, 3}};
    std::vector<std::array<unsigned, 3>> allElementConn{Tri1Conn, Tri2Conn};
};

struct TwoQuads
{
    TwoQuads() = default;
    static constexpr stk::topology::topology_t TOPOLOGY = stk::topology::QUADRILATERAL_4_2D;
    std::vector<stk::math::Vector2d> nodeLocs
    {{
        { 0.0, 0.0 },
        { 1.0, 0.0 },
        { 1.0, 1.0 },
        { 0.0, 1.0 },
        {-1.0, 1.0 },
        {-1.0, 0.0 },
    }};

    std::array<unsigned,4> Quad1Conn{{0, 1, 2, 3}};
    std::array<unsigned,4> Quad2Conn{{0, 3, 4, 5}};
    std::vector<std::array<unsigned, 4>> allElementConn{Quad1Conn, Quad2Conn};
};

}



#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_MESHSPECS_HPP_ */
