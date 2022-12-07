#ifndef KRINO_KRINO_UNIT_TESTS_AKRI_MESHSPECS_HPP_
#define KRINO_KRINO_UNIT_TESTS_AKRI_MESHSPECS_HPP_
#include <array>
#include <vector>

#include <stk_math/StkVector.hpp>

namespace krino {

struct RegularTri
{
    RegularTri() = default;
    static constexpr int DIM = 2;
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
    static constexpr int DIM = 2;
    std::vector<stk::math::Vector2d> nodeLocs
    {{
        { 0.000,  0.000 },
        { 1.000,  0.000 },
        { 0.500,  std::sqrt(3.)/2. },

        { 0.500,  0.000 },
        { 0.750,  std::sqrt(3.)/4. },
        { 0.250,  std::sqrt(3.)/4. }
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
    static constexpr int DIM = 2;
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
    static constexpr int DIM = 2;
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
    static constexpr int DIM = 2;
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

struct RegularTet
{
    RegularTet() = default;
    static constexpr int DIM = 3;
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
    static constexpr int DIM = 3;
    std::vector<stk::math::Vector3d> nodeLocs
    {{
        { 0.5,  0.0, -0.5/std::sqrt(2.) },
        {-0.5,  0.0, -0.5/std::sqrt(2.) },
        { 0.0, -0.5,  0.5/std::sqrt(2.) },
        { 0.0,  0.5,  0.5/std::sqrt(2.) },

        { 0.0,  0.0, -std::sqrt(2.) },
        {-0.25,-0.25, 0.0 },
        { 0.25,-0.25, 0.0 },
        { 0.25, 0.25, 0.0 },
        {-0.25, 0.25, 0.0 },
        { 0.0,  0.0,  std::sqrt(2.) },
    }};

    std::array<unsigned,4> TetConn0{{0,4,6,7}};
    std::array<unsigned,4> TetConn1{{1,5,4,8}};
    std::array<unsigned,4> TetConn2{{2,6,5,9}};
    std::array<unsigned,4> TetConn3{{3,8,7,9}};
    std::array<unsigned,4> TetConn4{{6,8,7,4}};
    std::array<unsigned,4> TetConn5{{5,8,6,4}};
    std::array<unsigned,4> TetConn6{{6,7,8,9}};
    std::array<unsigned,4> TetConn7{{5,6,8,9}};
    std::vector<std::array<unsigned, 4>> allElementConn{TetConn0,TetConn1,TetConn2,TetConn3,TetConn4,TetConn5,TetConn6,TetConn7};
};

struct RightTet
{
    RightTet() = default;
    static constexpr int DIM = 3;
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
    static constexpr int DIM = 3;
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
    static constexpr int DIM = 3;
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
    static constexpr int DIM = 3;
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
    static constexpr int DIM = 3;
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
    static constexpr int DIM = 2;
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

}



#endif /* KRINO_KRINO_UNIT_TESTS_AKRI_MESHSPECS_HPP_ */
