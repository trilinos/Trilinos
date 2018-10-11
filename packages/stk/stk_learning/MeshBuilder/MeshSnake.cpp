#include "MeshSnake.hpp"
#include <ctime>                        // for time, NULL
#include <iostream>
#include <vector>                       // for vector
#include "MeshBuilder/MeshBuilder.hpp"  // for MeshBuilder, HexMeshBuilder, etc
#include "stk_util/util/ReportHandler.hpp"  // for ThrowRequire

/*
 * MeshSnake
 */
//public
MeshSnake::MeshSnake(MeshBuilder* mesh)
:m_mesh(mesh), m_moveCounter(0), m_moveMax(5)
{
}
void MeshSnake::begin_snake()
{
    srand(time(NULL));
    get_starting_direction();
    update_move_counter();
}
void MeshSnake::crawl(unsigned numSteps)
{
    m_mesh->begin_modification();
    for (unsigned stepNum = 0; stepNum < numSteps; stepNum++)
    {
        m_mesh->write_mesh();
        move();
        std::cout << stepNum << std::endl;
    }
    m_mesh->end_modification();
}

//private
void MeshSnake::update_move_counter()
{
    m_moveMax = rand()%10 + 2;
    m_moveCounter = 0;
}
void MeshSnake::move()
{
    if (facing_wall() || time_to_switch_directions())
    {
        get_new_direction();
        update_move_counter();
    }
    move_one();
    increment_move_counter();
}
bool MeshSnake::time_to_switch_directions()
{
    return m_moveCounter == m_moveMax;
}
void MeshSnake::increment_move_counter()
{
   m_moveCounter++;
}
/*
 * QuadMeshSnake
 */
//public
QuadMeshSnake::QuadMeshSnake(QuadMeshBuilder& mesh)
:MeshSnake(&mesh), m_quadMesh(mesh), m_xLower(1), m_xUpper(1), m_yLower(1), m_yUpper(1),
 m_xPos(1), m_yPos(1), m_dir(INVALID_DIR)
{
}
void QuadMeshSnake::set_x_bounds(unsigned xLower, unsigned xUpper)
{
    m_xLower = xLower;
    m_xUpper = xUpper;
}
void QuadMeshSnake::set_y_bounds(unsigned yLower, unsigned yUpper)
{
    m_yLower = yLower;
    m_yUpper = yUpper;
}
void QuadMeshSnake::set_x_pos(unsigned x)
{
    if (m_xLower <= x && x <= m_xUpper)
        m_xPos = x;
}
void QuadMeshSnake::set_y_pos(unsigned y)
{
    if (m_yLower <= y && y <= m_yUpper)
        m_yPos = y;
}

void QuadMeshSnake::set_pos(unsigned x, unsigned y)
{
    set_x_pos(x);
    set_y_pos(y);
}
//private
void QuadMeshSnake::get_starting_direction()
{
    switch(rand()%4)
    {
        case 0:
            m_dir = X_POSITIVE;
            break;
        case 1:
            m_dir = X_NEGATIVE;
            break;
        case 2:
            m_dir = Y_POSITIVE;
            break;
        case 3:
            m_dir = Y_NEGATIVE;
            break;
    }
}
bool QuadMeshSnake::wall_on_right()
{
   return m_xPos == m_xUpper;
}
bool QuadMeshSnake::wall_on_left()
{
   return m_xPos == m_xLower;
}
bool QuadMeshSnake::wall_on_top()
{
   return m_yPos == m_yUpper;
}
bool QuadMeshSnake::wall_on_bottom()
{
   return m_yPos == m_yLower;
}
bool QuadMeshSnake::facing_wall()
{
    bool result = false;
    switch (m_dir)
    {
        case X_POSITIVE:
            result = wall_on_right();
            break;
        case X_NEGATIVE:
            result = wall_on_left();
            break;
        case Y_POSITIVE:
            result = wall_on_top();
            break;
        case Y_NEGATIVE:
            result = wall_on_bottom();
            break;
        default:
            ThrowRequire(false);
    }
    return result;
}
void QuadMeshSnake::get_new_direction()
{
    switch (m_dir)
    {
        case X_POSITIVE:
        case X_NEGATIVE:
            get_new_y_direction();
            break;
        case Y_POSITIVE:
        case Y_NEGATIVE:
            get_new_x_direction();
            break;
        default:
            ThrowRequire(false);
    }
}
void QuadMeshSnake::get_new_y_direction()
{
    bool top = wall_on_top();
    bool bottom = wall_on_bottom();
    ThrowRequire(!(top && bottom));

    if (top)
        m_dir = Y_NEGATIVE;
    else if (bottom)
        m_dir = Y_POSITIVE;
    else
        m_dir = rand()%2 ? Y_NEGATIVE : Y_POSITIVE;

}
void QuadMeshSnake::get_new_x_direction()
{
    bool right = wall_on_right();
    bool left = wall_on_left();
    ThrowRequire(!(right && left));

    if (right)
        m_dir = X_NEGATIVE;
    else if (left)
        m_dir = X_POSITIVE;
    else
        m_dir = rand()%2 ? X_NEGATIVE : X_POSITIVE;

}
void QuadMeshSnake::move_one()
{
    m_quadMesh.remove_element(m_xPos, m_yPos);
    move_forward();
}
void QuadMeshSnake::move_forward()
{
    switch (m_dir)
    {
        case X_POSITIVE:
            m_xPos++;
            break;
        case X_NEGATIVE:
            m_xPos--;
            break;
        case Y_POSITIVE:
            m_yPos++;
            break;
        case Y_NEGATIVE:
            m_yPos--;
            break;
        default:
            ThrowRequire(false);
    }
}

/*
 * HexMeshSnake (yaaay)
 */
HexMeshSnake::HexMeshSnake(HexMeshBuilder& mesh)
:MeshSnake(&mesh), m_hexMesh(mesh), m_xLower(1), m_xUpper(1), m_yLower(1), m_yUpper(1),
m_zLower(1), m_zUpper(1), m_xPos(1), m_yPos(1), m_zPos(1), m_dir(INVALID_DIR)
{

}
void HexMeshSnake::set_x_bounds(unsigned xLower, unsigned xUpper)
{
    m_xLower = xLower;
    m_xUpper= xUpper;
}
void HexMeshSnake::set_y_bounds(unsigned yLower, unsigned yUpper)
{
    m_yLower = yLower;
    m_yUpper= yUpper;
}
void HexMeshSnake::set_z_bounds(unsigned zLower, unsigned zUpper)
{
    m_zLower = zLower;
    m_zUpper = zUpper;
}
void HexMeshSnake::set_x_pos(unsigned x)
{
    if (m_xLower <= x && x <= m_xUpper)
        m_xPos = x;
}
void HexMeshSnake::set_y_pos(unsigned y)
{
    if (m_yLower <= y && y <= m_yUpper)
        m_yPos = y;
}
void HexMeshSnake::set_z_pos(unsigned z)
{
    if (m_zLower <= z && z <= m_zUpper)
        m_zPos = z;
}
void HexMeshSnake::set_pos(unsigned x, unsigned y, unsigned z)
{
    set_x_pos(x);
    set_y_pos(y);
    set_z_pos(z);
}
//private
void HexMeshSnake::get_starting_direction()
{
    switch(rand()%6)
    {
        case 0:
            m_dir = X_POSITIVE;
            break;
        case 1:
            m_dir = X_NEGATIVE;
            break;
        case 2:
            m_dir = Y_POSITIVE;
            break;
        case 3:
            m_dir = Y_NEGATIVE;
            break;
        case 4:
            m_dir = Z_POSITIVE;
            break;
        case 5:
            m_dir = Z_NEGATIVE;
            break;
        default:
            ThrowRequire(false);
    }
}
bool HexMeshSnake::wall_on_right()
{
    return m_xPos == m_xUpper;
}
bool HexMeshSnake::wall_on_left()
{
    return m_xPos == m_xLower;
}
bool HexMeshSnake::wall_on_top()
{
    return m_yPos == m_yUpper;
}
bool HexMeshSnake::wall_on_bottom()
{
    return m_yPos == m_yLower;
}
bool HexMeshSnake::wall_on_front()
{
    return m_zPos == m_zUpper;
}
bool HexMeshSnake::wall_on_back()
{
    return m_zPos == m_zLower;
}
bool HexMeshSnake::facing_wall()
{
   bool result = false;
    switch (m_dir)
    {
        case X_POSITIVE:
            result = wall_on_right();
            break;
        case X_NEGATIVE:
            result = wall_on_left();
            break;
        case Y_POSITIVE:
            result = wall_on_top();
            break;
        case Y_NEGATIVE:
            result = wall_on_bottom();
            break;
        case Z_POSITIVE:
            result = wall_on_front();
            break;
        case Z_NEGATIVE:
            result = wall_on_back();
            break;
        default:
            ThrowRequire(false);
    }
    return result;
}
void HexMeshSnake::get_new_direction()
{
    switch(m_dir)
    {
        case X_POSITIVE:
        case X_NEGATIVE:
            get_new_yz_direction();
            break;
        case Y_POSITIVE:
        case Y_NEGATIVE:
            get_new_xz_direction();
            break;
        case Z_POSITIVE:
        case Z_NEGATIVE:
            get_new_xy_direction();
            break;
        default:
            ThrowRequire(false);
    }
}
void HexMeshSnake::get_new_xy_direction()
{
    std::vector<direction> availableDir;

    if (!wall_on_top())
        availableDir.push_back(Y_POSITIVE);
    if (!wall_on_bottom())
        availableDir.push_back(Y_NEGATIVE);
    if (!wall_on_right())
        availableDir.push_back(X_POSITIVE);
    if (!wall_on_left())
        availableDir.push_back(X_NEGATIVE);

    m_dir = availableDir[rand()%availableDir.size()];
}
void HexMeshSnake::get_new_xz_direction()
{
    std::vector<direction> availableDir;

    if (!wall_on_right())
        availableDir.push_back(X_POSITIVE);
    if (!wall_on_left())
        availableDir.push_back(X_NEGATIVE);
    if (!wall_on_front())
        availableDir.push_back(Z_POSITIVE);
    if (!wall_on_back())
        availableDir.push_back(Z_NEGATIVE);

    m_dir = availableDir[rand()%availableDir.size()];

}
void HexMeshSnake::get_new_yz_direction()
{
    std::vector<direction> availableDir;

    if (!wall_on_top())
        availableDir.push_back(Y_POSITIVE);
    if (!wall_on_bottom())
        availableDir.push_back(Y_NEGATIVE);
    if (!wall_on_front())
        availableDir.push_back(Z_POSITIVE);
    if (!wall_on_back())
        availableDir.push_back(Z_NEGATIVE);

    m_dir = availableDir[rand()%availableDir.size()];

}
void HexMeshSnake::move_one()
{
    m_hexMesh.remove_element(m_xPos, m_yPos, m_zPos);
    move_forward();
}
void HexMeshSnake::move_forward()
{
    switch (m_dir)
    {
        case X_POSITIVE:
            m_xPos++;
            break;
        case X_NEGATIVE:
            m_xPos--;
            break;
        case Y_POSITIVE:
            m_yPos++;
            break;
        case Y_NEGATIVE:
            m_yPos--;
            break;
        case Z_POSITIVE:
            m_zPos++;
            break;
        case Z_NEGATIVE:
            m_zPos--;
            break;
        default:
            ThrowRequire(false);
    }
}
