#ifndef COORDINATESETS_H
#define COORDINATESETS_H

#include <stk_mesh/base/Types.hpp>

inline stk::mesh::EntityId generate_two_dim_elem_id(unsigned x, unsigned y)
//returns a unique result for each coordinate in positive 2D space
{
    return (x*x + 2*x*y - 3*x + y*y - y + 2)/2;
}

inline stk::mesh::EntityId generate_three_dim_elem_id(unsigned x, unsigned y, unsigned z)
//returns a unique result for each coordinate in positive 3D space
{
    return (x*x*x + 3*x*x*y + 3*x*x*z - 6*x*x + 3*x*y*y + 6*x*y*z - 12*x*y + 3*x*z*z - 6*x*z +
            5*x + y*y*y + 3*y*y*z - 6*y*y + 3*y*z*z - 6*y*z + 11*y + z*z*z - 3*z*z + 2*z)/6;
}

struct ElemCoordPair
{
    ElemCoordPair(unsigned xCoord, unsigned yCoord)
    {
       x = xCoord;
       y = yCoord;
       generate_element_id();
       generate_node_ids();
    }
    unsigned x;
    unsigned y;
    stk::mesh::EntityId elemId;
    stk::mesh::EntityIdVector nodeIds;
private:
    void  generate_element_id()
    {
        elemId = generate_two_dim_elem_id(x, y);
    }
    void generate_node_ids()
    {
        nodeIds.resize(4);
        unsigned XYOffset = x + y;
        nodeIds = {elemId, elemId+ XYOffset-1, elemId+ 2*XYOffset, elemId+XYOffset};
    }
};

struct ElemCoordTriple
{
    ElemCoordTriple(unsigned xCoord, unsigned yCoord, unsigned zCoord)
    {
       x = xCoord;
       y = yCoord;
       z = zCoord;
       generate_element_id();
       generate_node_ids();
    }
    unsigned x;
    unsigned y;
    unsigned z;
    stk::mesh::EntityId elemId;
    stk::mesh::EntityIdVector nodeIds;
private:
    void generate_element_id()
    {
        elemId = generate_three_dim_elem_id(x, y, z);
    }
    void generate_node_ids()
    {
        nodeIds.resize(8);
        stk::mesh::EntityId one = elemId;
        stk::mesh::EntityId two = generate_three_dim_elem_id(x+1, y, z);
        stk::mesh::EntityId three = generate_three_dim_elem_id(x+1, y+1, z);
        stk::mesh::EntityId four = generate_three_dim_elem_id(x, y+1, z);
        stk::mesh::EntityId five = generate_three_dim_elem_id(x, y, z+1);
        stk::mesh::EntityId six = generate_three_dim_elem_id(x+1, y, z+1);
        stk::mesh::EntityId seven = generate_three_dim_elem_id(x+1, y+1, z+1);
        stk::mesh::EntityId eight = generate_three_dim_elem_id(x, y+1, z+1);
        nodeIds = {one, two, three, four, five, six, seven, eight};
    }
};

#endif // COORDINATESETS_H
