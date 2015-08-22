#ifndef COORDINATESETS_H
#define COORDINATESETS_H

#include <stk_mesh/base/Types.hpp>

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
        elemId = (y * (y+1))/2 + (y-1)*(x-1) + (x)*(x-1)/2;
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
       elemId = calculate_element_id(x, y, z);
       generate_node_ids();
    }
    unsigned x;
    unsigned y;
    unsigned z;
   stk::mesh::EntityId elemId;
   stk::mesh::EntityIdVector nodeIds;

private:
    stk::mesh::EntityId calculate_element_id(unsigned x, unsigned y, unsigned z)
    {
        unsigned xPrime = x+y-2;
        return (z*(z+1)*(z+2)+xPrime*(xPrime*xPrime+3*z*xPrime+3*z*z+6*z-7))/6+(y-1);
    }
    void generate_node_ids()
    {
        nodeIds.resize(8);
        stk::mesh::EntityId one = elemId;
        stk::mesh::EntityId two = calculate_element_id(x+1, y, z);
        stk::mesh::EntityId three = calculate_element_id(x+1, y+1, z);
        stk::mesh::EntityId four = calculate_element_id(x, y+1, z);
        stk::mesh::EntityId five = calculate_element_id(x, y, z+1);
        stk::mesh::EntityId six = calculate_element_id(x+1, y, z+1);
        stk::mesh::EntityId seven = calculate_element_id(x+1, y+1, z+1);
        stk::mesh::EntityId eight = calculate_element_id(x, y+1, z+1);
        nodeIds = {one, two, three, four, five, six, seven, eight};
    }
};

#endif // COORDINATESETS_H
