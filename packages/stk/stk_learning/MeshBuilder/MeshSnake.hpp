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

#ifndef _MESHSNAKE_HPP_
#define _MESHSNAKE_HPP_

#include <stdlib.h>                     // for rand
class HexMeshBuilder;
class MeshBuilder;
class QuadMeshBuilder;


enum direction
{
    X_POSITIVE,
    X_NEGATIVE,
    Y_POSITIVE,
    Y_NEGATIVE,
    Z_POSITIVE,
    Z_NEGATIVE,
    INVALID_DIR
};
inline direction rand_2D_dir()
{
    direction dir;
    switch (rand()%4)
    {
        case 0:
            dir = X_POSITIVE;
            break;
        case 1:
            dir = X_NEGATIVE;
            break;
        case 2:
            dir = Y_POSITIVE;
            break;
        case 3:
            dir = Y_NEGATIVE;
            break;
    }
    return dir;
}
inline direction rand_3D_dir()
{
    direction dir;
    switch (rand()%6)
    {
        case 0:
            dir = X_POSITIVE;
            break;
        case 1:
            dir = X_NEGATIVE;
            break;
        case 2:
            dir = Y_POSITIVE;
            break;
        case 3:
            dir = Y_NEGATIVE;
            break;
        case 4:
            dir = Z_POSITIVE;
            break;
        case 5:
            dir = Z_NEGATIVE;
            break;
    }
    return dir;
}

class MeshSnake
{
public:
   MeshSnake(MeshBuilder* mesh);

   virtual ~MeshSnake() {}

   virtual void begin_snake();

   virtual void crawl(unsigned numSteps);


private:
   MeshBuilder* m_mesh;

   unsigned m_moveCounter;
   unsigned m_moveMax;

   void update_move_counter();

   //begin snake
   virtual void get_starting_direction() = 0;

   //crawl
   void move();
       bool time_to_switch_directions();

       void increment_move_counter();

       virtual bool facing_wall() = 0;

       virtual void get_new_direction() = 0;

       virtual void move_one() = 0;

};

class QuadMeshSnake : public MeshSnake
{
public:
    QuadMeshSnake(QuadMeshBuilder& mesh);

    virtual ~QuadMeshSnake() {}

    void set_x_bounds(unsigned xLower, unsigned xUpper);

    void set_y_bounds(unsigned yLower, unsigned yUpper);

    void set_x_pos(unsigned x);

    void set_y_pos(unsigned y);

    void set_pos(unsigned x, unsigned y);

    //test functions
    inline unsigned x_lower_bound() const;

    inline unsigned x_upper_bound() const;

    inline unsigned y_lower_bound() const;

    inline unsigned y_upper_bound() const;

    inline unsigned x_pos() const;

    inline unsigned y_pos() const;

    inline direction dir() const;

private:
    QuadMeshBuilder& m_quadMesh;

    unsigned m_xLower;
    unsigned m_xUpper;
    unsigned m_yLower;
    unsigned m_yUpper;

    unsigned m_xPos;
    unsigned m_yPos;

    direction m_dir;

    //begin snake
    virtual void get_starting_direction();

    //move
    bool wall_on_right();

    bool wall_on_left();

    bool wall_on_top();

    bool wall_on_bottom();

    virtual bool facing_wall();

    virtual void get_new_direction();

        void get_new_y_direction();

        void get_new_x_direction();

    virtual void move_one();

        void move_forward();


};

//test
inline unsigned QuadMeshSnake::x_lower_bound() const
{
    return m_xLower;
}
inline unsigned QuadMeshSnake::x_upper_bound() const
{
    return m_xUpper;
}
inline unsigned QuadMeshSnake::y_lower_bound() const
{
    return m_yLower;
}
inline unsigned QuadMeshSnake::y_upper_bound() const
{
    return m_yUpper;
}
inline unsigned QuadMeshSnake::x_pos() const
{
    return m_xPos;
}
inline unsigned QuadMeshSnake::y_pos() const
{
    return m_yPos;
}
inline direction QuadMeshSnake::dir() const
{
    return m_dir;
}

class HexMeshSnake : public MeshSnake
{
public:
    HexMeshSnake(HexMeshBuilder& mesh);

    virtual ~HexMeshSnake() {}

    void set_x_bounds(unsigned xLower, unsigned xUpper);

    void set_y_bounds(unsigned yLower, unsigned yUpper);

    void set_z_bounds(unsigned zLower, unsigned zUpper);

    void set_x_pos(unsigned x);

    void set_y_pos(unsigned y);

    void set_z_pos(unsigned z);

    void set_pos(unsigned x, unsigned y, unsigned z);

    //test
    inline unsigned x_lower_bound() const;

    inline unsigned x_upper_bound() const;

    inline unsigned y_lower_bound() const;

    inline unsigned y_upper_bound() const;

    inline unsigned z_lower_bound() const;

    inline unsigned z_upper_bound() const;

    inline unsigned x_pos() const;

    inline unsigned y_pos() const;

    inline unsigned z_pos() const;

    inline direction dir() const;

private:
    HexMeshBuilder& m_hexMesh;

    unsigned m_xLower;
    unsigned m_xUpper;
    unsigned m_yLower;
    unsigned m_yUpper;
    unsigned m_zLower;
    unsigned m_zUpper;

    unsigned m_xPos;
    unsigned m_yPos;
    unsigned m_zPos;

    direction m_dir;

    //begin snake
    virtual void get_starting_direction();

    //move
    bool wall_on_right();

    bool wall_on_left();

    bool wall_on_top();

    bool wall_on_bottom();

    bool wall_on_front();

    bool wall_on_back();

    virtual bool facing_wall();

    virtual void get_new_direction();

        void get_new_xy_direction();

        void get_new_xz_direction();

        void get_new_yz_direction();

    virtual void move_one();

        void move_forward();
};

inline unsigned HexMeshSnake::x_lower_bound() const
{
    return m_xLower;
}
inline unsigned HexMeshSnake::x_upper_bound() const
{
    return m_xUpper;
}
inline unsigned HexMeshSnake::y_lower_bound() const
{
    return m_yLower;
}
inline unsigned HexMeshSnake::y_upper_bound() const
{
    return m_yUpper;
}
inline unsigned HexMeshSnake::z_lower_bound() const
{
    return m_zLower;
}
inline unsigned HexMeshSnake::z_upper_bound() const
{
    return m_zUpper;
}
inline unsigned HexMeshSnake::x_pos() const
{
    return m_xPos;
}
inline unsigned HexMeshSnake::y_pos() const
{
    return m_yPos;
}
inline unsigned HexMeshSnake::z_pos() const
{
    return m_zPos;
}
inline direction HexMeshSnake::dir() const
{
    return m_dir;
}

#endif /*_MESHSNAKE_HPP_*/
