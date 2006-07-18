//@HEADER
// ***********************************************************************
//
//                     Rythmos Package
//                 Copyright (2006) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Todd S. Coffey (tscoffe@sandia.gov)
//
// ***********************************************************************
//@HEADER

class Bar
{
  public:
    Bar() {};
    ~Bar() {};
  protected:
    double t_;
};

class Foo
{
  public:
    Foo() { };
    ~Foo() { };
    const Bar &getBar() { return(B_); };
//    const Bar &getBar() const { return(B_); }; // this allows the line below to work.
  protected:
    Bar B_;
};

const Bar &somefunction(const Foo &F)
{
  const Bar &B = F.getBar(); // This doesn't work.  Why?
  return(B);
};

int main(int argc, char *argv[])
{
  Foo F;
  const Bar &B1 = F.getBar(); // This works.
  const Bar &B2 = somefunction(F);

  return(0);
}
