// @HEADER
// ************************************************************************
//
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//                    Copyright 2008 Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov), Sandia
// National Laboratories.
//
// ************************************************************************
// @HEADER


#include "Phalanx_config.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"

// From test/Utilities directory
#include "Traits.hpp"

#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_DimTag.hpp"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
struct Spatial : public PHX::DimTag {
  Spatial(){};
  const char * name() const ;
  static const Spatial & tag();
};

struct Quadrature : public PHX::DimTag {
  Quadrature(){};
  const char * name() const ;
  static const Quadrature & tag();
};

struct Node : public PHX::DimTag {
  Node(){};
  const char * name() const ;
  static const Node & tag();
};

struct Cell : public PHX::DimTag {
  Cell(){};
  const char * name() const ;
  static const Cell & tag();
};

struct Ordinal1 : public PHX::DimTag {
  Ordinal1(){};
  const char * name() const ;
  static const Ordinal1 & tag();
};

struct Ordinal2 : public PHX::DimTag {
  Ordinal2(){};
  const char * name() const ;
  static const Ordinal2 & tag();
};

struct Ordinal3 : public PHX::DimTag {
  Ordinal3(){};
  const char * name() const ;
  static const Ordinal3 & tag();
};

struct Ordinal4 : public PHX::DimTag {
  Ordinal4(){};
  const char * name() const ;
  static const Ordinal4 & tag();
};

struct Ordinal5 : public PHX::DimTag {
  Ordinal5(){};
  const char * name() const ;
  static const Ordinal5 & tag();
};

struct Ordinal6 : public PHX::DimTag {
  Ordinal6(){};
  const char * name() const ;
  static const Ordinal6 & tag();
};

struct Ordinal7 : public PHX::DimTag {
  Ordinal7(){};
  const char * name() const ;
  static const Ordinal7 & tag();
};

struct Ordinal8 : public PHX::DimTag {
  Ordinal8(){};
  const char * name() const ;
  static const Ordinal8 & tag();
};

const char * Spatial::name() const 
{ static const char n[] = "Spatial" ; return n ; }
const Spatial & Spatial::tag() 
{ static const Spatial myself ; return myself ; }

const char * Quadrature::name() const 
{ static const char n[] = "Quadrature" ; return n ; }
const Quadrature & Quadrature::tag() 
{ static const Quadrature myself ; return myself ; }

const char * Node::name() const 
{ static const char n[] = "Node" ; return n ; }
const Node & Node::tag() 
{ static const Node myself ; return myself ; }

const char * Cell::name() const 
{ static const char n[] = "Cell" ; return n ; }
const Cell & Cell::tag() 
{ static const Cell myself ; return myself ; }

const char * Ordinal1::name() const 
{ static const char n[] = "Ordinal1" ; return n ; }
const Ordinal1 & Ordinal1::tag() 
{ static const Ordinal1 myself ; return myself ; }

const char * Ordinal2::name() const 
{ static const char n[] = "Ordinal2" ; return n ; }
const Ordinal2 & Ordinal2::tag() 
{ static const Ordinal2 myself ; return myself ; }

const char * Ordinal3::name() const 
{ static const char n[] = "Ordinal3" ; return n ; }
const Ordinal3 & Ordinal3::tag() 
{ static const Ordinal3 myself ; return myself ; }

const char * Ordinal4::name() const 
{ static const char n[] = "Ordinal4" ; return n ; }
const Ordinal4 & Ordinal4::tag() 
{ static const Ordinal4 myself ; return myself ; }

const char * Ordinal5::name() const 
{ static const char n[] = "Ordinal5" ; return n ; }
const Ordinal5 & Ordinal5::tag() 
{ static const Ordinal5 myself ; return myself ; }

const char * Ordinal6::name() const 
{ static const char n[] = "Ordinal6" ; return n ; }
const Ordinal6 & Ordinal6::tag() 
{ static const Ordinal6 myself ; return myself ; }

const char * Ordinal7::name() const 
{ static const char n[] = "Ordinal7" ; return n ; }
const Ordinal7 & Ordinal7::tag() 
{ static const Ordinal7 myself ; return myself ; }

const char * Ordinal8::name() const 
{ static const char n[] = "Ordinal8" ; return n ; }
const Ordinal8 & Ordinal8::tag() 
{ static const Ordinal8 myself ; return myself ; }

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TEUCHOS_UNIT_TEST(DataLayout, DataLayout)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ctor
  cout << "\nTesting constructor...";    
  MDALayout<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5,Ordinal6,Ordinal7,Ordinal8> rank8(1,2,3,4,5,6,7,8);
  MDALayout<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5,Ordinal6,Ordinal7> rank7(1,2,3,4,5,6,7);
  MDALayout<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5,Ordinal6> rank6(1,2,3,4,5,6);
  MDALayout<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5> rank5(1,2,3,4,5);
  MDALayout<Ordinal1,Ordinal2,Ordinal3,Ordinal4> rank4(1,2,3,4);
  MDALayout<Ordinal1,Ordinal2,Ordinal3> rank3(1,2,3);
  MDALayout<Ordinal1,Ordinal2> rank2(1,2);
  MDALayout<Ordinal1> rank1(1);
  MDALayout<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5,Ordinal6,Ordinal7,Ordinal8> prank8("DL8",1,2,3,4,5,6,7,8);
  MDALayout<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5,Ordinal6,Ordinal7> prank7("DL7",1,2,3,4,5,6,7);
  MDALayout<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5,Ordinal6> prank6("DL6",1,2,3,4,5,6);
  MDALayout<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5> prank5("DL5",1,2,3,4,5);
  MDALayout<Ordinal1,Ordinal2,Ordinal3,Ordinal4> prank4("DL4",1,2,3,4);
  MDALayout<Ordinal1,Ordinal2,Ordinal3> prank3("DL3",1,2,3);
  MDALayout<Ordinal1,Ordinal2> prank2("DL2",1,2);
  MDALayout<Ordinal1> prank1("DL1",1);
  MDALayout<Cell,Node,Spatial,Spatial> n_mat(100,4,2,2);
  cout << "passed!" << endl;
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // identifier()
  TEST_EQUALITY(rank8.identifier(),"<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5,Ordinal6,Ordinal7,Ordinal8>(1,2,3,4,5,6,7,8)");
  TEST_EQUALITY(rank7.identifier(),"<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5,Ordinal6,Ordinal7>(1,2,3,4,5,6,7)");
  TEST_EQUALITY(rank6.identifier(),"<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5,Ordinal6>(1,2,3,4,5,6)");
  TEST_EQUALITY(rank5.identifier(),"<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5>(1,2,3,4,5)");
  TEST_EQUALITY(rank4.identifier(),"<Ordinal1,Ordinal2,Ordinal3,Ordinal4>(1,2,3,4)");
  TEST_EQUALITY(rank3.identifier(),"<Ordinal1,Ordinal2,Ordinal3>(1,2,3)");
  TEST_EQUALITY(rank2.identifier(),"<Ordinal1,Ordinal2>(1,2)");
  TEST_EQUALITY(rank1.identifier(),"<Ordinal1>(1)");
  TEST_EQUALITY(prank8.identifier(),"DL8<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5,Ordinal6,Ordinal7,Ordinal8>(1,2,3,4,5,6,7,8)");
  TEST_EQUALITY(prank7.identifier(),"DL7<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5,Ordinal6,Ordinal7>(1,2,3,4,5,6,7)");
  TEST_EQUALITY(prank6.identifier(),"DL6<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5,Ordinal6>(1,2,3,4,5,6)");
  TEST_EQUALITY(prank5.identifier(),"DL5<Ordinal1,Ordinal2,Ordinal3,Ordinal4,Ordinal5>(1,2,3,4,5)");
  TEST_EQUALITY(prank4.identifier(),"DL4<Ordinal1,Ordinal2,Ordinal3,Ordinal4>(1,2,3,4)");
  TEST_EQUALITY(prank3.identifier(),"DL3<Ordinal1,Ordinal2,Ordinal3>(1,2,3)");
  TEST_EQUALITY(prank2.identifier(),"DL2<Ordinal1,Ordinal2>(1,2)");
  TEST_EQUALITY(prank1.identifier(),"DL1<Ordinal1>(1)");
  TEST_EQUALITY(n_mat.identifier(),"<Cell,Node,Spatial,Spatial>(100,4,2,2)");
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // rank()
  TEST_EQUALITY(rank8.rank(),8);
  TEST_EQUALITY(rank7.rank(),7);
  TEST_EQUALITY(rank6.rank(),6);
  TEST_EQUALITY(rank5.rank(),5);
  TEST_EQUALITY(rank4.rank(),4);
  TEST_EQUALITY(rank3.rank(),3);
  TEST_EQUALITY(rank2.rank(),2);
  TEST_EQUALITY(rank1.rank(),1);
  TEST_EQUALITY(prank8.rank(),8);
  TEST_EQUALITY(prank7.rank(),7);
  TEST_EQUALITY(prank6.rank(),6);
  TEST_EQUALITY(prank5.rank(),5);
  TEST_EQUALITY(prank4.rank(),4);
  TEST_EQUALITY(prank3.rank(),3);
  TEST_EQUALITY(prank2.rank(),2);
  TEST_EQUALITY(prank1.rank(),1);
  TEST_EQUALITY(n_mat.rank(),4);
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // dimensions()
  {
    std::vector<PHX::Device::size_type> dims;
    n_mat.dimensions(dims);
    TEST_EQUALITY(dims[0],100);
    TEST_EQUALITY(dims[1],4);
    TEST_EQUALITY(dims[2],2);
    TEST_EQUALITY(dims[3],2);
  }
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // names()
  {
    std::vector<std::string> names;
    n_mat.names(names);
    TEST_EQUALITY(names[0],"Cell");
    TEST_EQUALITY(names[1],"Node");
    TEST_EQUALITY(names[2],"Spatial");
    TEST_EQUALITY(names[3],"Spatial");
  }
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // size()
  TEST_EQUALITY(n_mat.size(),1600);
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // operator==()
  {
    MDALayout<Cell,Node,Spatial> n_vec_a(100,4,2);
    MDALayout<Cell,Node,Spatial> n_vec_b(50,4,2);
    MDALayout<Cell,Node,Spatial> n_vec_c(50,4,2);
    MDALayout<Cell,Quadrature,Spatial> qp_vec(100,4,2);
    MDALayout<Cell,Quadrature,Spatial,Spatial> qp_mat(100,4,2,2);
    
    // Same data layout, different objects
    TEST_EQUALITY(n_vec_b,n_vec_c);
    
    // Same ordinal types, different dim sizes
    TEST_INEQUALITY(n_vec_a,n_vec_c);
    
    // Different ordinal types, same rank, same dim sizes
    TEST_INEQUALITY(n_vec_a,qp_vec);
    
    // Different ordianl types, different dim sizes
    TEST_INEQUALITY(n_vec_a,qp_mat);
  }
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // dimension()
  {
    TEST_EQUALITY(rank1.dimension(0),1);
    
    TEST_EQUALITY(rank2.dimension(0),1);
    TEST_EQUALITY(rank2.dimension(1),2);
    
    TEST_EQUALITY(rank3.dimension(0),1);
    TEST_EQUALITY(rank3.dimension(1),2);
    TEST_EQUALITY(rank3.dimension(2),3);
    
    TEST_EQUALITY(rank4.dimension(0),1);
    TEST_EQUALITY(rank4.dimension(1),2);
    TEST_EQUALITY(rank4.dimension(2),3);
    TEST_EQUALITY(rank4.dimension(3),4);
    
    TEST_EQUALITY(rank5.dimension(0),1);
    TEST_EQUALITY(rank5.dimension(1),2);
    TEST_EQUALITY(rank5.dimension(2),3);
    TEST_EQUALITY(rank5.dimension(3),4);
    TEST_EQUALITY(rank5.dimension(4),5);
    
    TEST_EQUALITY(rank6.dimension(0),1);
    TEST_EQUALITY(rank6.dimension(1),2);
    TEST_EQUALITY(rank6.dimension(2),3);
    TEST_EQUALITY(rank6.dimension(3),4);
    TEST_EQUALITY(rank6.dimension(4),5);
    TEST_EQUALITY(rank6.dimension(5),6);
    
    TEST_EQUALITY(rank7.dimension(0),1);
    TEST_EQUALITY(rank7.dimension(1),2);
    TEST_EQUALITY(rank7.dimension(2),3);
    TEST_EQUALITY(rank7.dimension(3),4);
    TEST_EQUALITY(rank7.dimension(4),5);
    TEST_EQUALITY(rank7.dimension(5),6);
    TEST_EQUALITY(rank7.dimension(6),7);
    
    TEST_EQUALITY(rank8.dimension(0),1);
    TEST_EQUALITY(rank8.dimension(1),2);
    TEST_EQUALITY(rank8.dimension(2),3);
    TEST_EQUALITY(rank8.dimension(3),4);
    TEST_EQUALITY(rank8.dimension(4),5);
    TEST_EQUALITY(rank8.dimension(5),6);
    TEST_EQUALITY(rank8.dimension(6),7);
    TEST_EQUALITY(rank8.dimension(7),8);

    TEST_EQUALITY(prank1.dimension(0),1);
    
    TEST_EQUALITY(prank2.dimension(0),1);
    TEST_EQUALITY(prank2.dimension(1),2);
    
    TEST_EQUALITY(prank3.dimension(0),1);
    TEST_EQUALITY(prank3.dimension(1),2);
    TEST_EQUALITY(prank3.dimension(2),3);
    
    TEST_EQUALITY(prank4.dimension(0),1);
    TEST_EQUALITY(prank4.dimension(1),2);
    TEST_EQUALITY(prank4.dimension(2),3);
    TEST_EQUALITY(prank4.dimension(3),4);
    
    TEST_EQUALITY(prank5.dimension(0),1);
    TEST_EQUALITY(prank5.dimension(1),2);
    TEST_EQUALITY(prank5.dimension(2),3);
    TEST_EQUALITY(prank5.dimension(3),4);
    TEST_EQUALITY(prank5.dimension(4),5);
    
    TEST_EQUALITY(prank6.dimension(0),1);
    TEST_EQUALITY(prank6.dimension(1),2);
    TEST_EQUALITY(prank6.dimension(2),3);
    TEST_EQUALITY(prank6.dimension(3),4);
    TEST_EQUALITY(prank6.dimension(4),5);
    TEST_EQUALITY(prank6.dimension(5),6);
    
    TEST_EQUALITY(prank7.dimension(0),1);
    TEST_EQUALITY(prank7.dimension(1),2);
    TEST_EQUALITY(prank7.dimension(2),3);
    TEST_EQUALITY(prank7.dimension(3),4);
    TEST_EQUALITY(prank7.dimension(4),5);
    TEST_EQUALITY(prank7.dimension(5),6);
    TEST_EQUALITY(prank7.dimension(6),7);
    
    TEST_EQUALITY(prank8.dimension(0),1);
    TEST_EQUALITY(prank8.dimension(1),2);
    TEST_EQUALITY(prank8.dimension(2),3);
    TEST_EQUALITY(prank8.dimension(3),4);
    TEST_EQUALITY(prank8.dimension(4),5);
    TEST_EQUALITY(prank8.dimension(5),6);
    TEST_EQUALITY(prank8.dimension(6),7);
    TEST_EQUALITY(prank8.dimension(7),8);

  }
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // name() 
  {
    cout << "\n" << rank1.name(0) << std::endl;
    
    TEST_EQUALITY(rank1.name(0),"Ordinal1");
    
    TEST_EQUALITY(rank2.name(0),"Ordinal1");
    TEST_EQUALITY(rank2.name(1),"Ordinal2");
    
    TEST_EQUALITY(rank3.name(0),"Ordinal1");
    TEST_EQUALITY(rank3.name(1),"Ordinal2");
    TEST_EQUALITY(rank3.name(2),"Ordinal3");
    
    TEST_EQUALITY(rank4.name(0),"Ordinal1");
    TEST_EQUALITY(rank4.name(1),"Ordinal2");
    TEST_EQUALITY(rank4.name(2),"Ordinal3");
    TEST_EQUALITY(rank4.name(3),"Ordinal4");
    
    TEST_EQUALITY(rank5.name(0),"Ordinal1");
    TEST_EQUALITY(rank5.name(1),"Ordinal2");
    TEST_EQUALITY(rank5.name(2),"Ordinal3");
    TEST_EQUALITY(rank5.name(3),"Ordinal4");
    TEST_EQUALITY(rank5.name(4),"Ordinal5");

    TEST_EQUALITY(rank6.name(0),"Ordinal1");
    TEST_EQUALITY(rank6.name(1),"Ordinal2");
    TEST_EQUALITY(rank6.name(2),"Ordinal3");
    TEST_EQUALITY(rank6.name(3),"Ordinal4");
    TEST_EQUALITY(rank6.name(4),"Ordinal5");
    TEST_EQUALITY(rank6.name(5),"Ordinal6");
    
    TEST_EQUALITY(rank7.name(0),"Ordinal1");
    TEST_EQUALITY(rank7.name(1),"Ordinal2");
    TEST_EQUALITY(rank7.name(2),"Ordinal3");
    TEST_EQUALITY(rank7.name(3),"Ordinal4");
    TEST_EQUALITY(rank7.name(4),"Ordinal5");
    TEST_EQUALITY(rank7.name(5),"Ordinal6");
    TEST_EQUALITY(rank7.name(6),"Ordinal7");
    
    TEST_EQUALITY(rank8.name(0),"Ordinal1");
    TEST_EQUALITY(rank8.name(1),"Ordinal2");
    TEST_EQUALITY(rank8.name(2),"Ordinal3");
    TEST_EQUALITY(rank8.name(3),"Ordinal4");
    TEST_EQUALITY(rank8.name(4),"Ordinal5");
    TEST_EQUALITY(rank8.name(5),"Ordinal6");
    TEST_EQUALITY(rank8.name(6),"Ordinal7");
    TEST_EQUALITY(rank8.name(7),"Ordinal8");
    
    TEST_EQUALITY(prank1.name(0),"Ordinal1");
    
    TEST_EQUALITY(prank2.name(0),"Ordinal1");
    TEST_EQUALITY(prank2.name(1),"Ordinal2");
    
    TEST_EQUALITY(prank3.name(0),"Ordinal1");
    TEST_EQUALITY(prank3.name(1),"Ordinal2");
    TEST_EQUALITY(prank3.name(2),"Ordinal3");
    
    TEST_EQUALITY(prank4.name(0),"Ordinal1");
    TEST_EQUALITY(prank4.name(1),"Ordinal2");
    TEST_EQUALITY(prank4.name(2),"Ordinal3");
    TEST_EQUALITY(prank4.name(3),"Ordinal4");
    
    TEST_EQUALITY(prank5.name(0),"Ordinal1");
    TEST_EQUALITY(prank5.name(1),"Ordinal2");
    TEST_EQUALITY(prank5.name(2),"Ordinal3");
    TEST_EQUALITY(prank5.name(3),"Ordinal4");
    TEST_EQUALITY(prank5.name(4),"Ordinal5");

    TEST_EQUALITY(prank6.name(0),"Ordinal1");
    TEST_EQUALITY(prank6.name(1),"Ordinal2");
    TEST_EQUALITY(prank6.name(2),"Ordinal3");
    TEST_EQUALITY(prank6.name(3),"Ordinal4");
    TEST_EQUALITY(prank6.name(4),"Ordinal5");
    TEST_EQUALITY(prank6.name(5),"Ordinal6");
    
    TEST_EQUALITY(prank7.name(0),"Ordinal1");
    TEST_EQUALITY(prank7.name(1),"Ordinal2");
    TEST_EQUALITY(prank7.name(2),"Ordinal3");
    TEST_EQUALITY(prank7.name(3),"Ordinal4");
    TEST_EQUALITY(prank7.name(4),"Ordinal5");
    TEST_EQUALITY(prank7.name(5),"Ordinal6");
    TEST_EQUALITY(prank7.name(6),"Ordinal7");
    
    TEST_EQUALITY(prank8.name(0),"Ordinal1");
    TEST_EQUALITY(prank8.name(1),"Ordinal2");
    TEST_EQUALITY(prank8.name(2),"Ordinal3");
    TEST_EQUALITY(prank8.name(3),"Ordinal4");
    TEST_EQUALITY(prank8.name(4),"Ordinal5");
    TEST_EQUALITY(prank8.name(5),"Ordinal6");
    TEST_EQUALITY(prank8.name(6),"Ordinal7");
    TEST_EQUALITY(prank8.name(7),"Ordinal8");
  }
  
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ostream
  {
    cout << "Testing ostream...";
    ostringstream output;
    output << rank8 << endl;
    output << rank7 << endl;
    output << rank6 << endl;   
    output << rank5 << endl;
    output << rank4 << endl;
    output << rank3 << endl;
    output << rank2 << endl;
    output << rank1 << endl;
    output << prank8 << endl;
    output << prank7 << endl;
    output << prank6 << endl;   
    output << prank5 << endl;
    output << prank4 << endl;
    output << prank3 << endl;
    output << prank2 << endl;
    output << prank1 << endl;
    output << n_mat << endl;
    cout << "...passed:\n" << output.str() << endl;
  }
  
}
