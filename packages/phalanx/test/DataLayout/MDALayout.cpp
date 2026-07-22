// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_RCP.hpp"
#include "Teuchos_UnitTestHarness.hpp"

// From test/Utilities directory
#include "Traits.hpp"

#include "Phalanx_DataLayout.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_ExtentTraits.hpp"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
PHX_EXTENT(Spatial)
PHX_EXTENT(Quadrature)
PHX_EXTENT(Node)
PHX_EXTENT(Cell)
PHX_EXTENT(Ordinal1)
PHX_EXTENT(Ordinal2)
PHX_EXTENT(Ordinal3)
PHX_EXTENT(Ordinal4)
PHX_EXTENT(Ordinal5)
PHX_EXTENT(Ordinal6)
PHX_EXTENT(Ordinal7)
PHX_EXTENT(Ordinal8)

namespace PHX {
  template<> std::string print<Spatial>(){return "Spatial";}
  template<> std::string print<Quadrature>(){return "Quadrature";}
  template<> std::string print<Node>(){return "Node";}
  template<> std::string print<Cell>(){return "Cell";}
  template<> std::string print<Ordinal1>(){return "Ordinal1";}
  template<> std::string print<Ordinal2>(){return "Ordinal2";}
  template<> std::string print<Ordinal3>(){return "Ordinal3";}
  template<> std::string print<Ordinal4>(){return "Ordinal4";}
  template<> std::string print<Ordinal5>(){return "Ordinal5";}
  template<> std::string print<Ordinal6>(){return "Ordinal6";}
  template<> std::string print<Ordinal7>(){return "Ordinal7";}
  template<> std::string print<Ordinal8>(){return "Ordinal8";}
}

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TEUCHOS_UNIT_TEST(DataLayout, basic)
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
  // extent()
  {
    TEST_EQUALITY(rank1.extent(0),1);

    TEST_EQUALITY(rank2.extent(0),1);
    TEST_EQUALITY(rank2.extent(1),2);

    TEST_EQUALITY(rank3.extent(0),1);
    TEST_EQUALITY(rank3.extent(1),2);
    TEST_EQUALITY(rank3.extent(2),3);

    TEST_EQUALITY(rank4.extent(0),1);
    TEST_EQUALITY(rank4.extent(1),2);
    TEST_EQUALITY(rank4.extent(2),3);
    TEST_EQUALITY(rank4.extent(3),4);

    TEST_EQUALITY(rank5.extent(0),1);
    TEST_EQUALITY(rank5.extent(1),2);
    TEST_EQUALITY(rank5.extent(2),3);
    TEST_EQUALITY(rank5.extent(3),4);
    TEST_EQUALITY(rank5.extent(4),5);

    TEST_EQUALITY(rank6.extent(0),1);
    TEST_EQUALITY(rank6.extent(1),2);
    TEST_EQUALITY(rank6.extent(2),3);
    TEST_EQUALITY(rank6.extent(3),4);
    TEST_EQUALITY(rank6.extent(4),5);
    TEST_EQUALITY(rank6.extent(5),6);

    TEST_EQUALITY(rank7.extent(0),1);
    TEST_EQUALITY(rank7.extent(1),2);
    TEST_EQUALITY(rank7.extent(2),3);
    TEST_EQUALITY(rank7.extent(3),4);
    TEST_EQUALITY(rank7.extent(4),5);
    TEST_EQUALITY(rank7.extent(5),6);
    TEST_EQUALITY(rank7.extent(6),7);

    TEST_EQUALITY(rank8.extent(0),1);
    TEST_EQUALITY(rank8.extent(1),2);
    TEST_EQUALITY(rank8.extent(2),3);
    TEST_EQUALITY(rank8.extent(3),4);
    TEST_EQUALITY(rank8.extent(4),5);
    TEST_EQUALITY(rank8.extent(5),6);
    TEST_EQUALITY(rank8.extent(6),7);
    TEST_EQUALITY(rank8.extent(7),8);

    TEST_EQUALITY(prank1.extent(0),1);

    TEST_EQUALITY(prank2.extent(0),1);
    TEST_EQUALITY(prank2.extent(1),2);

    TEST_EQUALITY(prank3.extent(0),1);
    TEST_EQUALITY(prank3.extent(1),2);
    TEST_EQUALITY(prank3.extent(2),3);

    TEST_EQUALITY(prank4.extent(0),1);
    TEST_EQUALITY(prank4.extent(1),2);
    TEST_EQUALITY(prank4.extent(2),3);
    TEST_EQUALITY(prank4.extent(3),4);

    TEST_EQUALITY(prank5.extent(0),1);
    TEST_EQUALITY(prank5.extent(1),2);
    TEST_EQUALITY(prank5.extent(2),3);
    TEST_EQUALITY(prank5.extent(3),4);
    TEST_EQUALITY(prank5.extent(4),5);

    TEST_EQUALITY(prank6.extent(0),1);
    TEST_EQUALITY(prank6.extent(1),2);
    TEST_EQUALITY(prank6.extent(2),3);
    TEST_EQUALITY(prank6.extent(3),4);
    TEST_EQUALITY(prank6.extent(4),5);
    TEST_EQUALITY(prank6.extent(5),6);

    TEST_EQUALITY(prank7.extent(0),1);
    TEST_EQUALITY(prank7.extent(1),2);
    TEST_EQUALITY(prank7.extent(2),3);
    TEST_EQUALITY(prank7.extent(3),4);
    TEST_EQUALITY(prank7.extent(4),5);
    TEST_EQUALITY(prank7.extent(5),6);
    TEST_EQUALITY(prank7.extent(6),7);

    TEST_EQUALITY(prank8.extent(0),1);
    TEST_EQUALITY(prank8.extent(1),2);
    TEST_EQUALITY(prank8.extent(2),3);
    TEST_EQUALITY(prank8.extent(3),4);
    TEST_EQUALITY(prank8.extent(4),5);
    TEST_EQUALITY(prank8.extent(5),6);
    TEST_EQUALITY(prank8.extent(6),7);
    TEST_EQUALITY(prank8.extent(7),8);

  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // extent_int()
  {
    TEST_EQUALITY(rank1.extent_int(0),1);

    TEST_EQUALITY(rank2.extent_int(0),1);
    TEST_EQUALITY(rank2.extent_int(1),2);

    TEST_EQUALITY(rank3.extent_int(0),1);
    TEST_EQUALITY(rank3.extent_int(1),2);
    TEST_EQUALITY(rank3.extent_int(2),3);

    TEST_EQUALITY(rank4.extent_int(0),1);
    TEST_EQUALITY(rank4.extent_int(1),2);
    TEST_EQUALITY(rank4.extent_int(2),3);
    TEST_EQUALITY(rank4.extent_int(3),4);

    TEST_EQUALITY(rank5.extent_int(0),1);
    TEST_EQUALITY(rank5.extent_int(1),2);
    TEST_EQUALITY(rank5.extent_int(2),3);
    TEST_EQUALITY(rank5.extent_int(3),4);
    TEST_EQUALITY(rank5.extent_int(4),5);

    TEST_EQUALITY(rank6.extent_int(0),1);
    TEST_EQUALITY(rank6.extent_int(1),2);
    TEST_EQUALITY(rank6.extent_int(2),3);
    TEST_EQUALITY(rank6.extent_int(3),4);
    TEST_EQUALITY(rank6.extent_int(4),5);
    TEST_EQUALITY(rank6.extent_int(5),6);

    TEST_EQUALITY(rank7.extent_int(0),1);
    TEST_EQUALITY(rank7.extent_int(1),2);
    TEST_EQUALITY(rank7.extent_int(2),3);
    TEST_EQUALITY(rank7.extent_int(3),4);
    TEST_EQUALITY(rank7.extent_int(4),5);
    TEST_EQUALITY(rank7.extent_int(5),6);
    TEST_EQUALITY(rank7.extent_int(6),7);

    TEST_EQUALITY(rank8.extent_int(0),1);
    TEST_EQUALITY(rank8.extent_int(1),2);
    TEST_EQUALITY(rank8.extent_int(2),3);
    TEST_EQUALITY(rank8.extent_int(3),4);
    TEST_EQUALITY(rank8.extent_int(4),5);
    TEST_EQUALITY(rank8.extent_int(5),6);
    TEST_EQUALITY(rank8.extent_int(6),7);
    TEST_EQUALITY(rank8.extent_int(7),8);

    TEST_EQUALITY(prank1.extent_int(0),1);

    TEST_EQUALITY(prank2.extent_int(0),1);
    TEST_EQUALITY(prank2.extent_int(1),2);

    TEST_EQUALITY(prank3.extent_int(0),1);
    TEST_EQUALITY(prank3.extent_int(1),2);
    TEST_EQUALITY(prank3.extent_int(2),3);

    TEST_EQUALITY(prank4.extent_int(0),1);
    TEST_EQUALITY(prank4.extent_int(1),2);
    TEST_EQUALITY(prank4.extent_int(2),3);
    TEST_EQUALITY(prank4.extent_int(3),4);

    TEST_EQUALITY(prank5.extent_int(0),1);
    TEST_EQUALITY(prank5.extent_int(1),2);
    TEST_EQUALITY(prank5.extent_int(2),3);
    TEST_EQUALITY(prank5.extent_int(3),4);
    TEST_EQUALITY(prank5.extent_int(4),5);

    TEST_EQUALITY(prank6.extent_int(0),1);
    TEST_EQUALITY(prank6.extent_int(1),2);
    TEST_EQUALITY(prank6.extent_int(2),3);
    TEST_EQUALITY(prank6.extent_int(3),4);
    TEST_EQUALITY(prank6.extent_int(4),5);
    TEST_EQUALITY(prank6.extent_int(5),6);

    TEST_EQUALITY(prank7.extent_int(0),1);
    TEST_EQUALITY(prank7.extent_int(1),2);
    TEST_EQUALITY(prank7.extent_int(2),3);
    TEST_EQUALITY(prank7.extent_int(3),4);
    TEST_EQUALITY(prank7.extent_int(4),5);
    TEST_EQUALITY(prank7.extent_int(5),6);
    TEST_EQUALITY(prank7.extent_int(6),7);

    TEST_EQUALITY(prank8.extent_int(0),1);
    TEST_EQUALITY(prank8.extent_int(1),2);
    TEST_EQUALITY(prank8.extent_int(2),3);
    TEST_EQUALITY(prank8.extent_int(3),4);
    TEST_EQUALITY(prank8.extent_int(4),5);
    TEST_EQUALITY(prank8.extent_int(5),6);
    TEST_EQUALITY(prank8.extent_int(6),7);
    TEST_EQUALITY(prank8.extent_int(7),8);

  }
  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // setExtents() from base class and derived class
  // NOTE: the derived class does not implement this method
  {
    DataLayout& b1 = rank1;
    b1.setExtents(10);
    TEST_EQUALITY(b1.extent(0),10);

    DataLayout& b2 = rank2;
    b2.setExtents(10,20);
    TEST_EQUALITY(b2.extent(0),10);
    TEST_EQUALITY(b2.extent(1),20);

    DataLayout& b3 = rank3;
    b3.setExtents(10,20,30);
    TEST_EQUALITY(b3.extent(0),10);
    TEST_EQUALITY(b3.extent(1),20);
    TEST_EQUALITY(b3.extent(2),30);

    DataLayout& b4 = rank4;
    b4.setExtents(10,20,30,40);
    TEST_EQUALITY(b4.extent(0),10);
    TEST_EQUALITY(b4.extent(1),20);
    TEST_EQUALITY(b4.extent(2),30);
    TEST_EQUALITY(b4.extent(3),40);

    DataLayout& b5 = rank5;
    b5.setExtents(10,20,30,40,50);
    TEST_EQUALITY(b5.extent(0),10);
    TEST_EQUALITY(b5.extent(1),20);
    TEST_EQUALITY(b5.extent(2),30);
    TEST_EQUALITY(b5.extent(3),40);
    TEST_EQUALITY(b5.extent(4),50);

    DataLayout& b6 = rank6;
    b6.setExtents(10,20,30,40,50,60);
    TEST_EQUALITY(b6.extent(0),10);
    TEST_EQUALITY(b6.extent(1),20);
    TEST_EQUALITY(b6.extent(2),30);
    TEST_EQUALITY(b6.extent(3),40);
    TEST_EQUALITY(b6.extent(4),50);
    TEST_EQUALITY(b6.extent(5),60);

    DataLayout& b7 = rank7;
    b7.setExtents(10,20,30,40,50,60,70);
    TEST_EQUALITY(b7.extent(0),10);
    TEST_EQUALITY(b7.extent(1),20);
    TEST_EQUALITY(b7.extent(2),30);
    TEST_EQUALITY(b7.extent(3),40);
    TEST_EQUALITY(b7.extent(4),50);
    TEST_EQUALITY(b7.extent(5),60);
    TEST_EQUALITY(b7.extent(6),70);

    DataLayout& b8 = rank8;
    // Set on the derived class instead of base
    rank8.setExtents(10,20,30,40,50,60,70,80);
    TEST_EQUALITY(b8.extent(0),10);
    TEST_EQUALITY(b8.extent(1),20);
    TEST_EQUALITY(b8.extent(2),30);
    TEST_EQUALITY(b8.extent(3),40);
    TEST_EQUALITY(b8.extent(4),50);
    TEST_EQUALITY(b8.extent(5),60);
    TEST_EQUALITY(b8.extent(6),70);
    TEST_EQUALITY(b8.extent(7),80);
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

  TEST_ASSERT(prank7.kokkosLayout() == PHX::DataLayout::KokkosLayoutType::Default);
}
