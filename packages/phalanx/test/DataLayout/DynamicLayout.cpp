// @HEADER
// *****************************************************************************
//        Phalanx: A Partial Differential Equation Field Evaluation 
//       Kernel for Flexible Management of Complex Dependency Chains
//
// Copyright 2008 NTESS and the Phalanx contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Teuchos_UnitTestHarness.hpp"

#include "Traits.hpp" // From test/Utilities directory
#include "Phalanx_DataLayout_DynamicLayout.hpp"

// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
TEUCHOS_UNIT_TEST(DynamicLayout, basic)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ctor with identifier and sizes
  Layout dl8("dl8",1,2,3,4,5,6,7,8);
  Layout dl7("dl7",1,2,3,4,5,6,7);
  Layout dl6("dl6",1,2,3,4,5,6);
  Layout dl5("dl5",1,2,3,4,5);
  Layout dl4("dl4",1,2,3,4);
  Layout dl3("dl3",1,2,3);
  Layout dl2("dl2",1,2);
  Layout dl1("dl1",1);

  // ctor with identifier and delayed sizing
  Layout e8("e8");
  Layout e7("e7");
  Layout e6("e6");
  Layout e5("e5");
  Layout e4("e4");
  Layout e3("e3");
  Layout e2("e2");
  Layout e1("e1");
  e8.setExtents(1,2,3,4,5,6,7,8);
  e7.setExtents(1,2,3,4,5,6,7);
  e6.setExtents(1,2,3,4,5,6);
  e5.setExtents(1,2,3,4,5);
  e4.setExtents(1,2,3,4);
  e3.setExtents(1,2,3);
  e2.setExtents(1,2);
  e1.setExtents(1);

  // ctor with sizing only not implemented yet.
  //Layout ni8(1,2,3,4,5,6,7,8);

  Layout n_mat("n_mat",100,4,2,2);

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // identifier()
  TEST_EQUALITY(dl1.identifier(),"dl1");
  TEST_EQUALITY(dl2.identifier(),"dl2");
  TEST_EQUALITY(dl3.identifier(),"dl3");
  TEST_EQUALITY(dl4.identifier(),"dl4");
  TEST_EQUALITY(dl5.identifier(),"dl5");
  TEST_EQUALITY(dl6.identifier(),"dl6");
  TEST_EQUALITY(dl7.identifier(),"dl7");
  TEST_EQUALITY(dl8.identifier(),"dl8");
  TEST_EQUALITY(e1.identifier(),"e1");
  TEST_EQUALITY(e2.identifier(),"e2");
  TEST_EQUALITY(e3.identifier(),"e3");
  TEST_EQUALITY(e4.identifier(),"e4");
  TEST_EQUALITY(e5.identifier(),"e5");
  TEST_EQUALITY(e6.identifier(),"e6");
  TEST_EQUALITY(e7.identifier(),"e7");
  TEST_EQUALITY(e8.identifier(),"e8");

  TEST_EQUALITY(n_mat.identifier(),"n_mat");

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // rank()
  TEST_EQUALITY(dl8.rank(),8);
  TEST_EQUALITY(dl7.rank(),7);
  TEST_EQUALITY(dl6.rank(),6);
  TEST_EQUALITY(dl5.rank(),5);
  TEST_EQUALITY(dl4.rank(),4);
  TEST_EQUALITY(dl3.rank(),3);
  TEST_EQUALITY(dl2.rank(),2);
  TEST_EQUALITY(dl1.rank(),1);
  TEST_EQUALITY(e8.rank(),8);
  TEST_EQUALITY(e7.rank(),7);
  TEST_EQUALITY(e6.rank(),6);
  TEST_EQUALITY(e5.rank(),5);
  TEST_EQUALITY(e4.rank(),4);
  TEST_EQUALITY(e3.rank(),3);
  TEST_EQUALITY(e2.rank(),2);
  TEST_EQUALITY(e1.rank(),1);
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
    TEST_EQUALITY(names[0],"EXT0");
    TEST_EQUALITY(names[1],"EXT1");
    TEST_EQUALITY(names[2],"EXT2");
    TEST_EQUALITY(names[3],"EXT3");
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // size()
  TEST_EQUALITY(n_mat.size(),1600);

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // operator==()
  {
    Layout a("a",50,4,2); // base
    Layout b("a",50,4,2); // matches a, but different object
    Layout c("c",50,4,2); // different identifier than a
    Layout d("a",50,4);   // different rank than a
    Layout e("a",50,4,3); // different extent than a

    // Returns true if identifier rank and extents match

    // Same data layout, different objects
    TEST_EQUALITY(a,b);

    // Different identifier
    TEST_INEQUALITY(a,c);

    // Different rank
    TEST_INEQUALITY(a,d);

    // Different extent
    TEST_INEQUALITY(a,e);
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // dimension()
  {
    TEST_EQUALITY(dl1.dimension(0),1);

    TEST_EQUALITY(dl2.dimension(0),1);
    TEST_EQUALITY(dl2.dimension(1),2);

    TEST_EQUALITY(dl3.dimension(0),1);
    TEST_EQUALITY(dl3.dimension(1),2);
    TEST_EQUALITY(dl3.dimension(2),3);

    TEST_EQUALITY(dl4.dimension(0),1);
    TEST_EQUALITY(dl4.dimension(1),2);
    TEST_EQUALITY(dl4.dimension(2),3);
    TEST_EQUALITY(dl4.dimension(3),4);

    TEST_EQUALITY(dl5.dimension(0),1);
    TEST_EQUALITY(dl5.dimension(1),2);
    TEST_EQUALITY(dl5.dimension(2),3);
    TEST_EQUALITY(dl5.dimension(3),4);
    TEST_EQUALITY(dl5.dimension(4),5);

    TEST_EQUALITY(dl6.dimension(0),1);
    TEST_EQUALITY(dl6.dimension(1),2);
    TEST_EQUALITY(dl6.dimension(2),3);
    TEST_EQUALITY(dl6.dimension(3),4);
    TEST_EQUALITY(dl6.dimension(4),5);
    TEST_EQUALITY(dl6.dimension(5),6);

    TEST_EQUALITY(dl7.dimension(0),1);
    TEST_EQUALITY(dl7.dimension(1),2);
    TEST_EQUALITY(dl7.dimension(2),3);
    TEST_EQUALITY(dl7.dimension(3),4);
    TEST_EQUALITY(dl7.dimension(4),5);
    TEST_EQUALITY(dl7.dimension(5),6);
    TEST_EQUALITY(dl7.dimension(6),7);

    TEST_EQUALITY(dl8.dimension(0),1);
    TEST_EQUALITY(dl8.dimension(1),2);
    TEST_EQUALITY(dl8.dimension(2),3);
    TEST_EQUALITY(dl8.dimension(3),4);
    TEST_EQUALITY(dl8.dimension(4),5);
    TEST_EQUALITY(dl8.dimension(5),6);
    TEST_EQUALITY(dl8.dimension(6),7);
    TEST_EQUALITY(dl8.dimension(7),8);

    TEST_EQUALITY(e1.dimension(0),1);

    TEST_EQUALITY(e2.dimension(0),1);
    TEST_EQUALITY(e2.dimension(1),2);

    TEST_EQUALITY(e3.dimension(0),1);
    TEST_EQUALITY(e3.dimension(1),2);
    TEST_EQUALITY(e3.dimension(2),3);

    TEST_EQUALITY(e4.dimension(0),1);
    TEST_EQUALITY(e4.dimension(1),2);
    TEST_EQUALITY(e4.dimension(2),3);
    TEST_EQUALITY(e4.dimension(3),4);

    TEST_EQUALITY(e5.dimension(0),1);
    TEST_EQUALITY(e5.dimension(1),2);
    TEST_EQUALITY(e5.dimension(2),3);
    TEST_EQUALITY(e5.dimension(3),4);
    TEST_EQUALITY(e5.dimension(4),5);

    TEST_EQUALITY(e6.dimension(0),1);
    TEST_EQUALITY(e6.dimension(1),2);
    TEST_EQUALITY(e6.dimension(2),3);
    TEST_EQUALITY(e6.dimension(3),4);
    TEST_EQUALITY(e6.dimension(4),5);
    TEST_EQUALITY(e6.dimension(5),6);

    TEST_EQUALITY(e7.dimension(0),1);
    TEST_EQUALITY(e7.dimension(1),2);
    TEST_EQUALITY(e7.dimension(2),3);
    TEST_EQUALITY(e7.dimension(3),4);
    TEST_EQUALITY(e7.dimension(4),5);
    TEST_EQUALITY(e7.dimension(5),6);
    TEST_EQUALITY(e7.dimension(6),7);

    TEST_EQUALITY(e8.dimension(0),1);
    TEST_EQUALITY(e8.dimension(1),2);
    TEST_EQUALITY(e8.dimension(2),3);
    TEST_EQUALITY(e8.dimension(3),4);
    TEST_EQUALITY(e8.dimension(4),5);
    TEST_EQUALITY(e8.dimension(5),6);
    TEST_EQUALITY(e8.dimension(6),7);
    TEST_EQUALITY(e8.dimension(7),8);
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // extent()
  {
    TEST_EQUALITY(dl1.extent(0),1);

    TEST_EQUALITY(dl2.extent(0),1);
    TEST_EQUALITY(dl2.extent(1),2);

    TEST_EQUALITY(dl3.extent(0),1);
    TEST_EQUALITY(dl3.extent(1),2);
    TEST_EQUALITY(dl3.extent(2),3);

    TEST_EQUALITY(dl4.extent(0),1);
    TEST_EQUALITY(dl4.extent(1),2);
    TEST_EQUALITY(dl4.extent(2),3);
    TEST_EQUALITY(dl4.extent(3),4);

    TEST_EQUALITY(dl5.extent(0),1);
    TEST_EQUALITY(dl5.extent(1),2);
    TEST_EQUALITY(dl5.extent(2),3);
    TEST_EQUALITY(dl5.extent(3),4);
    TEST_EQUALITY(dl5.extent(4),5);

    TEST_EQUALITY(dl6.extent(0),1);
    TEST_EQUALITY(dl6.extent(1),2);
    TEST_EQUALITY(dl6.extent(2),3);
    TEST_EQUALITY(dl6.extent(3),4);
    TEST_EQUALITY(dl6.extent(4),5);
    TEST_EQUALITY(dl6.extent(5),6);

    TEST_EQUALITY(dl7.extent(0),1);
    TEST_EQUALITY(dl7.extent(1),2);
    TEST_EQUALITY(dl7.extent(2),3);
    TEST_EQUALITY(dl7.extent(3),4);
    TEST_EQUALITY(dl7.extent(4),5);
    TEST_EQUALITY(dl7.extent(5),6);
    TEST_EQUALITY(dl7.extent(6),7);

    TEST_EQUALITY(dl8.extent(0),1);
    TEST_EQUALITY(dl8.extent(1),2);
    TEST_EQUALITY(dl8.extent(2),3);
    TEST_EQUALITY(dl8.extent(3),4);
    TEST_EQUALITY(dl8.extent(4),5);
    TEST_EQUALITY(dl8.extent(5),6);
    TEST_EQUALITY(dl8.extent(6),7);
    TEST_EQUALITY(dl8.extent(7),8);

    TEST_EQUALITY(e1.extent(0),1);

    TEST_EQUALITY(e2.extent(0),1);
    TEST_EQUALITY(e2.extent(1),2);

    TEST_EQUALITY(e3.extent(0),1);
    TEST_EQUALITY(e3.extent(1),2);
    TEST_EQUALITY(e3.extent(2),3);

    TEST_EQUALITY(e4.extent(0),1);
    TEST_EQUALITY(e4.extent(1),2);
    TEST_EQUALITY(e4.extent(2),3);
    TEST_EQUALITY(e4.extent(3),4);

    TEST_EQUALITY(e5.extent(0),1);
    TEST_EQUALITY(e5.extent(1),2);
    TEST_EQUALITY(e5.extent(2),3);
    TEST_EQUALITY(e5.extent(3),4);
    TEST_EQUALITY(e5.extent(4),5);

    TEST_EQUALITY(e6.extent(0),1);
    TEST_EQUALITY(e6.extent(1),2);
    TEST_EQUALITY(e6.extent(2),3);
    TEST_EQUALITY(e6.extent(3),4);
    TEST_EQUALITY(e6.extent(4),5);
    TEST_EQUALITY(e6.extent(5),6);

    TEST_EQUALITY(e7.extent(0),1);
    TEST_EQUALITY(e7.extent(1),2);
    TEST_EQUALITY(e7.extent(2),3);
    TEST_EQUALITY(e7.extent(3),4);
    TEST_EQUALITY(e7.extent(4),5);
    TEST_EQUALITY(e7.extent(5),6);
    TEST_EQUALITY(e7.extent(6),7);

    TEST_EQUALITY(e8.extent(0),1);
    TEST_EQUALITY(e8.extent(1),2);
    TEST_EQUALITY(e8.extent(2),3);
    TEST_EQUALITY(e8.extent(3),4);
    TEST_EQUALITY(e8.extent(4),5);
    TEST_EQUALITY(e8.extent(5),6);
    TEST_EQUALITY(e8.extent(6),7);
    TEST_EQUALITY(e8.extent(7),8);
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // extent_int()
  {
    TEST_EQUALITY(dl1.extent_int(0),1);

    TEST_EQUALITY(dl2.extent_int(0),1);
    TEST_EQUALITY(dl2.extent_int(1),2);

    TEST_EQUALITY(dl3.extent_int(0),1);
    TEST_EQUALITY(dl3.extent_int(1),2);
    TEST_EQUALITY(dl3.extent_int(2),3);

    TEST_EQUALITY(dl4.extent_int(0),1);
    TEST_EQUALITY(dl4.extent_int(1),2);
    TEST_EQUALITY(dl4.extent_int(2),3);
    TEST_EQUALITY(dl4.extent_int(3),4);

    TEST_EQUALITY(dl5.extent_int(0),1);
    TEST_EQUALITY(dl5.extent_int(1),2);
    TEST_EQUALITY(dl5.extent_int(2),3);
    TEST_EQUALITY(dl5.extent_int(3),4);
    TEST_EQUALITY(dl5.extent_int(4),5);

    TEST_EQUALITY(dl6.extent_int(0),1);
    TEST_EQUALITY(dl6.extent_int(1),2);
    TEST_EQUALITY(dl6.extent_int(2),3);
    TEST_EQUALITY(dl6.extent_int(3),4);
    TEST_EQUALITY(dl6.extent_int(4),5);
    TEST_EQUALITY(dl6.extent_int(5),6);

    TEST_EQUALITY(dl7.extent_int(0),1);
    TEST_EQUALITY(dl7.extent_int(1),2);
    TEST_EQUALITY(dl7.extent_int(2),3);
    TEST_EQUALITY(dl7.extent_int(3),4);
    TEST_EQUALITY(dl7.extent_int(4),5);
    TEST_EQUALITY(dl7.extent_int(5),6);
    TEST_EQUALITY(dl7.extent_int(6),7);

    TEST_EQUALITY(dl8.extent_int(0),1);
    TEST_EQUALITY(dl8.extent_int(1),2);
    TEST_EQUALITY(dl8.extent_int(2),3);
    TEST_EQUALITY(dl8.extent_int(3),4);
    TEST_EQUALITY(dl8.extent_int(4),5);
    TEST_EQUALITY(dl8.extent_int(5),6);
    TEST_EQUALITY(dl8.extent_int(6),7);
    TEST_EQUALITY(dl8.extent_int(7),8);

    TEST_EQUALITY(e1.extent_int(0),1);

    TEST_EQUALITY(e2.extent_int(0),1);
    TEST_EQUALITY(e2.extent_int(1),2);

    TEST_EQUALITY(e3.extent_int(0),1);
    TEST_EQUALITY(e3.extent_int(1),2);
    TEST_EQUALITY(e3.extent_int(2),3);

    TEST_EQUALITY(e4.extent_int(0),1);
    TEST_EQUALITY(e4.extent_int(1),2);
    TEST_EQUALITY(e4.extent_int(2),3);
    TEST_EQUALITY(e4.extent_int(3),4);

    TEST_EQUALITY(e5.extent_int(0),1);
    TEST_EQUALITY(e5.extent_int(1),2);
    TEST_EQUALITY(e5.extent_int(2),3);
    TEST_EQUALITY(e5.extent_int(3),4);
    TEST_EQUALITY(e5.extent_int(4),5);

    TEST_EQUALITY(e6.extent_int(0),1);
    TEST_EQUALITY(e6.extent_int(1),2);
    TEST_EQUALITY(e6.extent_int(2),3);
    TEST_EQUALITY(e6.extent_int(3),4);
    TEST_EQUALITY(e6.extent_int(4),5);
    TEST_EQUALITY(e6.extent_int(5),6);

    TEST_EQUALITY(e7.extent_int(0),1);
    TEST_EQUALITY(e7.extent_int(1),2);
    TEST_EQUALITY(e7.extent_int(2),3);
    TEST_EQUALITY(e7.extent_int(3),4);
    TEST_EQUALITY(e7.extent_int(4),5);
    TEST_EQUALITY(e7.extent_int(5),6);
    TEST_EQUALITY(e7.extent_int(6),7);

    TEST_EQUALITY(e8.extent_int(0),1);
    TEST_EQUALITY(e8.extent_int(1),2);
    TEST_EQUALITY(e8.extent_int(2),3);
    TEST_EQUALITY(e8.extent_int(3),4);
    TEST_EQUALITY(e8.extent_int(4),5);
    TEST_EQUALITY(e8.extent_int(5),6);
    TEST_EQUALITY(e8.extent_int(6),7);
    TEST_EQUALITY(e8.extent_int(7),8);
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // setExtents() from base class implementation
  {
    DataLayout& b1 = e1;
    b1.setExtents(10);
    TEST_EQUALITY(b1.extent(0),10);

    DataLayout& b2 = e2;
    b2.setExtents(10,20);
    TEST_EQUALITY(b2.extent(0),10);
    TEST_EQUALITY(b2.extent(1),20);

    DataLayout& b3 = e3;
    b3.setExtents(10,20,30);
    TEST_EQUALITY(b3.extent(0),10);
    TEST_EQUALITY(b3.extent(1),20);
    TEST_EQUALITY(b3.extent(2),30);

    DataLayout& b4 = e4;
    b4.setExtents(10,20,30,40);
    TEST_EQUALITY(b4.extent(0),10);
    TEST_EQUALITY(b4.extent(1),20);
    TEST_EQUALITY(b4.extent(2),30);
    TEST_EQUALITY(b4.extent(3),40);

    DataLayout& b5 = e5;
    b5.setExtents(10,20,30,40,50);
    TEST_EQUALITY(b5.extent(0),10);
    TEST_EQUALITY(b5.extent(1),20);
    TEST_EQUALITY(b5.extent(2),30);
    TEST_EQUALITY(b5.extent(3),40);
    TEST_EQUALITY(b5.extent(4),50);

    DataLayout& b6 = e6;
    b6.setExtents(10,20,30,40,50,60);
    TEST_EQUALITY(b6.extent(0),10);
    TEST_EQUALITY(b6.extent(1),20);
    TEST_EQUALITY(b6.extent(2),30);
    TEST_EQUALITY(b6.extent(3),40);
    TEST_EQUALITY(b6.extent(4),50);
    TEST_EQUALITY(b6.extent(5),60);

    DataLayout& b7 = e7;
    b7.setExtents(10,20,30,40,50,60,70);
    TEST_EQUALITY(b7.extent(0),10);
    TEST_EQUALITY(b7.extent(1),20);
    TEST_EQUALITY(b7.extent(2),30);
    TEST_EQUALITY(b7.extent(3),40);
    TEST_EQUALITY(b7.extent(4),50);
    TEST_EQUALITY(b7.extent(5),60);
    TEST_EQUALITY(b7.extent(6),70);

    DataLayout& b8 = e8;
    b8.setExtents(10,20,30,40,50,60,70,80);
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
    cout << "\n" << e1.name(0) << std::endl;

    TEST_EQUALITY(e1.name(0),"EXT0");

    TEST_EQUALITY(e2.name(0),"EXT0");
    TEST_EQUALITY(e2.name(1),"EXT1");

    TEST_EQUALITY(e3.name(0),"EXT0");
    TEST_EQUALITY(e3.name(1),"EXT1");
    TEST_EQUALITY(e3.name(2),"EXT2");

    TEST_EQUALITY(e4.name(0),"EXT0");
    TEST_EQUALITY(e4.name(1),"EXT1");
    TEST_EQUALITY(e4.name(2),"EXT2");
    TEST_EQUALITY(e4.name(3),"EXT3");

    TEST_EQUALITY(e5.name(0),"EXT0");
    TEST_EQUALITY(e5.name(1),"EXT1");
    TEST_EQUALITY(e5.name(2),"EXT2");
    TEST_EQUALITY(e5.name(3),"EXT3");
    TEST_EQUALITY(e5.name(4),"EXT4");

    TEST_EQUALITY(e6.name(0),"EXT0");
    TEST_EQUALITY(e6.name(1),"EXT1");
    TEST_EQUALITY(e6.name(2),"EXT2");
    TEST_EQUALITY(e6.name(3),"EXT3");
    TEST_EQUALITY(e6.name(4),"EXT4");
    TEST_EQUALITY(e6.name(5),"EXT5");

    TEST_EQUALITY(e7.name(0),"EXT0");
    TEST_EQUALITY(e7.name(1),"EXT1");
    TEST_EQUALITY(e7.name(2),"EXT2");
    TEST_EQUALITY(e7.name(3),"EXT3");
    TEST_EQUALITY(e7.name(4),"EXT4");
    TEST_EQUALITY(e7.name(5),"EXT5");
    TEST_EQUALITY(e7.name(6),"EXT6");

    TEST_EQUALITY(e8.name(0),"EXT0");
    TEST_EQUALITY(e8.name(1),"EXT1");
    TEST_EQUALITY(e8.name(2),"EXT2");
    TEST_EQUALITY(e8.name(3),"EXT3");
    TEST_EQUALITY(e8.name(4),"EXT4");
    TEST_EQUALITY(e8.name(5),"EXT5");
    TEST_EQUALITY(e8.name(6),"EXT6");
    TEST_EQUALITY(e8.name(7),"EXT7");
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // ostream
  {
    cout << "Testing ostream...";
    ostringstream output;
    output << dl8 << endl;
    output << dl7 << endl;
    output << dl6 << endl;
    output << dl5 << endl;
    output << dl4 << endl;
    output << dl3 << endl;
    output << dl2 << endl;
    output << dl1 << endl;
    output << e8 << endl;
    output << e7 << endl;
    output << e6 << endl;
    output << e5 << endl;
    output << e4 << endl;
    output << e3 << endl;
    output << e2 << endl;
    output << e1 << endl;
    output << n_mat << endl;
    cout << "...passed:\n" << output.str() << endl;
  }

  // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  // Arbitrary layout
  TEST_ASSERT(e7.kokkosLayout() == PHX::DataLayout::KokkosLayoutType::Default);
  e7.setKokkosLayout(PHX::DataLayout::KokkosLayoutType::Right);
  TEST_ASSERT(e7.kokkosLayout() == PHX::DataLayout::KokkosLayoutType::Right);
  e7.setKokkosLayout(PHX::DataLayout::KokkosLayoutType::Left);
  TEST_ASSERT(e7.kokkosLayout() == PHX::DataLayout::KokkosLayoutType::Left);
}
