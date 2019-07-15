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


#include <map>

#include "Phalanx_FieldTag.hpp"
#include "Phalanx_FieldTag_Tag.hpp"
#include "Phalanx_FieldTag_STL_Functors.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"
#include "Phalanx_ExtentTraits.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_UnitTestHarness.hpp"

// From test/Utilities directory
#include "Traits.hpp"

PHX_EXTENT(Cell)
PHX_EXTENT(Node)
PHX_EXTENT(QP)
PHX_EXTENT(DIM)

TEUCHOS_UNIT_TEST(FieldTag,FieldTag)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  // Dummy data layouts (same size different name/type)
  RCP<DataLayout> node4 = rcp(new MDALayout<Cell,Node>(100,4));
  RCP<DataLayout> quad4 = rcp(new MDALayout<Cell,QP>(100,4));
  RCP<DataLayout> quad4grad = rcp(new MDALayout<Cell,QP,DIM>(100,4,4));
  
  // Allocate tags with same name but different topology
  RCP<FieldTag> rcp_nodal_density = rcp(new Tag<double>("density", node4));
  RCP<FieldTag> rcp_qp_density = rcp(new Tag<double>("density", quad4));
  RCP<FieldTag> rcp_grad_qp_density = 
    rcp(new Tag<double>("density", quad4grad));
  
  // Get references to field tags
  FieldTag& nodal_density = *rcp_nodal_density;
  FieldTag& qp_density = *rcp_qp_density;
  FieldTag& grad_qp_density = *rcp_grad_qp_density;
  
  // test ostream
  cout << "Printing field tags" << endl;
  cout << nodal_density << endl;
  cout << qp_density << endl;
  cout << grad_qp_density << endl;
  cout << endl;
  
  // test operator ==
  TEST_INEQUALITY(nodal_density, qp_density);
  
  // New constructor that should be same as nodal_density
  RCP<FieldTag> rcp_nodal_density_copy = 
    rcp(new Tag<double>("density", node4));
  FieldTag& nodal_density_copy = *rcp_nodal_density_copy;
  TEST_EQUALITY(nodal_density, nodal_density_copy);
  
  // Different layout (used to test scalar vs vector)
  TEST_INEQUALITY(qp_density, grad_qp_density);
  
  // test assignment operator =
  Tag<double> op_eq("Garbage", node4);
  TEST_INEQUALITY(op_eq, nodal_density);
  op_eq = dynamic_cast< Tag<double>& >(nodal_density);
  TEST_EQUALITY(op_eq, nodal_density);
  
  // name() accessor
  TEST_EQUALITY(nodal_density.name(),"density");
  
  // dataLayout() accessor
  const DataLayout& tmp = *node4;
  TEST_EQUALITY(nodal_density.dataLayout(), tmp);
  
  // clone()
  RCP<FieldTag> my_copy = rcp_nodal_density->clone();
  TEST_EQUALITY(*my_copy, *rcp_nodal_density);
  
  // Comparison for map key operations  
  map<RCP<FieldTag>, int, FTComp> my_map;
  my_map[rcp_nodal_density] = 0;
  my_map[rcp_qp_density] = 1;
  my_map[rcp_grad_qp_density] = 2;
  
  RCP<FieldTag> tmp_rcp_nodal_density = 
    rcp(new Tag<double>("density", node4));
  RCP<FieldTag> tmp_rcp_qp_density = 
    rcp(new Tag<double>("density", quad4));
  RCP<FieldTag> tmp_rcp_grad_qp_density = 
    rcp(new Tag<double>("density", quad4grad));
  
  // Make sure we can create a field tag and access matching map entry
  TEST_EQUALITY(my_map[tmp_rcp_nodal_density], 0);
  TEST_EQUALITY(my_map[tmp_rcp_qp_density], 1);
  TEST_EQUALITY(my_map[tmp_rcp_grad_qp_density], 2);
  
  // Comparison for vector search operations
  std::vector< Teuchos::RCP<FieldTag> > vec(0);
  vec.push_back(rcp_nodal_density);
  vec.push_back(rcp_qp_density);
  vec.push_back(rcp_grad_qp_density);
  
  { // reference version
    PHX::FTPredRef pred(*rcp_grad_qp_density);
    std::vector< Teuchos::RCP<FieldTag> >::iterator test = 
      std::find_if(vec.begin(), vec.end(), pred);
    
    TEST_EQUALITY(*(*test), *tmp_rcp_grad_qp_density);
  }
  
  { // RCP version
    PHX::FTPred pred(rcp_grad_qp_density);
    std::vector< Teuchos::RCP<FieldTag> >::iterator test = 
      std::find_if(vec.begin(), vec.end(), pred);
    
    TEST_EQUALITY(*(*test), *tmp_rcp_grad_qp_density);
  }
}

TEUCHOS_UNIT_TEST(FieldTag,ConstCorrectness)
{
  using namespace std;
  using namespace Teuchos;
  using namespace PHX;

  // Dummy data layouts (same size different name/type)
  RCP<DataLayout> node4 = rcp(new MDALayout<Cell,Node>(100,4));
  
  // Tags ignore const in the scalar type. Should be equal.
  Tag<double> nonconst_a("a", node4);
  Tag<const double> const_a("a", node4);
  TEST_ASSERT(nonconst_a == const_a);

  // Tags should ignore const - test all combinations of constructors
  Tag<double> nonconst_b(nonconst_a);      // nonconst from nonconst
  Tag<const double> const_c(const_a);      // const    from const
  Tag<const double> const_b(nonconst_a); // const    from nonconst
  Tag<double> nonconst_c(const_a);       // nonconst from const

  // This should fail to compile: uncomment to test static_assert
  // warning message
  //Tag<int> d(nonconst_a); // create int tag from double
}

TEUCHOS_UNIT_TEST(FieldTag,EmptyCtor)
{
  using namespace PHX;

  Tag<double> empty_tag;
  TEST_THROW(empty_tag.dataLayout(),std::logic_error);
  TEST_EQUALITY(empty_tag.name(),std::string("TAG_NAME_NOT_SET"));

  const Tag<double> const_empty_tag;
  TEST_THROW(const_empty_tag.dataLayout(),std::logic_error);
}
