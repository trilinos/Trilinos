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


#include "Cell.hpp"

//**********************************************************************
MyCell::MyCell() :
  m_phi_mem(Teuchos::arcp<double>(4*4)),
  m_grad_phi_mem(Teuchos::arcp<double>(4*4*3)),
  m_phi(&(m_phi_mem[0]),4,4),
  m_grad_phi(&(m_grad_phi_mem[0]),4,4,3)
{ 
  using namespace Teuchos;

  for (ArrayRCP<double>::Ordinal i = 0; i < m_phi_mem.size(); ++i)
    m_phi_mem[i] = 0.25;

  for (ArrayRCP<double>::Ordinal i = 0; i < m_grad_phi_mem.size(); ++i)
    m_grad_phi_mem[i] = 0.25;
}

//**********************************************************************
shards::Array<double,shards::NaturalOrder,Node,Dim>& 
MyCell::getNodeCoordinates()
{
  return m_coords;
}

//**********************************************************************
shards::Array<double,shards::NaturalOrder,QuadPoint,Node>& 
MyCell::getBasisFunctions()
{
  return m_phi;
}

//**********************************************************************
shards::Array<double,shards::NaturalOrder,QuadPoint,Node,Dim>& 
MyCell::getBasisFunctionGradients()
{
  return m_grad_phi;
}

//**********************************************************************
std::size_t MyCell::localIndex()
{
  return local_index_;
}

//**********************************************************************
void MyCell::setLocalIndex(std::size_t index)
{
  local_index_ = index;
}

//**********************************************************************
