
//@HEADER
// ************************************************************************
// 
//               ShyLU: Hybrid preconditioner package
//                 Copyright 2012 Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Clark R. Dohrmann (crdohrm@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#ifndef OPERATORSTANDARDBDDC_H
#define OPERATORSTANDARDBDDC_H
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string.h>

#include "shylu_OperatorBase.hpp"

namespace bddc {
  
template <class SX> class OperatorStandard : 
  public OperatorBase<SX>
{
public:

  OperatorStandard(const int numNode,
		   const int* nodeBegin,
		   const std::vector< std::vector<int> > & subNodes,
		   const std::vector< std::vector<int> > & rowBeginSub,
		   const std::vector< std::vector<int> > & columnsSub,
		   const std::vector< std::vector<SX> > & valuesSub) :
    m_numNode(numNode),
    m_nodeBegin(nodeBegin),
    m_subNodes(subNodes),
    m_rowBeginSub(rowBeginSub),
    m_columnsSub(columnsSub),
    m_valuesSub(valuesSub)
  {
    determineSubDofs();
  }

  ~OperatorStandard()
  {
  }

  void Apply(const SX* x,
	     SX* Ax)
  {
    const int numSub = m_subNodes.size();
    const int numRows = m_nodeBegin[m_numNode];
    memset(Ax, 0, numRows*sizeof(SX));
    for (int m=0; m<numSub; m++) {
      const std::vector<int> & rowBegin = m_rowBeginSub[m];
      const std::vector<int> & columns = m_columnsSub[m];
      const std::vector<SX> & values = m_valuesSub[m];
      const std::vector<int> & subDofs = m_subDofs[m];
      const int numDofs = subDofs.size();
      for (int i=0; i<numDofs; i++) {
	SX sum(0);
	for (int j=rowBegin[i]; j<rowBegin[i+1]; j++) {
	  const int col = subDofs[columns[j]];
	  sum += values[j]*x[col];
	}
	const int row = subDofs[i];
	Ax[row] += sum;
      }
    }
  }

private:
  const int m_numNode{0}, *m_nodeBegin{nullptr};
  const std::vector< std::vector<int> > & m_subNodes, m_rowBeginSub, m_columnsSub;
  const std::vector< std::vector<SX> > & m_valuesSub;
  std::vector< std::vector<int> > m_subDofs;

  void determineSubDofs()
  {
    const int numSub = m_subNodes.size();
    m_subDofs.resize(numSub);
    for (int k=0; k<numSub; k++) {
      const std::vector<int> subNodes = m_subNodes[k];
      const int numNodes = subNodes.size();
      int numDofs(0);
      for (int i=0; i<numNodes; i++) {
	const int node = subNodes[i];
	numDofs += m_nodeBegin[node+1] - m_nodeBegin[node];
      }
      std::vector<int> & subDofs = m_subDofs[k];
      subDofs.resize(numDofs);
      numDofs = 0;
      for (int i=0; i<numNodes; i++) {
	const int node = subNodes[i];
	for (int j=m_nodeBegin[node]; j<m_nodeBegin[node+1]; j++) {
	  subDofs[numDofs++] = j;
	}
      }
    }
  }

};
  
} // namespace bddc

#endif // OPERATORSTANDARDBDDC_H
  
