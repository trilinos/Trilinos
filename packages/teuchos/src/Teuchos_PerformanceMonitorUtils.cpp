// @HEADER
// ***********************************************************************
//
//                    Teuchos: Common Tools Package
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ***********************************************************************
// @HEADER

#include "Teuchos_PerformanceMonitorUtils.hpp"
#include "Teuchos_MPIContainerComm.hpp"
#include "Teuchos_ConfigDefs.hpp"
using namespace Teuchos;


void PerformanceMonitorUtils::synchNames(const MPIComm& comm,
                                        const Array<std::string>& localNames,
                                        Array<std::string>& allNames)
{
  if (comm.getNProc() > 1)
    {
      /* gather names of counters from all processors */
      int root = 0;
      std::set<std::string> nameSet;
      Array<Array<std::string> > namesForAllProcs;
      MPIContainerComm<std::string>::gatherv(localNames, namesForAllProcs, 
                                        root, comm);
      
      /* on the root processor, compile the set union of all names */
      if (comm.getRank()==0)
        {
          for (Array<Array<std::string> >::size_type p=0; p<namesForAllProcs.size(); p++)
            {
              for (Array<std::string>::size_type i=0; i<namesForAllProcs[p].size(); i++)
                {
                  nameSet.insert(namesForAllProcs[p][i]);
                }
            }
        }
      
      /* convert the set to an array so we can send it out by MPI */
      allNames.resize(0);
      for (std::set<std::string>::const_iterator i=nameSet.begin(); i!=nameSet.end(); i++)
        {
          allNames.append(*i);
        }
      /* broadcast the union of all names to all processors */
      MPIContainerComm<std::string>::bcast(allNames, root, comm); 
    }
  else
    {
      allNames = localNames;
    }
}

void PerformanceMonitorUtils
::synchValues(const MPIComm& comm,
              const Array<std::string>& localNames,
              const Array<Array<double> >& localValues,
              Array<std::string>& allNames,
              Array<Array<double> >& allValues)
{
  std::map<std::string, Array<double> > localNameToValMap;

  for (Array<std::string>::size_type i=0; i<localNames.size(); i++)
    {
      Array<double> tmp(localValues.size());
      for (Array<Array<double> >::size_type j=0; j<localValues.size(); j++)
        {
          tmp[j] = localValues[j][i];
        }
      localNameToValMap[localNames[i]] = tmp;
    }

  synchNames(comm, localNames, allNames);
  
  allValues.resize(localValues.size());
  for (Array<Array<double> >::size_type i=0; i<allValues.size(); i++)
    {
      allValues[i].resize(allNames.size());
    }

  for (Array<std::string>::size_type i=0; i<allNames.size(); i++)
    {
      const std::string& name = allNames[i];
      if (localNameToValMap.find(name) != localNameToValMap.end()) 
        {
          const Array<double>& tmp = localNameToValMap[name];
          for (Array<double>::size_type j=0; j<tmp.size(); j++)
            {
              allValues[j][i] = tmp[j];
            }
        }
      else 
        {
          for (Array<Array<double> >::size_type j=0; j<allValues.size(); j++)
            {
              allValues[j][i] = 0.0;
            }
        }
    }
}


void  PerformanceMonitorUtils::reduce(const MPIComm& comm,
                                     const EMetricReduction& reductionType,
                                     const Array<double>& localVals,
                                     Array<double>& reducedVals)
{
  /* if we're asking for local values, do nothing but copy the local array
   * to the reduced array. */
  if (comm.getNProc()==1 || reductionType==ELocal)
    {
      reducedVals = localVals;
      return;
    }

  /* If we're to this point we must do a reduction */
  reducedVals.resize(localVals.size());

  int op = MPIComm::SUM;
  if (reductionType==EMax) op = MPIComm::MAX;
  if (reductionType==EMin) op = MPIComm::MIN;
  
  int sendCount = localVals.size();

  if (sendCount==0) return;

  double* sendBuf = const_cast<double*>(&localVals[0]);
  double* recvBuf = const_cast<double*>(&reducedVals[0]);

  comm.allReduce( (void*) sendBuf, (void*) recvBuf, sendCount, MPIComm::DOUBLE, op);

  if (reductionType==EAvg)
    {
      for (Array<double>::size_type i=0; i<reducedVals.size(); i++)
        {
          reducedVals[i] /= ((double) comm.getNProc());
        }
    }
}




