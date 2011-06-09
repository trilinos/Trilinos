/*
 * Copyright(C) 2010 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
 * certain rights in this software
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 * 
 *     * Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 * 
 *     * Redistributions in binary form must reproduce the above
 *       copyright notice, this list of conditions and the following
 *       disclaimer in the documentation and/or other materials provided
 *       with the distribution.
 * 
 *     * Neither the name of Sandia Corporation nor the names of its
 *       contributors may be used to endorse or promote products derived
 *       from this software without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 */
#ifndef SEACAS_Variables_H
#define SEACAS_Variables_H
#include <smart_assert.h>
#include <string>
#include <cstring>
#include <vector>
#include <ObjectType.h>

namespace Excn {
  enum InOut {IN = 1, OUT = 2};
  
  typedef std::vector<int>  IntVector;
  
  struct Variables {
    Variables(ObjectType otype) : objectType(otype), outputCount(0), addProcessorId(false)
    {
      SMART_ASSERT(otype == EBLK || otype == NSET || otype == SSET || otype == NODE || otype == GLOBAL);
    }
    
    int count(InOut in_out = IN) const
    {
      int ret_val = 0;
      switch (in_out) {
      case IN:
	ret_val = index_.size() - (addProcessorId ? 1 : 0);
	break;
      case OUT:
	ret_val = outputCount;
	break;
      }
      return ret_val;
    }

    const char *label() const
    {
      switch(objectType) {
      case EBLK:
	return "element";
      case NSET:
	return "nodeset";
      case GLOBAL:
	return "global";
      case NODE:
	return "nodal";
      case SSET:
	return "sideset";
      default:
	return "UNKNOWN";
      }
    }

    ex_entity_type type() const
    {
      switch(objectType) {
      case EBLK:
	return EX_ELEM_BLOCK;
      case NSET:
	return EX_NODE_SET;
      case SSET:
	return EX_SIDE_SET;
      case NODE:
	return EX_NODAL;
      case GLOBAL:
	return EX_GLOBAL;
      default:
	return EX_INVALID;
      }
    }

    bool add_processor_id() const {return addProcessorId;}
    
    ObjectType objectType;
    int outputCount;
    bool addProcessorId;
    IntVector index_;
    std::string type_;
  };
}

#endif
