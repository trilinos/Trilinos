// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef __Panzer_WorksetDescriptor_hpp__
#define __Panzer_WorksetDescriptor_hpp__

#include <string>
#include <ostream>

#include <boost/functional/hash.hpp>

namespace panzer {

/** Class provides a simple description of the types of worksets
  * that need to be constructed and used. This is a generalization of
  * using strings and pairs of string to represent the element blocks
  * and sidesets. It is primarily used in specifying the "domain" of the 
  * assembly algorithm, that is which elements will be used in the assembly
  * process.
  */
class WorksetDescriptor {
public:

  /** Constructor specifying a lone element block.
    *
    * \param[in] eBlock Name of the element block
    */
  WorksetDescriptor(const std::string & eBlock)
    : elementBlock_(eBlock), 
      useSideset_(false), 
      sideset_(""), 
      sideAssembly_(false)
  { }

  /** Constructor that defines a side set. Note that the
    * specified sideset must be a non-empty string. 
    *
    * \param[in] eBlock Element block that includes the side
    * \param[in] eBlock Side set that is being used
    * \param[in] sideAssembly Are integration rules and 
    *                         basis functions evaluated on the
    *                         side or on the volume of the element.
    */
  WorksetDescriptor(const std::string & eBlock,
                    const std::string & sideset,
                    bool sideAssembly)
    : elementBlock_(eBlock), 
      useSideset_(true), 
      sideset_(sideset), 
      sideAssembly_(sideAssembly)
  {
    TEUCHOS_TEST_FOR_EXCEPTION(sideset_=="",std::runtime_error,
                               "WorksetDescriptor constr: Side set name must be non-empty!");
  }

  //! Copy constructor
  WorksetDescriptor(const WorksetDescriptor & src)
    : elementBlock_(src.elementBlock_), 
      useSideset_(src.useSideset_), 
      sideset_(src.sideset_), 
      sideAssembly_(src.sideAssembly_)
  {}
 
  //! Get element block 
  std::string getElementBlock() const
  { return elementBlock_; }

  //! Get the side set
  std::string getSideset() const
  { return sideset_; }

  //! Expects side set assembly on volume
  bool sideAssembly() const
  { return sideAssembly_; }
 
  //! This descriptor is for a side set.
  bool useSideset() const
  { return useSideset_; } 

private:

  //! Element block, required to be non-empty
  std::string elementBlock_;
 
  //! Use the side set information or not
  bool useSideset_;
  
  //! Side set, must be non-empty if <code>useSideset_</code> is true
  std::string sideset_;

  /** This indicates if side quadrature rules are constructed
    * or volume rules are constructued. Ignored if useSideset_
    * is false.
    */
  bool sideAssembly_;
};

//! Equality operation for use with hash tables and maps
inline bool operator==(const WorksetDescriptor & a,const WorksetDescriptor & b)
{
  if(a.useSideset())
    // if side set is in use, check all fields
    return    a.getElementBlock()==b.getElementBlock()
           && a.getSideset()==b.getSideset()
           && a.sideAssembly()==b.sideAssembly()
           && a.useSideset()==b.useSideset();
  else
    // otherwise check that both descriptor don't use side sets
    // and check the element block (the remaining fields are allowed
    // to be unset)
    return    a.getElementBlock()==b.getElementBlock()
           && a.useSideset()==b.useSideset();
}

//! Hash function that satisifies the boost hash interface
inline std::size_t hash_value(const WorksetDescriptor & wd)
{
  std::size_t seed = 0;

  boost::hash_combine(seed,wd.getElementBlock());
  if(wd.useSideset()) {
    // optionally hash on side set and side assembly
    boost::hash_combine(seed,wd.getSideset());
    boost::hash_combine(seed,wd.sideAssembly());
  }
 
  return seed;
}

//! I/O utility
inline std::ostream & operator<<(std::ostream & os,const WorksetDescriptor & wd)
{
  if(wd.useSideset())
    os << "Side descriptor: "
       << "eblock = \"" << wd.getElementBlock() << "\", "
       << "ss = \"" << wd.getSideset() << "\", "
       << "side assembly = " << (wd.sideAssembly() ? "on" : "off");
  else
    os << "Block descriptor: "
       << "eblock = \"" << wd.getElementBlock() << "\"";

  return os;
}

/** Builds a descriptor specifying an element block.
  */
inline WorksetDescriptor blockDescriptor(const std::string & eBlock)
{ return WorksetDescriptor(eBlock); }

/** Builds a descriptor specifying a sideset, specify surface terms.
  */
inline WorksetDescriptor sidesetDescriptor(const std::string & eBlock,const std::string & sideset)
{ return WorksetDescriptor(eBlock,sideset,false); }

/** Builds a descriptor specifying a sideset, however specify volumetric terms not
  * surface terms.
  */
inline WorksetDescriptor sidesetVolumeDescriptor(const std::string & eBlock,const std::string & sideset)
{ return WorksetDescriptor(eBlock,sideset,true); }

}

#endif
