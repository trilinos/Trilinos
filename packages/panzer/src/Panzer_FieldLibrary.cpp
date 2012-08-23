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

#include "Panzer_FieldLibrary.hpp"

namespace panzer {

Teuchos::RCP<const panzer::PureBasis> FieldLayoutLibrary::lookupBasis(const std::string & fieldName) const
{
   Teuchos::RCP<panzer::BasisIRLayout> layout = lookupLayout(fieldName);

   return layout->getBasis();
}

void FieldLayoutLibrary::uniqueBases(std::vector<Teuchos::RCP<const panzer::PureBasis> > & bases) const
{
   bases.clear();
   
   // simply loop over map of basis name to pointers and add them to the vector
   std::map<std::string,Teuchos::RCP<const panzer::PureBasis> >::const_iterator itr;
   for(itr=basisNameToPointer_.begin();itr!=basisNameToPointer_.end();++itr) 
      bases.push_back(itr->second);
}

void FieldLayoutLibrary::addFieldAndLayout(const std::string & fieldName,
                                      const Teuchos::RCP<panzer::BasisIRLayout> & layout)
{
   fieldToLayout_[fieldName] = layout; 
   basisNameToPointer_[layout->getBasis()->name()] = layout->getBasis();
}

Teuchos::RCP<panzer::BasisIRLayout> FieldLayoutLibrary::lookupLayout(const std::string & fieldName) const
{
   typedef std::map<std::string,Teuchos::RCP<panzer::BasisIRLayout> > Map;
   Map::const_iterator itr = fieldToLayout_.find(fieldName);
   if(itr!=fieldToLayout_.end())
      return itr->second;
 
   return Teuchos::null;
}

void FieldLayoutLibrary::print(std::ostream & os) const
{
   typedef std::map<std::string,Teuchos::RCP<panzer::BasisIRLayout> > Map;
   
   for(Map::const_iterator itr=fieldToLayout_.begin();itr!=fieldToLayout_.end();++itr) {
      std::string fieldName = itr->first; 
      Teuchos::RCP<BasisIRLayout> basis = itr->second;

      os << "\"" << fieldName << "\"" << " {" << basis->name() 
         << "(dim=" << basis->getDimension() 
         << ",cells=" << basis->getNumCells() 
         << ",points=" << basis->getNumPoints() << ")} ";
   }
}

void FieldLayoutLibrary::basisPairs(std::vector<std::pair<std::string,Teuchos::RCP<const panzer::PureBasis> > > & bases) const
{
   typedef std::map<std::string,Teuchos::RCP<panzer::BasisIRLayout> > Map;
   bases.clear();
   
   for(Map::const_iterator itr=fieldToLayout_.begin();itr!=fieldToLayout_.end();++itr) {
      std::string fieldName = itr->first; 
      Teuchos::RCP<const PureBasis> basis = itr->second->getBasis();

      bases.push_back(std::make_pair(fieldName,basis));
   }
}

///////////////////////////////////////////////////////////////////

Teuchos::RCP<const panzer::PureBasis> FieldLibrary::lookupBasis(const std::string & fieldName) const
{
   typedef std::map<std::string,Teuchos::RCP<panzer::PureBasis> > Map;
   Map::const_iterator itr = fieldToBasis_.find(fieldName);
   if(itr!=fieldToBasis_.end())
      return itr->second;
 
   return Teuchos::null;
}

void FieldLibrary::uniqueBases(std::vector<Teuchos::RCP<const panzer::PureBasis> > & bases) const
{
   bases.clear();
   
   // simply loop over map of basis name to pointers and add them to the vector
   std::map<std::string,Teuchos::RCP<const panzer::PureBasis> >::const_iterator itr;
   for(itr=basisNameToPointer_.begin();itr!=basisNameToPointer_.end();++itr) 
      bases.push_back(itr->second);
}

void FieldLibrary::addFieldAndBasis(const std::string & fieldName,
                                    const Teuchos::RCP<panzer::PureBasis> & basis)
{
   fieldToBasis_[fieldName] = basis;
   basisNameToPointer_[basis->name()] = basis;
}

Teuchos::RCP<const FieldLayoutLibrary> FieldLibrary::buildFieldLayoutLibrary(panzer::PointRule & ir) const
{
   typedef std::map<std::string,Teuchos::RCP<panzer::PureBasis> > Map;

   Teuchos::RCP<FieldLayoutLibrary> layoutLibrary = Teuchos::rcp(new FieldLayoutLibrary);

   // loop over each member of the map, create a new FieldLayout and addit to layout library
   for(Map::const_iterator itr=fieldToBasis_.begin();itr!=fieldToBasis_.end();++itr) {
      Teuchos::RCP<BasisIRLayout> layout = Teuchos::rcp(new BasisIRLayout(itr->second,ir));
      layoutLibrary->addFieldAndLayout(itr->first,layout);
   }
 
   return layoutLibrary;
}

void FieldLibrary::print(std::ostream & os) const
{
   typedef std::map<std::string,Teuchos::RCP<panzer::PureBasis> > Map;
   
   for(Map::const_iterator itr=fieldToBasis_.begin();itr!=fieldToBasis_.end();++itr) {
      std::string fieldName = itr->first; 
      Teuchos::RCP<PureBasis> basis = itr->second;

      os << "\"" << fieldName << "\"" << " {" << basis->name() 
         << "(dim=" << basis->getDimension() 
         << ",cells=" << basis->getNumCells() << ") ";
   }
}

//! Get vector of unique bases contained in this field library
void FieldLibrary::basisPairs(std::vector<std::pair<std::string,Teuchos::RCP<const panzer::PureBasis> > > & bases) const
{
   typedef std::map<std::string,Teuchos::RCP<panzer::PureBasis> > Map;
   bases.clear();
   
   for(Map::const_iterator itr=fieldToBasis_.begin();itr!=fieldToBasis_.end();++itr) {
      std::string fieldName = itr->first; 
      Teuchos::RCP<PureBasis> basis = itr->second;

      bases.push_back(std::make_pair(fieldName,basis.getConst()));
   }
}

}
