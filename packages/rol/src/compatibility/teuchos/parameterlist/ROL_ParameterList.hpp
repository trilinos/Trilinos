// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#pragma once
#include "ROL_Ptr.hpp"
#include "ROL_ParameterList.hpp"
#include "Teuchos_XMLParameterListCoreHelpers.hpp"

namespace ROL {
using ParameterList = Teuchos::ParameterList;

inline
void readParametersFromXml( const std::string &filename,
                            ParameterList& parlist ) {
  Teuchos::Ptr<ParameterList> p{&parlist};
  Teuchos::updateParametersFromXmlFile( filename, p );
}

  inline void writeParameterListToXmlFile( ParameterList& parlist,
                                           const std::string& filename ) {
    Teuchos::writeParameterListToXmlFile(parlist, filename);
  }

  inline void updateParametersFromXmlFile(const std::string& filename ,
                                          ParameterList& parlist) {
    Teuchos::Ptr<ParameterList> p{&parlist};
    Teuchos::updateParametersFromXmlFile( filename, p );
  }



inline
ROL::Ptr<ParameterList> getParametersFromXmlFile( const std::string& filename )
  {
    auto parlist = ROL::makePtr<ParameterList>();
    Teuchos::Ptr<ParameterList> p{&*parlist};
    Teuchos::updateParametersFromXmlFile( filename, p );
    return parlist;
  }


template <class T>
inline std::vector<T> getArrayFromStringParameter( const ParameterList& parlist,
                                                  const std::string& name ) {
  auto teuchos_arr = Teuchos::getArrayFromStringParameter<T>( parlist, name );

    return teuchos_arr.toVector();
}


} // namespace ROL
