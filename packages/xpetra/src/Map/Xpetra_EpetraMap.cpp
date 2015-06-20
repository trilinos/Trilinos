// @HEADER
//
// ***********************************************************************
//
//             Xpetra: A linear algebra interface package
//                  Copyright 2012 Sandia Corporation
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
// Questions? Contact
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include "Xpetra_ConfigDefs.hpp"

#ifdef HAVE_XPETRA_EPETRA

#include "Xpetra_EpetraMap.hpp"

namespace Xpetra {

  template<class GlobalOrdinal>
  const Epetra_Map & toEpetra(const Map<int,GlobalOrdinal> &map) {
    // TODO: throw exception
    const EpetraMapT<GlobalOrdinal> & epetraMap = dynamic_cast<const EpetraMapT<GlobalOrdinal> &>(*map.getMap());
    return epetraMap.getEpetra_Map();
  }

  template<class GlobalOrdinal>
  const Epetra_Map & toEpetra(const RCP< const Map<int,GlobalOrdinal> > &map) {
    XPETRA_RCP_DYNAMIC_CAST(const EpetraMapT<GlobalOrdinal>, map->getMap(), epetraMap, "toEpetra");
    return epetraMap->getEpetra_Map();
  }


//  template<class GlobalOrdinal>
//   const RCP<const Map<int,GlobalOrdinal> > toXpetra(const RCP<const Epetra_Map>& map) {
//     return rcp(new EpetraMapT(map));
//   }

  template<class GlobalOrdinal>
  const RCP< const Map<int, GlobalOrdinal> > toXpetra(const Epetra_BlockMap &map) {
    RCP<const Epetra_BlockMap> m = rcp(new Epetra_BlockMap(map));
    return rcp( new EpetraMapT<GlobalOrdinal>(m) );
  }
  //

#ifndef XPETRA_EPETRA_NO_32BIT_GLOBAL_INDICES
template const Epetra_Map & toEpetra<int>(const Map<int,int> &map);
template const Epetra_Map & toEpetra<int>(const RCP< const Map<int, int> > &);
template const RCP< const Map<int, int> > toXpetra<int>(const Epetra_BlockMap &map);
#endif

#ifndef XPETRA_EPETRA_NO_64BIT_GLOBAL_INDICES
template const Epetra_Map & toEpetra<long long>(const Map<int,long long> &map);
template const Epetra_Map & toEpetra<long long>(const RCP< const Map<int, long long> > &);
template const RCP< const Map<int, long long> > toXpetra<long long>(const Epetra_BlockMap &map);
#endif

}


#endif
