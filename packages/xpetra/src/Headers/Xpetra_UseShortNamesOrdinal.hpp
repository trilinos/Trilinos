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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
// Get rid of template parameters

// New definition of types using the types LocalOrdinal, GlobalOrdinal, Node, LocalMatOps of the current context.
#ifdef XPETRA_MAP_SHORT
typedef Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node> Map;
#endif

#ifdef XPETRA_MAPFACTORY_SHORT
typedef Xpetra::MapFactory<LocalOrdinal, GlobalOrdinal, Node> MapFactory;
#endif

#ifdef XPETRA_CRSGRAPH_SHORT
typedef Xpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsGraph;
#endif

#ifdef XPETRA_CRSGRAPHFACTORY_SHORT
typedef Xpetra::CrsGraphFactory<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> CrsGraphFactory;
#endif

#ifdef XPETRA_VECTOR_SHORT
typedef Xpetra::Vector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> LocalOrdinalVector;
typedef Xpetra::Vector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> GlobalOrdinalVector;
typedef LocalOrdinalVector  LOVector;
typedef GlobalOrdinalVector GOVector;
#endif

#ifdef XPETRA_MULTIVECTOR_SHORT
typedef Xpetra::MultiVector<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> LocalOrdinalMultiVector;
typedef Xpetra::MultiVector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> GlobalOrdinalMultiVector;
typedef LocalOrdinalMultiVector  LOMultiVector;
typedef GlobalOrdinalMultiVector GOMultiVector;
#endif

#ifdef XPETRA_VECTORFACTORY_SHORT
typedef Xpetra::VectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> LocalOrdinalVectorFactory;
typedef Xpetra::VectorFactory<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> GlobalOrdinalVectorFactory;
typedef LocalOrdinalVectorFactory  LOVectorFactory;
typedef GlobalOrdinalVectorFactory GOVectorFactory;
#endif

#ifdef XPETRA_MULTIVECTORFACTORY_SHORT
typedef Xpetra::MultiVectorFactory<LocalOrdinal, LocalOrdinal, GlobalOrdinal, Node> LocalOrdinalMultiVectorFactory;
typedef Xpetra::MultiVectorFactory<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node> GlobalOrdinalMultiVectorFactory;
typedef LocalOrdinalMultiVectorFactory  LOMultiVectorFactory;
typedef GlobalOrdinalMultiVectorFactory GOMultiVectorFactory;
#endif

#ifdef XPETRA_IMPORT_SHORT
typedef Xpetra::Import<LocalOrdinal, GlobalOrdinal, Node> Import;
#endif

#ifdef XPETRA_EXPORT_SHORT
typedef Xpetra::Export<LocalOrdinal, GlobalOrdinal, Node> Export;
#endif

#ifdef XPETRA_IMPORTFACTORY_SHORT
typedef Xpetra::ImportFactory<LocalOrdinal, GlobalOrdinal, Node> ImportFactory;
#endif

#ifdef XPETRA_EXPORTFACTORY_SHORT
typedef Xpetra::ExportFactory<LocalOrdinal, GlobalOrdinal, Node> ExportFactory;
#endif

#ifdef XPETRA_TPETRAMAP_SHORT
typedef Xpetra::TpetraMap<LocalOrdinal, GlobalOrdinal, Node> TpetraMap;
#endif

#ifdef XPETRA_TPETRACRSGRAPH_SHORT
typedef Xpetra::TpetraCrsGraph<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> TpetraCrsGraph;
#endif

#ifdef XPETRA_STRIDEDMAP_SHORT
typedef Xpetra::StridedMap<LocalOrdinal, GlobalOrdinal, Node> StridedMap;
#endif

#ifdef XPETRA_STRIDEDMAPFACTORY_SHORT
typedef Xpetra::StridedMapFactory<LocalOrdinal, GlobalOrdinal, Node> StridedMapFactory;
#endif

// Note: There is no #ifndef/#define/#end in this header file because it can be included more than once (it can be included in methods templated by Scalar, LocalOrdinal, GlobalOrdinal, Node).

// TODO: add namespace {} for shortcut types

// Define convenient shortcut for data types
typedef LocalOrdinal  LO;
typedef GlobalOrdinal GO;
typedef Node          NO;
typedef LocalMatOps   LMO;

// TODO: do the same for Epetra object (problem of namespace)
