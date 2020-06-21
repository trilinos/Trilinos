/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
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
// ************************************************************************
// @HEADER
*/
#include <numeric>
#include <Tpetra_Map.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_MultiVector.hpp>
#include "Tpetra_TestingUtilities.hpp"

#include <Teuchos_Array.hpp>
#include <Teuchos_CommHelpers.hpp>
#include <Teuchos_Tuple.hpp>

#include "TpetraCore_ETIHelperMacros.h"

namespace {

using Tpetra::TestingUtilities::getDefaultComm;

TEUCHOS_UNIT_TEST_TEMPLATE_3_DECL(Directory, AllMinGIDs, SC, LO, GO)
{
  /*
   * This issue is described in
   * [Issue 6987](https://github.com/trilinos/Trilinos/issues/6987)
   * and summarized here
   *
   * Map::myMinGID_ is the minimum GID held by a map on its owning process. For
   * a map that has no entries on a process, its value is std::numeric_limits<GO>::max()
   *
   * Directory::allMinGIDs_ array holds the minimum GID on each participating
   * process. This array is used, eg, in the directory's getEntriesImpl method.
   * For ranks whose maps are empty allMinGIDs_[rank] is std::numeric_limits<GO>::max()
   * and any work using it violates the assumptions of the search algorithm that assumes
   * the GIDs are uniformly distributed and monotonically increasing
   *
   * @rppawlo discovered this in a simple test wherein a map was empty on every
   * process except the last. That map is used to build a MultiVector and then
   * import that MultiVector to another MultiVector with a locally replicated
   * map. The result is failure to import the values from the original vector to
   * the locally replicated vector (because of comparisons with uninitialized
   * values in getEntriesImpl.
   */

  using map_type = Tpetra::Map<LO, GO>;
  using vector_type = Tpetra::MultiVector<SC, LO, GO>;
  using import_type = Tpetra::Import<LO, GO>;

  auto comm = getDefaultComm();
  const auto my_rank = comm->getRank();
  const auto num_procs = comm->getSize();

  const GO index_base = 0;
  const size_t num_vecs = 1;

  const size_t num_non_zero = 5;
  //const int non_empty_rank = num_procs - 1;
  const SC expected_value = Teuchos::as<SC>(5.0);
  std::vector<int> non_empty_ranks(num_procs);
  std::iota(non_empty_ranks.begin(), non_empty_ranks.end(), 0);
  for (auto && non_empty_rank : non_empty_ranks)
  {

    Teuchos::RCP<vector_type> vec1;
    Teuchos::RCP<map_type> map1;
    {
      const auto num_global_elements = Teuchos::as<Tpetra::global_size_t>(num_non_zero);
      const size_t num_local_elements = (my_rank == non_empty_rank) ? num_non_zero : 0;
      map1 = Teuchos::rcp(new map_type(num_global_elements, num_local_elements, index_base, comm));
      vec1 = Teuchos::rcp(new vector_type(map1, num_vecs, true));
      if (my_rank == non_empty_rank) vec1->putScalar(expected_value);
    }

    Teuchos::RCP<vector_type> vec2;
    Teuchos::RCP<const map_type> map2;
    {
      map2 = Tpetra::createLocalMap<LO,GO>(num_non_zero, comm);
      vec2 = Teuchos::rcp(new vector_type(map2, num_vecs, true));
    }

    auto import = import_type(map1, map2);
    vec2->doImport(*vec1, import, Tpetra::INSERT);

    auto data = vec2->getData(0);
    TEUCHOS_TEST_FOR_EXCEPTION(
      data.size() != num_non_zero,
      std::logic_error,
      "Vector data.size should be " << num_non_zero << " but is " << data.size()
    );

    std::vector<size_t> bad;
    for (size_t i=0; i<num_non_zero; i++)
    {
      if (data[0] != expected_value) bad.push_back(i);
    }

    if (bad.size() > 0)
    {
      out << "The following vector entries are incorrect after import:\n";
      for (auto && i : bad)
        out << "data[" << i << "] = " << data[i] << " != " << expected_value << "\n";
      TEST_ASSERT(false);
    }
  }
}

//
// INSTANTIATIONS
//

#define THIS_TEST_GROUP(SC, LO, GO) \
  TEUCHOS_UNIT_TEST_TEMPLATE_3_INSTANT(Directory, AllMinGIDs, SC, LO, GO)

TPETRA_ETI_MANGLING_TYPEDEFS()
TPETRA_INSTANTIATE_SLG(THIS_TEST_GROUP)

} // namespace (anonymous)


int
main (int argc, char* argv[])
{
  Tpetra::ScopeGuard tpetra_scope(&argc, &argv);
  const int err_code =
    Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
  return err_code;
}
