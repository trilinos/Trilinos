// Copyright(C) 1999-2017 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
//
//     * Neither the name of NTESS nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "Ioss_DBUsage.h" // for DatabaseUsage
#include <Ioss_IOFactory.h>
#include <Ioss_ParallelUtils.h>
#include <Ioss_Utils.h> // for IOSS_ERROR
#include <cstddef>      // for nullptr
#include <map>          // for _Rb_tree_iterator, etc
#include <ostream>      // for operator<<, basic_ostream, etc
#include <string>       // for char_traits, string, etc
#include <utility>      // for pair
namespace {
#if defined(IOSS_THREADSAFE)
  std::mutex m_;
#endif

  int describe__(Ioss::IOFactoryMap *registry, Ioss::NameList *names)
  {
    int                                count = 0;
    Ioss::IOFactoryMap::const_iterator I;
    for (I = registry->begin(); I != registry->end(); ++I) {
      names->push_back((*I).first);
      ++count;
    }
    return count;
  }
} // namespace

namespace Ioss {
  class DatabaseIO;
  class PropertyManager;
  using IOFactoryValuePair = IOFactoryMap::value_type;
} // namespace Ioss

/** \brief Create an IO database.
 *
 *  This is the public interface method for creating an Ioss::DatabaseIO object.
 *  If Ioss is invoked with more than one process, type == "exodus", db_usage == Ioss::READ_MODEL,
 *  and properties contains the property DECOMPOSITION_METHOD, which is not set to EXTERNAL,
 *  then a parallel decompose-on-the-fly (dof) style Exodus database will be created. The mesh in a
 * single
 *  Exodus file will be decomposed according to DECOMPOSITION_METHOD. Otherwise, if the number of
 *  processes is greater than one, then a file-per-process (fpp) style Exodus database will be
 * created.
 *  In this case, Ioss expects the mesh in p Exodus files, where p is the number of processes in
 *  communicator.
 *
 *  \param[in] type The database file format. Use Ioss::IOFactory::describe to list formats known to
 * Ioss.
 *  \param[in] filename The name of the database file to read from or write to.
 *  \param[in] db_usage Specifies whether the database will be used for input, normal output,
 * restart output, etc.
 *  \param[in] communicator The MPI communicator.
 *  \param[in] properties The property manager associated with the database.
 *  \returns A pointer to the newly-constructed Ioss::DatabaseIO object, or NULL if unsuccessful.
 */
Ioss::DatabaseIO *Ioss::IOFactory::create(const std::string &type, const std::string &filename,
                                          Ioss::DatabaseUsage db_usage, MPI_Comm communicator,
                                          const Ioss::PropertyManager &properties)
{
  IOSS_FUNC_ENTER(m_);
  Ioss::DatabaseIO *db   = nullptr;
  auto              iter = registry()->find(type);
  if (iter == registry()->end()) {
    if (registry()->empty()) {
      std::ostringstream errmsg;
      errmsg << "ERROR: No database types have been registered.\n"
             << "       Was Ioss::Init::Initializer() called?\n\n";
      IOSS_ERROR(errmsg);
    }
    else {
      std::ostringstream errmsg;
      errmsg << "ERROR: The database type '" << type << "' is not supported.\n";
      NameList db_types;
      describe__(registry(), &db_types);
      errmsg << "\nSupported database types:\n\t";
      for (Ioss::NameList::const_iterator IF = db_types.begin(); IF != db_types.end(); ++IF) {
        errmsg << *IF << "  ";
      }
      errmsg << "\n\n";
      IOSS_ERROR(errmsg);
    }
  }
  else {
    auto                my_props(properties);
    Ioss::ParallelUtils pu(communicator);
    pu.add_environment_properties(my_props, false);
    Ioss::IOFactory *factory = (*iter).second;
    db                       = factory->make_IO(filename, db_usage, communicator, my_props);
  }
  return db;
}

/** \brief Get the names of database formats known to IOSS.
 *
 *  \param[out] names The list of known database format names.
 *  \returns The number of known database formats.
 */
int Ioss::IOFactory::describe(NameList *names)
{
  IOSS_FUNC_ENTER(m_);
  return describe__(registry(), names);
}

Ioss::IOFactory::IOFactory(const std::string &type)
{
  registry()->insert(IOFactoryValuePair(type, this));
}

void Ioss::IOFactory::alias(const std::string &base, const std::string &syn)
{
  Ioss::IOFactory *factory = (*registry()->find(base)).second;
  registry()->insert(IOFactoryValuePair(syn, factory));
}

Ioss::IOFactoryMap *Ioss::IOFactory::registry()
{
  static IOFactoryMap registry_;
  return &registry_;
}

/** \brief Empty method.
 */
void Ioss::IOFactory::clean() {}
