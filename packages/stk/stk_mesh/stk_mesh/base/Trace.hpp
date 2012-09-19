#ifndef stk_mesh_Trace_hpp
#define stk_mesh_Trace_hpp

///////////////////////////////////////////////////////////////////////////////
// Macros/functions for tracing. This file contains the "layer" between
// stk_mesh and the stk-util tracing interface; this layer is necessary because
// one of our goal is to ensure that all tracing code will be compiled-out
// unless tracing is enabled. We also wanted to support the ability to watch
// certain objects/values while tracing and this capability is not supported
// natively in the stk_util interface.
//
// No tracing will occur unless STK_MESH_TRACE_ENABLED is defined. You rebuild
// with tracing with: 'bake -c ; bake <product> cxxflags=-DSTK_MESH_TRACE_ENABLED'
// You'll need to be sure you link with stk_util/use_cases.
//
// If a macro is setting up your main function, be sure to use STKUNIT_WITH_TRACING_MAIN.
// If you're setting up the main function yourself, be sure to set up the
// environment with UseCaseEnvironment.
//
// You'll probably want to add the "--pout" argument so that each processor
// produces it's own trace file.
//
// (Diag|Trace)If will produce diag/trace output if the given PrintMask is
// activated.
//
// (Diag|Trace)IfWatching will produce diag/trace output if the given key is in
// the watch list AND the given PrintMask is activated.
//
// A common pattern for code that wants tracing (put this code somewhere before
// the code you want to trace).
//   stk::mesh::setStream(use_case::dwout());
//   meshlog.setPrintMask(stk::mesh::LOG_ENTITY | stk::mesh::LOG_TRACE | stk::mesh::LOG_TRACE_SUB_CALLS);
//   stk::mesh::watch(stk::mesh::EntityKey(0, 11)); // Node 11
//   stk::diag::Trace::addTraceFunction("stk::mesh::");
//
// Other common items to watch are Parts, Buckets, and Fields
//
// For unit-testing, all trace-related libraries should be header-only so that
// we can define STK_MESH_TRACE_ENABLED and #include the necessary headers and
// be able to trace.
//
// TODO
// * Describe the tracing/diagnostic command-line interface
// * Describe the PrintMask system.
// * Command-line interface to key watching system?
// * Would like to see a "watch-all" operation based on type
// * What if the id's of two different types of objects are the same?
// * Would be nice if Trace("func") supported printing of arg values too
// * Would be nice to have some automated way of checking that the method names
//   are up to date.
///////////////////////////////////////////////////////////////////////////////

#include <stk_mesh/base/DiagWriter.hpp>
#include <stk_mesh/base/EntityKey.hpp>

#include <stk_util/diag/WriterExt.hpp>

#include <string>
#include <typeinfo>
#include <vector>

namespace stk {
namespace mesh {

////////////////////// INTERNAL ////////////////////

class Watch
{
 public:
  virtual ~Watch() {}
  virtual const std::type_info& type() const = 0;
  virtual bool match(const void* item_) const = 0;
  virtual void* item() = 0;
};

std::vector<Watch*>& watch_vector();

template <typename T>
bool internal_is_watching(const T& item)
{
  for (std::vector<Watch*>::const_iterator
       itr = watch_vector().begin(); itr != watch_vector().end(); ++itr) {
    if ((*itr)->type() == typeid(T) &&
        (*itr)->match(&item)) {
      return true;
    }
  }
  return false;
}

template <typename T>
class WatchClass : public Watch
{
 public:
  WatchClass(const T& watch_item) : m_watch_item(watch_item),
                                    m_type_info(&typeid(T)) {
    watch_vector().push_back(this);
  }

  virtual const std::type_info& type() const { return *m_type_info; }

  virtual bool match(const void* item_) const {
    return *(static_cast<const T*>(item_)) == m_watch_item;
  }

  virtual void* item() { return &m_watch_item; }

 private:
  T m_watch_item;
  const std::type_info* m_type_info;
};

/////////////////// API /////////////////////////

template <typename T>
void watch(const T& watch_item)
{
  // leaks, but who cares
  new WatchClass<T>(watch_item);
}

#ifdef STK_MESH_TRACE_ENABLED

inline void setStream(std::ostream& stream)
{
  initDiagWriter(stream);
}

#define Trace_(location) stk::mesh::Trace trace__(location)

#define TraceIf(location, mask) stk::mesh::Trace trace__(location, mask)

#define TraceIfWatching(location, mask, item) \
stk::mesh::Trace trace__(location, mask, stk::mesh::internal_is_watching(item)); \
DiagIfWatching(mask, item, "Watched item is: " << item << stk::diag::dendl)

// Useful if you need two traces in the same scope, dec is used to modify the
// name of the trace object to avoid conflicting.
#define TraceIfWatchingDec(location, mask, item, dec) \
stk::mesh::Trace trace##dec__(location, mask, stk::mesh::internal_is_watching(item)); \
DiagIfWatching(mask, item, "Watched item is: " << item << stk::diag::dendl)

#define DiagIfWatching(mask, item, message)                             \
meshlog.w(stk::mesh::internal_is_watching(item), mask) << message << stk::diag::dendl

#define DiagIf(mask, message)                   \
meshlog.m(mask) << message << stk::diag::dendl

/////////////////////////////////////////////////

#else

inline void setStream(std::ostream& stream) { }

#define Trace_(location)                              ((void) (0))
#define TraceIf(location, mask)                       ((void) (0))
#define TraceIfWatching(location, mask, item)         ((void) (0))
#define TraceIfWatchingDec(location, mask, item, dec) ((void) (0))
#define DiagIf(mask, message)                         ((void) (0))
#define DiagIfWatching(mask, item, message)           ((void) (0))

#endif

} // namespace mesh
} // namespace stk

#endif
