// @HEADER
// *****************************************************************************
//               Rapid Optimization Library (ROL) Package
//
// Copyright 2014 NTESS and the ROL contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// Standalone teuchos-lite shim: no-op timers + the few Teuchos/ROL symbols the
// ROL PinT adapter headers use that are otherwise absent without full Trilinos.
// This file is #included (by name) from ROL_PinTVector.hpp, ROL_PinTConstraint.hpp,
// and ROL_PinTHierarchy.hpp, BEFORE their first use of Teuchos::as / Teuchos::null /
// ROL::is_null / TEUCHOS_ASSERT, so defining them here reaches every call site and
// keeps the adapter headers byte-identical to upstream Trilinos ROL.
#ifndef TEUCHOS_TIME_HPP
#define TEUCHOS_TIME_HPP

#include <ostream>
#include <string>
#include <memory>

#include "ROL_Ptr.hpp"        // ROL::Ptr, ROL::makePtr, ROL::nullPtr
#include "Teuchos_Assert.hpp" // TEUCHOS_ASSERT (not transitively reachable otherwise)

namespace Teuchos {

// ---- no-op stacked timer (return value is dereferenced via ->, so it must be live) ----
class StackedTimer {
public:
  void start( const std::string& ) {}
  void stop ( const std::string& ) {}
  void report( std::ostream& )     {} // unused by adapters/mpi/src; kept for parity
};

struct TimeMonitor {
  static ROL::Ptr<StackedTimer> getStackedTimer() {
    static ROL::Ptr<StackedTimer> t = ROL::makePtr<StackedTimer>();
    return t;
  }
};

// ---- Teuchos::as<Out>(in)  (only instantiated as as<int>(size_type)) ----
template<class Out, class In>
inline Out as( const In& x ) { return static_cast<Out>(x); }

// ---- Teuchos::null sentinel, comparable to ROL::Ptr (== std::shared_ptr) ----
struct ENull {};
static const ENull null = ENull{};

template<class T>
inline bool operator==( const std::shared_ptr<T>& p, ENull ) { return p == nullptr; }
template<class T>
inline bool operator!=( const std::shared_ptr<T>& p, ENull ) { return p != nullptr; }
template<class T>
inline bool operator==( ENull, const std::shared_ptr<T>& p ) { return p == nullptr; }
template<class T>
inline bool operator!=( ENull, const std::shared_ptr<T>& p ) { return p != nullptr; }

} // namespace Teuchos

namespace ROL {
// bare is_null(...) in ROL_PinTConstraint.hpp:749 resolves here (call site is inside
// namespace ROL). ROL_Ptr.hpp ships only is_nullPtr; provide the is_null spelling used.
template<class T>
inline bool is_null( const Ptr<T>& x ) { return x == nullPtr; }
} // namespace ROL

#endif // TEUCHOS_TIME_HPP
