// Teko_KrylovSurrogateInit.cpp
//
// Registers Teko::KrylovSurrogate::adaptiveLoop as the global adaptive hook
// with Belos::AdaptiveHook.  The registration happens via a file-scope
// static object whose constructor runs before main().
//
// Because Teko is a shared library (libteko.so) all translation units are
// linked, so this constructor always executes when teko_ext.so loads Teko.

#include "Teko_KrylovSurrogate.hpp"

namespace {

struct TekoAdaptiveRegistration {
    TekoAdaptiveRegistration() {
        Belos::AdaptiveHook::registerHook(&Teko::KrylovSurrogate::adaptiveLoop);
    }
};

// One static instance — constructor fires at dynamic-library load time.
static TekoAdaptiveRegistration g_tekoAdaptiveRegistration;

} // anonymous namespace
