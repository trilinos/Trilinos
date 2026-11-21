#include <type_traits>
#include <Compadre_KokkosParser.hpp>

using namespace Compadre;

// KokkosParser() is the form for code without access to the command line; it once failed to compile, being
// ambiguous with a private default constructor
static_assert(std::is_default_constructible_v<KokkosParser>, "KokkosParser() must compile");

TEST (KokkosParser, NoArgumentsJoinsRunningKokkos) {
    // the unit test driver initializes Kokkos, so KokkosParser() joins it, and must not finalize it
    {
        KokkosParser kp;
        EXPECT_TRUE(Kokkos::is_initialized());
    }
    EXPECT_TRUE(Kokkos::is_initialized());
}
