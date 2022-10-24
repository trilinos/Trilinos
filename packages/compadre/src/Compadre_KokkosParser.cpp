#include "Compadre_KokkosParser.hpp"

using namespace Compadre;

// for command line arguments, pass them directly in to Kokkos
KokkosParser::KokkosParser(int narg, char* args[], bool print_status) {
    this->initialize(narg, args, print_status);
}

KokkosParser::KokkosParser(std::vector<std::string> stdvec_args, bool print_status) {
    std::vector<char*> char_args;
    for (const auto& arg : stdvec_args) {
        char_args.push_back((char*)arg.data());
    }
    char_args.push_back(nullptr);
    int narg = (int)stdvec_args.size();
    this->initialize(narg, char_args.data(), print_status);
}

KokkosParser::KokkosParser(bool print_status) {
    std::vector<std::string> stdvec_args;
    stdvec_args.push_back("placeholder");
    std::vector<char*> char_args;
    for (const auto& arg : stdvec_args) {
        char_args.push_back((char*)arg.data());
    }
    char_args.push_back(nullptr);
    this->initialize(1, char_args.data(), print_status);
}

int KokkosParser::initialize(int narg, char* argv[], bool print_status) {
    // return codes:
    // 1  - success
    // 0  - already initialized
    // -1 - failed for some other reason

    // determine if Kokkos is already initialized
    // if it has been, get the parameters needed from it
    bool preinitialized = Kokkos::is_initialized();
    
    // if already initialized, return
    if (preinitialized) {
        if (print_status) printf("Previously initialized.\n");
        return 0;
    } else {
        compadre_assert_release((narg!=0 && argv!=NULL) && "Invalid input to initialize()\n");
        try {
            Kokkos::initialize(narg, argv);
            bool success = Kokkos::is_initialized();
            compadre_assert_release(success && "Kokkos did not initialize successfully.\n");
            _called_initialize = 1;
            if (print_status) {
                this->status();
            }
            return 1;
        } catch (const std::exception& e) {
            std::cout << e.what() << std::endl;
            throw e;
        } catch (...) {
            return -1;
        }
    }
}

int KokkosParser::finalize(bool hard_finalize) {
    if (hard_finalize || _called_initialize==1) {
        try {
            Kokkos::finalize();
            _called_initialize = 0; // reset since we finalized
            return 1;
        } catch (...) {
            return 0;
        }
    } else {
        return 1;
    }
}

void KokkosParser::status() const {
    Kokkos::print_configuration(std::cout, true);
}
