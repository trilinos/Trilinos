#ifndef INTREPID2_UNIT_TEST_DISCRETIZATION_MACROS_HPP
#define INTREPID2_UNIT_TEST_DISCRETIZATION_MACROS_HPP

#define INTREPID2_TEST_ERROR_EXPECTED( S )                                                                       \
    try {                                                                                                        \
      ++nthrow;                                                                                                  \
      S ;                                                                                                        \
    }                                                                                                            \
    catch (std::exception &err) {                                                                                \
      ++ncatch;                                                                                                  \
      *outStream << "Expected Error ----------------------------------------------------------------\n";         \
      *outStream << err.what() << '\n';                                                                          \
      *outStream << "-------------------------------------------------------------------------------" << "\n\n"; \
    }

#define ConstructWithLabel(obj, ...) obj(#obj, __VA_ARGS__)

#endif // INTREPID2_UNIT_TEST_DISCRETIZATION_MACROS_HPP
