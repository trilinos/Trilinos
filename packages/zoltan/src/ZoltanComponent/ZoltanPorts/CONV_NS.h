#ifndef CONV_NS_HSeen
#define CONV_NS_HSeen

/*
  10/02/02 Jaideep Ray, SNL CA

  After the framework got babelized, cca-spec got fragmented into cca-spec-classic
  where ::gov::cca::X became ::classic::gov::cca::X. Ben says this may change again,
  so using namespace classic may not be a good solution. hence I am macro-ing all
  the ::gov::cca::X out into CONV_NS(X).
*/

#ifdef OLD_CCA

#define CONV_NS(X) ::gov::cca::X

#else

#define CONV_NS(X) ::classic::gov::cca::X

#endif


#endif
