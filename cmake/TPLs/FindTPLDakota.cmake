INCLUDE(TPLDeclareLibraries)

TPL_DECLARE_LIBRARIES( Dakota
  REQUIRED_HEADERS sandia_rules.H
  REQUIRED_LIBS_NAMES dakota pecos fftw3 lhs evidence surfpack conmin ddace dot fsudace  jega nlpql npsol opt psuade newmat  ncsuopt gsl quadrature coliny colin pebbl utilib 3po nappspack appspack conveyor shared cdd amplsolver lhs
  )
