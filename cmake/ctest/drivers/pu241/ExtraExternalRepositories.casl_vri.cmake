
#
# List of extra external repositories for CASL.
#
# See documentation in Trilinos preCopyrightTrilinos/ExtraExternalRepositories.cmake
#

SET( Trilinos_EXTRAREPOS_DIR_REPOTYPE_REPOURL_PACKSTAT_CATEGORY
  #Dakota  packages/TriKota/Dakota SVN
  #   https://software.sandia.gov/svn/public/dakota/public/trunk  NOPACKAGES  Nightly  
  StarCCMExt      ""                              GIT  casl-dev.ornl.gov:/git-root/StarCCMClient  ""  Continuous
  DeCARTExt       ""                              GIT  casl-dev.ornl.gov:/git-root/casl_decart    ""  Continuous
  NeutronicsExt   "NeutronicsExt"                 GIT  casl-dev.ornl.gov:/git-root/denovoExt      ""  Continuous
  Denovo          "NeutronicsExt/denovo"          GIT  angmar.ornl.gov/data/git/denovo.git        NOPACKAGES  Continuous
  Nemesis         "NeutronicsExt/denovo/src/nemesis"  GIT  angmar.ornl.gov/data/git/nemesis.git   NOPACKAGES  Continuous
  CASLBOA         ""                              GIT  casl-dev.ornl.gov:/git-root/casl_boa       ""  Continuous
  CASLRAVE        ""                              GIT  casl-dev.ornl.gov:/git-root/casl_rave      ""  Continuous
  LIMEExt         ""                              GIT  software.sandia.gov:/space/git/LIMEExt     ""  Continuous
  PSSDriversExt   ""                              GIT  casl-dev.ornl.gov:/git-root/casl_vripss    ""  Continuous
  )
