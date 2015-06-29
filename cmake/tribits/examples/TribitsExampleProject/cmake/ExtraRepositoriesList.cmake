SET_DEFAULT_AND_FROM_ENV(GIT_URL_REPO_BASE  git@github.com:TriBITSPub/)

TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES(
  TribitsExampleProjectAddons   ""  GIT  ${GIT_URL_REPO_BASE}TribitsExampleProjectAddons
    ""   Continuous
  )
