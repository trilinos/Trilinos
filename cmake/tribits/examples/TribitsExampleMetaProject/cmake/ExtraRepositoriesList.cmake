SET_DEFAULT_AND_FROM_ENV(TribitsExMetaProj_GIT_URL_REPO_BASE  https://github.com/tribits/)

TRIBITS_PROJECT_DEFINE_EXTRA_REPOSITORIES(
  TribitsExampleProject   ""  GIT  ${TribitsExMetaProj_GIT_URL_REPO_BASE}TribitsExampleProject.git
    ""   Continuous
  TribitsExampleProjectAddons   ""  GIT  ${TribitsExMetaProj_GIT_URL_REPO_BASE}TribitsExampleProjectAddons.git
    ""   Continuous
  )
