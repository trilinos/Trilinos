set_default_and_from_env(TribitsExMetaProj_GIT_URL_REPO_BASE  https://github.com/tribits/)

tribits_project_define_extra_repositories(
  TribitsExampleProject   ""  GIT  ${TribitsExMetaProj_GIT_URL_REPO_BASE}TribitsExampleProject.git
    ""   Continuous
  TribitsExampleProjectAddons   ""  GIT  ${TribitsExMetaProj_GIT_URL_REPO_BASE}TribitsExampleProjectAddons.git
    ""   Continuous
  )
