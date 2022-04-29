set_default_and_from_env(GIT_URL_REPO_BASE  git@github.com:TriBITSPub/)

tribits_project_define_extra_repositories(
  TribitsExampleProjectAddons   ""  GIT  ${GIT_URL_REPO_BASE}TribitsExampleProjectAddons
    ""   Continuous
  )
