TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  SUBPACKAGES_DIRS_CLASSIFICATIONS_OPTREQS
    Core          core           PT  REQUIRED
    DofMgr        dof-mgr        PT  OPTIONAL
    DiscFE        disc-fe        PT  OPTIONAL
    AdaptersSTK   adapters-stk   PT  OPTIONAL
#    AdaptersIOSS  adapters-ioss  PT  OPTIONAL
    MiniEM        mini-em        PT  OPTIONAL
    ExprEval      expr-eval      EX  OPTIONAL
  )
