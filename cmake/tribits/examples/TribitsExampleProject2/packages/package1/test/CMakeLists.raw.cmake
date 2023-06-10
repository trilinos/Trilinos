add_test(NAME Package1_Prg COMMAND package1-prg)
set_tests_properties(Package1_Prg
  PROPERTIES PASS_REGULAR_EXPRESSION "Package1 Deps: tpl1")
