INCLUDE(TribitsSubPackageMacros)

TRIBITS_SUBPACKAGE(EpetraExtAdapters)

ADD_SUBDIRECTORY(src)

#TRIBITS_ADD_EXAMPLE_DIRECTORIES(example)

# 2008/10/14: rabartl: Had to comment out exapmle because it has a
# circular dependence on Stratimikos and we can't allow that right
# now!

TRIBITS_ADD_TEST_DIRECTORIES(test)

TRIBITS_SUBPACKAGE_POSTPROCESS()
