INCLUDE(SubPackageMacros)

SUBPACKAGE(EpetraExtAdapters)

ADD_SUBDIRECTORY(src)

#PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

# 2008/10/14: rabartl: Had to comment out exapmle because it has a
# circular dependence on Stratimikos and we can't allow that right
# now!

PACKAGE_ADD_TEST_DIRECTORIES(test)

SUBPACKAGE_POSTPROCESS()
