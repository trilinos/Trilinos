
defaults = [
  "-j12",
  "--ctest-timeout=180",
  "--st-extra-builds=",
  "--disable-packages=PyTrilinos,Pliris,STK,TriKota",
  "--skip-case-no-email",
  "--ctest-options=\"-E '(MueLu_ParameterListInterpreterEpetra|MueLu_ParameterListInterpreterTpetra)'\"",
  ]

# Above, I have given up on running --st-extra-builds because they are alway
# so broken.  Just the PT --default-builds are often broken.
