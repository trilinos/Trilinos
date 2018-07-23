The file ProjectCiFileChangeLogic.py in this directory is only a symlink to
the real file under:

```
  commonTools/framework/ProjectCiFileChangeLogic.py
```

That is because that is under the base directory for the
`TrilinosFrameworkTests` TriBITS package.  We want changes to the file
`ProjectCiFileChangeLogic.py` to trigger the enable of the
TrilinosFrameworkTests package so as to run the unit tests for that module.
But the file `ProjectCiFileChangeLogic.py` also needs to exist in this this
directory so that TriBITS tools will find it and use it automatically.  Using
a symbolici link satisifes both of these requirements.
