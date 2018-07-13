# Trilinos Git Push Logger

GitHub does not provide information on pushes to a git repo (i.e. one can't
see the reflog for a GitHub repo).  Therefore, one must implement one's own
git push logger.  That is what the scripts in this directory do.

This directory contains some simple scripts to log pushed commits to the main
Trilinos git repo on any arbitrary tracking branch.  To set this up to log
commits pushed to the 'develop' branch, for example, simply do the following:

```
$ cd <some-base-dir>/
$ touch TrilinosPushLog.txt
$ touch loop_log_pushed_commits.log
$ git clone git@github.com:trilinos/Trilinos.git
$ cd Trilinos/
$ git checkout --track origin/develop
$ git branch -d master
```

Then set up a cron entry (using `crontab -e`) that runs the script
loop_log_pushed_commits.sh as:

```
  1  0  *  *  *  cd <some-base-dir>/ \
    && cp loop_log_pushed_commits.log loop_log_pushed_commits.last.log \
    && ./Trilinos/commonTools/git/log_pushes/loop_log_pushed_commits.sh \
    &> loop_log_pushed_commits.log
```

But take out the line continuation lines `\` and just list this on one line.
(cron will not allow continuation lines.  They are just added above to improve
readability.)

This will create and continuously prepend the top commit log message pulled
from the Trilinos 'develop' branch off of GitHub every minute and log that in
the file:

```
  <some-base-dir>/TrilinosPushLog.txt
```

Therefore, if pushes to Trilinos don't overlap within the same minute since
the last logging period, then the file TrilinosPushLog.txt will contain a
record of all of the pushes to Trilinos.  (But if multiple pushes do occur in
the last minute, then only the top commit for the most recent push will get
logged.  But multiple pushes to Trilinos in the same minute are very rare.)

Also note that this file TrilinosPushLog.txt will only be seen on the local
system where it is created.  For a wider number of people to see this, it will
need to be published in some way.

NOTE: With this process, the file TrilinosPushLog.txt will get appended to
forever and therefore will grow forever.  But it will not grow in size very
fast since these scripts only log the commit message and not the list of files
that got updated.  It will log the full commit log for the top commit but will
only log the oneline commit summary for all of the other commtis pulled in the
last minutes.  One may want to archive the TrilinosPushLog.txt file and then
start over on a yearly basis if it starts to get too big.
