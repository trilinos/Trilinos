# Trilinos Git Push Logger

GitHub does not provid informatioin on pushes to a git repo (i.e. you can't se
the reflog for a GitHub repo).  Therefore, one must implement one's own git
push logger.  That is what the scripts in this directory do.

This directory contains some simple scripts to log pushed commits to the main
Trilinos git repo.  To set this up, simply do the following:

  $ cd <some-base-dir>/
  $ git clone git@github.com:trilinos/Trilinos.git
  $ cd Trilinos/
  $ git checkout --track origin/develop
  $ git branch -d master

Then set up a cron entry (using `crontab -e`) that runs the script
loop_log_pushed_commits.sh as:

  1  0  *  *  *  cd <some-base-dir>/ \
    && cp loop_log_pushed_commits.log loop_log_pushed_commits.last.log \
    && ./loop_log_pushed_commits.sh &> loop_log_pushed_commits.log

(but take out the line continuation lines, cron will not allow them).

This will create and continuously append the top commit log message pulled
from the Trilinos 'develop' branch off of GitHub every minute and log that in
the file:

  <some-base-dir>/TrilinosPushLog.txt

Therefore, if pushes to Trilinos don't overlap within one minute, then the the
file TrilinosPushLog.txt will contain a record of all of the pushes to
Trilinos.  (But if multiple pushes do occur in the same minute, then only the
top commit for the most recent push will get logged.  But multiple pushes to
Trilinos in the same day are very rare.)

Also note that this file TrilinosPushLog.txt will only be seen on the local
system where it is created.  For a wider number of people to see this, it will
need to be published in some way.

NOTE: With this process, the file TrilinosPushLog.txt will get appended to
forever and therefore will grow forever.  But it will not grow in size very
fast since these scripts only log the commit message and not the list of files
that got updated.  But one may want to archive the TrilinosPushLog.txt file on
a yearly basis if it starts to get too big.
