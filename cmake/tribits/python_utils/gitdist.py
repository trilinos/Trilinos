#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# Byte array / string / unicode support across Python 2 & 3
#
# Note that the str class in Python 2 is an ASCII string (byte) array and in
# Python 3 it is a Unicode object. For Python 3 code that is backward compatible
# with Python 2, we sometimes need version-specific conversion functions to give
# us the data type we desire. These functions are:
#
#     b(x)    return a byte array of str x, much like b'<string const>' in
#             Python 3
#     s(x)    return a version-specific str object equivalent to x
#
import sys
if sys.version_info < (3,):
  # Python 2
  def b(x): return x
  def s(x): return x
else:
  # Python 3
  import codecs
  def b(x): return codecs.latin_1_encode(x)[0]
  def s(x):
    try:
      return x.decode("utf-8")
    except AttributeError:
      return x

#
# Pieces of the --help documentation
#


distRepoStatusLegend = r"""Legend:
* ID: Repository ID, zero based (order git commands are run)
* Repo Dir: Relative to base repo (base repo shown first with '(Base)')
* Branch: Current branch (or detached HEAD)
* Tracking Branch: Tracking branch (or empty if no tracking branch exists)
* C: Number local commits w.r.t. tracking branch (empty if zero or no TB)
* M: Number of tracked modified (uncommitted) files (empty if zero)
* ?: Number of untracked, non-ignored files (empty if zero)
"""


helpTopics = [
  'overview',
  'repo-selection-and-setup',
  'dist-repo-status',
  'repo-versions',
  'dist-repo-versions-table',
  'aliases', 
  'default-branch',
  'move-to-base-dir',
  'usage-tips',
  'script-dependencies',
  ]

 
def getHelpTopicsStr():
  helpTopicStr = "" 
  for helpTopic in helpTopics:
    helpTopicStr += "* '" + helpTopic + "'\n"
  return helpTopicStr


# Look up help help string given keys from helpTopics array.
helpTopicDefaultIdx = 0;


helpTopicsDict = {}


helpUsageHeader = r"""gitdist [gitdist arguments] <raw-git-command> [git arguments]
       gitdist [gitdist arguments] dist-repo-status
       gitdist [gitdist arguments] dist-repo-versions-table

Run git over a set of git repos in a multi-repository git project (see
--dist-help=overview --help).  This script also includes other tools like
printing a compact repo status table (see --dist-help=dist-repo-status) and
tracking compatible versions through multi-repository SHA1 version files (see
--dist-help=repo-versions).

The options in [gitdist options] are prefixed with '--dist-' and are pulled
out before running 'git <raw-git-command> [git arguments]' in each local git
repo that is processed (see --dist-help=repo-selection-and-setup).
"""


overviewHelp = r"""
OVERVIEW:

Running:

  $ gitdist [gitdist options] <raw-git-command> [git arguments]

will distribute git commands specified by '<raw-git-command> [git arguments]'
across the current base git repo and the set of git repos listed in the file
./.gitdist (or the file ./.gitdist.default, or the argument
--dist-repos=<repo0>,<repo1>,..., see
--dist-help=repo-selection-and-setup).

For example, consider the following base git repo 'BaseRepo' with three other
"extra" git repos cloned under it:

  BaseRepo/
    .git/
    .gitdist
    ExtraRepo1/
      .git/
      ExtraRepo2/
        .git/
    ExtraRepo3/
      .git/

The file .gitdist shown above is created by the user and in this example
should have the contents (note the base repo entry '.'):

  .
  ExtraRepo1
  ExtraRepo1/ExtraRepo2
  ExtraRepo3

For this example, running the command:

  $ cd BaseRepo/
  $ gitdist status

results in the following commands:

  $ git status
  $ cd ExtraRepo1/ ; git status ; ..
  $ cd ExtraRepo1/ExtraRepo2/ ; git status ; ../..
  $ cd ExtraRepo3/ ; git status ; ..

which produces output like:

  *** Base Git Repo: BaseRepo
  On branch master
  Your branch is up-to-date with 'origin/master'.
  nothing to commit, working directory clean
  
  *** Git Repo: ExtraRepo1
  On branch master
  Your branch is up-to-date with 'origin/master'.
  nothing to commit, working directory clean
  
  *** Git Repo: ExtraRepo1/ExtraRepo2
  On branch master
  Your branch is up-to-date with 'origin/master'.
  nothing to commit, working directory clean
  
  *** Git Repo: ExtraRepo3
  On branch master
  Your branch is up-to-date with 'origin/master'.
  nothing to commit, working directory clean

The gitdist tool allows managing a set of git repos like one big integrated
git repo.  For example, after cloning a set of git repos, one can perform
basic operations like for single git repos such as creating a new release
branch and pushing it with:

  $ gitdist checkout master
  $ gitdist pull
  $ gitdist tag -a -m "Start of the 2.3 release" release-2.3-start
  $ gitdist checkout -b release-2.3 release-2.3-start
  $ gitdist push origin release-2.3-start
  $ gitdist push origin -u release 2.3
  $ gitdist checkout master

The above gitdist commands create the same tag 'release-2.3-start' and the
same branch 'release-2.3' in all of the local git repos and pushes these to
the remote 'origin' for each git repo.

For more information about a certain topic, use '--dist-help=<topic-name>
[--help]' for <topic-name>:

"""+getHelpTopicsStr()+r"""
To see full help with all topics, use '--dist-help=all [--help]'.

This script is self-contained and has no dependencies other than standard
python 2.6+ packages so it can be copied to anywhere and used.
"""
helpTopicsDict.update( { 'overview' : overviewHelp } )


repoSelectionAndSetupHelp = r"""
REPO SELECTION AND SETUP:

Before using the gitdist tool, one must first add the gitdist script to one's
default path.  On bash, the simplest way to do this is to source the
gitdist-setup.py script:

  $ source <some-base-dir>/TriBITS/tribits/python_utils/gitdist-setup.sh

This will set an alias to the gitdist script in that same directory by
default, will set up useful alias 'gitdist-status', 'gitdist-mod', and
'gitdist-mod-status', and 'gitdist-repo-versions', and will set up
command-line completion just like for raw git (assuming that
git-completion.bash has been sourced first).  The files 'gitdist' and
'gitdist-setup.sh' can also be copied to another directory (e.g. ~/bin) and
then 'gitdist-setup.sh' can be sourced from there (as a simple "install"):

  $ cp <some-base-dir>/TriBITS/tribits/python_utils/gitdist \
       <some-base-dir>/TriBITS/tribits/python_utils/gitdist-setup.sh \
      ~/bin/
  $ source ~/bin/gitdist-setup.sh
  $ export PATH=$HOME/bin:$PATH

This script can also be set up manually, for example, by copying the gitdist
script to one's ~/bin/ directory:

  $ cp <some-base-dir>/TriBITS/tribits/python_utils/gitdist ~/bin/
  $ chmod a+x ~/bin/gitdist

and then adding $HOME/bin to one's 'PATH' env var with:

  $ export PATH=$HOME/bin:$PATH

(i.e. in one's ~/.bash_profile file).  Then, one will want to set up some
useful shell aliases like 'gitdist-status', 'gitdist-mod', and
'gitdist-mod-status' and 'gitdist-repo-versions' (see --dist-help=aliases).

The set of git repos processed by gitdist is determined by the argument:

  --dist-repos=<repo0>,<repo1>,...

or the files .gitdist or .gitdist.default.  If --dist-repos="", then the list
of repos to process will be read from the file '.gitdist' in the current
directory.  If the file '.gitdist' does not exist, then the list of repos to
process will be read from the file '.gitdist.default' in the current
directory.  The format of this files '.gitdist' and '.gitdist.default' is to
have one repo relative directory per line, for example:

  $ cat .gitdist
  .
  ExtraRepo1
  ExtraRepo1/ExtraRepo2
  ExtraRepo3

where each line is the relative path under the base git repo (i.e. under
'BaseRepo/').  The file .gitdist.default is meant to be committed to the base
git repo (i.e. 'BaseRepo') so that gitdist is ready to use right away after
the base repo and the extra repos are cloned.

If an extra repository directory (i.e. listed in
--dist-repos=<repo0>,<repo1>,..., .gitdist, or .gitdist.default) does
not exist, then it will be ignored by the script.  Therefore, be careful to
manually verify that the script recognizes the repositories that you list.
The best way to do that is to run 'gitdist-status' and see which repos are
listed.

Certain git repos can also be selectively excluded using the option
'--dist-not-repos=<repox>,<repoy>,...'.

Setting up to use gitdist on a specific set of local git repos first requires
cloning and organizing the local git repo. For the example listed here, one
would clone the base repo 'BaseRepo' and the three extra git repos, set up a
.gitdist file, and then add ignores for the extra cloned repos like:

  # A) Clone and organize the git repos
  $ git clone git@some.url:BaseRepo.git
  $ cd BaseRepo/
  $ git clone git@some.url:ExtraRepo1.git
  $ cd ExtraRepo1/
  $ git clone git@some.url:ExtraRepo2.git
  $ cd ..
  $ git clone git@some.url:ExtraRepo3.git

  # B) Create .gitdist
  $ echo .                      > .gitdist
  $ echo ExtraRepo1             >> .gitdist
  $ echo ExtraRepo1/ExtraRepo2  >> .gitdist
  $ echo ExtraRepo3             >> .gitdist

  # C) Add ignores in base repo
  $ echo /ExtraRepo1/ >> .git/info/exclude
  $ echo /ExtraRepo3/ >> .git/info/exclude

  # D) Add ignore in nested extra repo
  $ echo /ExtraRepo2/ >> ExtraRepo1/.git/info/exclude

(Note that one may instead add the above ignores to the version-controlled
files BaseRepo/.gitignore and ExtraRepo1/.gitignore.)

This produces the local repo structure:

  BaseRepo/
    .git/
    .gitdist
    ExtraRepo1/
      .git/
      ExtraRepo2/
        .git/
    ExtraRepo3/
      .git/

After this setup, running:

  $ gitdist <raw-git-command> [git arguments]

in the 'BaseRepo/ 'directory will automatically distribute a given command
across the base repo 'BaseRepo/ and the extra repos ExtraRepo1/,
ExtraRepo1/ExtraRepo2/, and ExtraRepo3/, in that order.

To simplify the setup for the usage of gitdist for a given set of local git
repos, one may choose to instead create the file .gitdist.default in the base
repo (i.e. `BaseRepo/`') and add the ignores for the extra repos to the
.gitignore files and commit the files to the repo(s).  That way, one does not
have to manually do any extra setup for every new set of local clones of the
repos.  But if the file .gitdist is present, then it will override the file
.gitdist.default as described above (which allows customization of what git
repos are processed at any time).
"""
helpTopicsDict.update( { 'repo-selection-and-setup' : repoSelectionAndSetupHelp } )


distRepoStatusHelp = r"""
SUMMARY OF REPO STATUS:

The script gitdist also supports the special command 'dist-repo-status' which
prints a compact table showing the current status of all the repos (see alias
'gitdist-status' in --dist-help=aliases).  For the example set of repos shown
in OVERVIEW (see --dist-help=overview), running:

  $ gitdist dist-repo-status    # alias 'gitdist-status'

outputs a table like:

  ----------------------------------------------------------------------
  | ID | Repo Dir              | Branch | Tracking Branch | C | M  | ? |
  |----|-----------------------|--------|-----------------|---|----|---|
  |  0 | BaseRepo (Base)       | dummy  |                 |   |    |   |
  |  1 | ExtraRepo1            | master | origin/master   | 1 |  2 |   |
  |  2 | ExtraRepo1/ExtraRepo2 | HEAD   |                 |   | 25 | 4 |
  |  3 | ExtraRepo3            | master | origin/master   |   |    |   |
  ----------------------------------------------------------------------

If the option --dist-legend is also passed in, the output will include:

"""+distRepoStatusLegend+\
r"""
One can also show the status of only changed repos with the command:

  $ gitdist dist-repo-status --dist-mod-only  # alias 'gitdist-mod-status'

which produces a table like:

  ----------------------------------------------------------------------
  | ID | Repo Dir              | Branch | Tracking Branch | C | M  | ? |
  |----|-----------------------|--------|-----------------|---|----|---|
  |  1 | ExtraRepo1            | master | origin/master   | 1 |  2 |   |
  |  2 | ExtraRepo1/ExtraRepo2 | HEAD   |                 |   | 25 | 4 |
  ----------------------------------------------------------------------

(see the alias 'gitdist-mod-status' in --dist-help=aliases).

Note that rows for the repos BaseRepo and ExtraRepo2 were left out but the
repo indexes for the remaining repos are preserved.  This allows one to
compactly show the status of the changed local repos even when there are many
local git repos by filtering out rows for repos that have no changes
w.r.t. their tracking branches.  This allows one to get the status on a few
repos with changes out of a large number of local repos (i.e. 10s and even
100s of local git repos).
"""
helpTopicsDict.update( { 'dist-repo-status' : distRepoStatusHelp } )


repoVersionFilesHelp = r"""
REPO VERSION FILES:

The script gitdist also supports the options --dist-version-file=<versionFile>
and --dist-version-file2=<versionFile2> which are used to provide different
SHA1 versions for each local git repo.  Each of these version files is
expected to represent a compatible set of versions of the repos (e.g. in the
same style as .gitmodule files used by the 'git submodule' command).

The format of these repo version files is shown in the following example:

  *** Base Git Repo: BaseRepo
  e102e27 [Mon Sep 23 11:34:59 2013 -0400] <author0@someurl.com>
  First summary message
  *** Git Repo: ExtraRepo1
  b894b9c [Fri Aug 30 09:55:07 2013 -0400] <author1@someurl.com>
  Second summary message
  *** Git Repo: ExtraRepo1/ExtraRepo2
  97cf1ac [Thu Dec 1 23:34:06 2011 -0500] <author2@someurl.com>
  Third summary message
  *** Git Repo: ExtraRepo3
  6facf33 [Fri May 6 15:28:35 2013 -0400] <author3@someurl.com>
  Fourth summary message

Each repository entry can have a summary message or not (i.e. use two or three
lines per repo in the file).  A compatible repo version file can be generated
with this script listing three lines per repo (e.g. as shown above) using (for
example):

  $ gitdist --dist-no-color log -1 --pretty=format:"%h [%ad] <%ae>%n%s" \
    | grep -v "^$" &> RepoVersion.txt

(which is defined as the alias 'gitdist-repo-versions' in the file
'gitdist-setup.sh') or two lines per repo using (for example):

  $ gitdist --dist-no-color log -1 --pretty=format:"%h [%ad] <%ae>" \
    | grep -v "^$" &> RepoVersion.txt

This allows checking out consistent versions of the set git repos, diffing two
consistent versions of the set of git repos, etc.

To checkout an older set of consistent versions of the set of repos
represented by the set of versions given in a file RepoVersion.<date>.txt,
use:

  $ gitdist fetch origin
  $ gitdist --dist-version-file=RepoVersion.<date>.txt checkout _VERSION_

The string '_VERSION_' is replaced with the SHA1 for each of the repos listed
in the file 'RepoVersion.<date>.txt'.  (NOTE: this puts the repos into a
detached head state so one has to know what that means.)

To tag a set of repos using a consistent set of versions, use (for example):

  $ gitdist --dist-version-file=RepoVersion.<date>.txt \
      tag -a -m "<message>" <some_tag> _VERSION_

To create a branch off of a consistent set of versions, use (for example):

  $ gitdist --dist-version-file=RepoVersion.<date>.txt \
      checkout -b some-branch _VERSION_

To diff two sets of versions of the repos, use (for example):

  $ gitdist \
      --dist-version-file=RepoVersion.<new-date>.txt \
      --dist-version-file2=RepoVersion.<old-date>.txt \
      diff _VERSION_ ^_VERSION2_

Here, _VERSION_ is replaced by the SHA1s listed in the file
'RepoVersion.<new-date>.txt' and _VERSION2_ is replaced by the SHA1s listed in
'RepoVersion.<old-date>.txt'.

One can construct any git command taking one or two different repo version
arguments (SHA1s) using this approach (which covers a huge number of different
git operations).

Note that the set of git repos listed in the 'RepoVersion.txt' file must be a
super-set of those processed by this script or an error will occur and the
script will abort (before running any git commands).  If there are additional
repos RepoX, RepoY, etc. not listed in the 'RepVersion'.txt file, then one can
exclude them with:

  $ gitdist --dist-not-repos=RepoX,RepoY,... \
    --dist-version-file=RepoVersion.txt \
    <raw-git-comand> [git arguments]
"""
helpTopicsDict.update( { 'repo-versions' : repoVersionFilesHelp } )


distRepoVersionsTableHelp = r"""
REPO VERSION TABLE:

The script gitdist also supports the special command
'dist-repo-versions-table', which prints a Markdown-formatted table of
repositories and corresponding commit information for easy inclusion in an
issue tracking system.  For instance, running:

  $ gitdist dist-repo-versions-table

outputs a table like:

  | Repository     | SHA1    | Commit Date         | Author                 | Summary                                        |
  |:-------------- |:-------:|:------------------- |:---------------------- |:---------------------------------------------- |
  | MockProjectDir | e2dc488 | 2019-10-23 10:16:07 | user@domain.com        | Merge Pull Request #1234 from user/repo/branch |
  | ExtraRepo1     | f671414 | 2019-10-22 11:18:47 | wile.e.coyote@acme.com | Fixed a Bug                                    |
  | ExtraRepo2     | 50bbf3e | 2019-10-17 16:32:15 | someone@somewhere.com   | Did Some Work                                  |

If the option --dist-short is also passed in, the output will be limited to:

  | Repository     | SHA1    |
  |:-------------- |:-------:|
  | MockProjectDir | e2dc488 |
  | ExtraRepo1     | f671414 |
  | ExtraRepo2     | 50bbf3e |
"""
helpTopicsDict.update( { 'dist-repo-versions-table' : distRepoVersionsTableHelp } )


usefulAliasesHelp =r"""
USEFUL ALIASES:

A few very useful (bash) shell aliases and setup commands to use with gitdist
include:

  $ alias gitdist-status="gitdist dist-repo-status"
  $ alias gitdist-mod="gitdist --dist-mod-only"
  $ alias gitdist-mod-status="gitdist dist-repo-status --dist-mod-only"
  $ alias gitdist-repo-versions="gitdist --dist-no-color log -1 \
    --pretty=format:\"%h [%ad] <%ae>%n%s\" | grep -v \"^$\""

These are added by sourcing the provided file 'gitdist-setup.sh' (which should
be sourced in your ~/.bash_profile file.) which also adds some useful
commandline tab completions.

This avoids lots of extra typing as these gitdist arguments are used a lot.
For example, to see the compact status table of all your local git repos, do:

  $ gitdist-status

To just see a compact status table of only changed repos, do:

  $ gitdist-mod-status

To process only repos that have changes and see commits in these repos
w.r.t. their tracking branches, do (for example):

  $ gitdist-mod log --name-status HEAD ^@{u}

or

  $ gitdist-mod local-stat

(where 'local-stat' is a useful git alias defined in the script
'git-config-alias.sh' which adds these to your ~/.gitconf file).
"""
helpTopicsDict.update( { 'aliases' : usefulAliasesHelp } )

defaultBranchHelp = r"""
DEFAULT BRANCH SPECIFICATION:

When using any git command that accepts a reference (a SHA1, or branch or tag
name), it is possible to use _DEFAULT_BRANCH_ instead.  For instance,

    gitdist checkout _DEFAULT_BRANCH_

will check out the default development branch in each repository being managed
by gitdist.  You can specify the default branch for each repository in your
.gitdist[.default] file.  For instance, if your .gitdist file contains

    . master
    extraRepo1 develop
    extraRepo2 app-devel

then the command above would check out 'master' in the base repo, 'develop' in
extraRepo1, and 'app-devel' in extraRepo2.  This makes it convenient when
working with multiple repositories that have different names for their main
development branches.  For instance, you can do a topic branch workflow like:

    gitdist checkout _DEFAULT_BRANCH_
    gitdist pull
    gitdist checkout -b newFeatureBranch
    <create some commits>
    gitdist fetch
    gitdist merge origin/_DEFAULT_BRANCH_
    <create some commits>
    gitdist checkout _DEFAULT_BRANCH_
    gitdist pull
    gitdist merge newFeatureBranch

and not worry about this 'newFeatureBranch' being off of 'master' in the root
repo, off of 'develop' in extraRepo1, and off of 'app-devel' in extraRepo2.

If no branch name is specified for any given repository in the
.gitdist[.default] file, then 'master' is assumed.
"""
helpTopicsDict.update( { 'default-branch' : defaultBranchHelp } )


moveToBaseDirHelp = r"""
MOVE TO BASE DIRECTORY:

By default, when you run gitdist, it will look in your current working
directory for a .gitdist[.default] file.  If it fails to find one, it will
treat the current directory as the base git repository (as if there was a
.gitdist file in it, having a single line with only "." in it) and then run as
usual.  You have the ability to change this behavior by setting the
GITDIST_MOVE_TO_BASE_DIR environment variable.

To describe the behavior for the differ net options, consider the following set
of nested git repositories and directories:

    BaseRepo/
      .git
      .gitdist
      ...
      ExtraRepo/
        .git
        .gitdist
        ...
        path/
          ...
          to/
            ...
            some/
              ...
              directory/
                ...


The valid settings for GITDIST_MOVE_TO_BASE_DIR include:

  "" (Empty)

    This gives the default behavior where gitdist runs in the current working
    directory.

  IMMEDIATE_BASE

    In this case, gitdist will start moving up the directory tree until it
    finds a .gitdist[.default] file, and then run in the directory where it
    finds it.  In the above example, if you are in
    BaseRepo/ExtraRepo/path/to/some/directory/ when you run gitdist, it will
    move up to ExtraRepo to execute the command you give it from there.

  EXTREME_BASE:

    In this case, gitdist will continue moving up the directory tree until it
    finds the outer-most repository containing a .gitdist[.default] file, and
    then run in that directory.  Given the directory tree above, if you were
    in BaseRepo/ExtraRepo/path/to/some/directory, it will move up to BaseRepo
    to execute the command you give it.

With either of the settings above, when gitdist is finished running, it will
leave you in the same directory you were in when you executed command in the
first place.  Additionally, if no .gitdist[.default] file can be found, gitdist
will execute the command you give it in your current working directory, as if
GITDIST_MOVE_TO_BASE_DIR hadn't been set.
"""
helpTopicsDict.update( { 'move-to-base-dir' : moveToBaseDirHelp } )

usageTipsHelp = r"""
USAGE TIPS:

Since gitdist allows treating a set of git repos as one big git repo, almost
any git workflow that is used for a single git repo can be used for a set of
repos using gitdist.  The main difference is that one will typically need to
create commits individually for each repo.  Also, pulls and pushes are no
longer atomic like is guaranteed for a single git repo.

In general, the mapping between the commands for a single-repo git workflow
using raw git vs. a multi-repo git workflow using gitdist (using the shell
aliases 'gitdist-status', 'gitdist-mod-status', and 'gitdist-mod'; see
--dist-help=aliases) is given by:

  git pull                          =>  gitdist pull
  git checkout -b <branch> [<ref>]  =>  gitdist checkout -b <branch> [<ref>]
  git checkout <branch>             =>  gitdist checkout <branch>
  git tag -a -m "<message>" <tag>   =>  gitdist tag -a -m "<message>" <tag>
  git status                        =>  gitdist-mod status  # status details
                                    =>  gitdist-status      # table for all
                                    =>  gitdist-mod-status  # table for mod.
  git commit                        =>  gitdist-mod commit
  git log HEAD ^@{u}                =>  gitdist-mod log HEAD ^@{u} 
  git push                          =>  gitdist-mod push
  git push [-u] <remote> <branch>   =>  gitdist push [-u] <remote> <branch>
  git push <remote> <tag>           =>  gitdist push <remote> <tag>

NOTE: The usage of 'gitdist-mod' can be replaced with just 'gitdist' in all of
the above commands.  It is just that in these cases gitdist-mod produces more
compact output and avoids do-nothing commands for repos that have no changes
with respect to their tracking branch.  But when it doubt, just use raw
'gitdist' if you are not sure.

A typical development iteration of the centralized workflow using using
multiple git repos looks like the following:

1) Update the local branches from the remote tracking branches:

    $ cd BaseRepo/
    $ gitdist pull

2) Make local modifications for each repo:

    $ emacs <base-files>
    $ cd ExtraRepo1/
    $ emacs <files-in-extra-repo1>
    $ cd ..
    $ cd ExtraRepo1/ExtraRepo2/
    $ emacs <files-in-extra-repo2>
    $ cd ../..
    $ cd ExtraRepo3/
    $ emacs <files-in-extra-repo3>
    $ cd ..

3) Build and test local modifications:

    $ cd BUILD/
    $ make -j16
    $ make test  # hopefully all pass!
    $ cd ..

4) View the modifications before committing:

    $ gitdist-mod-status   # Produces a summary table
    $ gitdist-mod status   # See status details

5) Make commits to each repo:

    $ gitdist-mod commit -a  # Opens editor for each repo in order

  or use the same commit message for all repos:

    $ emacs commitmsg.txt
    $ echo /commitmsg.txt >> .git/info/exclude
    $ gitdist-mod commit -a -F $PWD/commitmsg.txt

  or manually create the commits in each repo separately with raw git:

    $ cd BaseRepo/
    $ git commit -a
    $ cd ExtraRepo1/
    $ git commit -a
    $ cd ..
    $ cd ExtraRepo1/ExtraRepo2/
    $ git commit -a
    $ cd ../..
    $ cd ExtraRepo3/
    $ git commit -a
    $ cd ..

6) Examine the local commits that are about to be pushed:

    $ gitdist-mod-status  # Should be no unmodified or untracked files!
    $ gitdist-mod log --name-status HEAD ^@{u}  # or ...
    $ gitdist-mod local-stat    # alias defined in 'git-config-alias.sh'

7) Rebase and push local commits to remote tracking branch:

    $ gitdist pull --rebase
    $ gitdist-mod push
    $ gitdist-mod-status   # Make sure all the pushes occurred!

Another example workflow is creating a new release branch as shown in the
OVERVIEW section (--dist-help=overview).

Other usage tips:

 - 'gitdist --help' will run gitdist help, not git help.  If you want raw git
   help, then run 'git --help'.

 - Be sure to run 'gitdist-status' to make sure that each repo is on the
   correct local branch and is tracking the correct remote branch.

 - In general, for most workflows, one should use the same local branch name,
   remote repo name, and remote tracking branch name in each local git repo.
   That allows commands like 'gitdist checkout --track <remote>/<branch>' and
   'gitdist checkout <branch>' to work correctly.

 - For many git commands, it is better to process only repos that are changed
   w.r.t. their tracking branch with 'gitdist-mod <raw-git-command> [git
   arguments]'.  For example, to see the status of only changed repos use
   'gitdist-mod status'.  This allows the usage of gitdist to scale well when
   there are even 100s of git repos.

 - As an exception to the last item, a few different types of git commands
   tend to be run on all the git repos like 'gitdist pull', 'gitdist
   checkout', and 'gitdist tag'.

 - If one is not sure whether to run 'gitdist' or 'gitdist-mod', then just
   run 'gitdist' to be safe.
"""
helpTopicsDict.update( { 'usage-tips' : usageTipsHelp } )


scriptDependenciesHelp = r"""
SCRIPT DEPENDENCIES:

The Python script gitdist only depends on the Python 2.6+ standard modules
'sys', 'os', 'subprocess', and 're'. Also, of course, it requires some
compatible version of 'git' in your path (but gitdist works with several
versions of git starting as far back as git 1.6+).
"""
helpTopicsDict.update( { 'script-dependencies' : scriptDependenciesHelp } )


#
# Functions to help Format a table
#


# Shrink a string to a given width by inserting an ellipsis (...) in the
# middle.
def shrinkString(string, width):
  if len(string) > width:
    start = int(width//2) - 1
    stop  = width - start - 3
    return string[:start] + "..." + string[-stop:]
  else:
    return string


# Fill in a field
def getTableField(field, width, just):
  if just == "R":
    return field.rjust(width)
  return field.ljust(width)


# Format an ASCII/UTF-8 table from a set of fields.
#
# The format is of tableData input is:
#
#   [ { "label":"<label0>:, "align":"<align0>, "fields":[<fld00>, ... ]}, 
#     { "label":"<label1>:, "align":"<align1>, "fields":[<fld10>, ... ]},
#     ...
#     ]
#
# The "align" field is either "R" for right, or "L" for left.
#
def createTable(tableData, utf8=False):

  # Table size
  numFields = len(tableData)
  numRows = len(tableData[0]["fields"])

  # a) Get the max field width for each column.
  tableFieldWidth = []
  for fieldDict in tableData:
    label = fieldDict["label"]
    maxFieldWidth = len(label)
    if len(fieldDict["fields"]) != numRows:
      raise Exception("Error: column '"+label+"' numfields = " + \
        str(len(fieldDict["fields"])) + " != numRows = "+str(numRows)+"\n" )
    for field in fieldDict["fields"]:
      fieldWidth = len(field)
      if fieldWidth > maxFieldWidth: maxFieldWidth = fieldWidth
    tableFieldWidth.append(maxFieldWidth)

  # b) Shrink the dist-repo-status table to fit in the terminal if needed.
  shrink = True
  for fieldDict in tableData:
    label = fieldDict["label"]
    if (label != "ID"              and
        label != "Repo Dir"        and
        label != "Branch"          and
        label != "Tracking Branch" and
        label != "C"               and
        label != "M"               and
        label != "?"):
      shrink = False
  if shrink:
    try:
      mockSttySize = os.environ.get("GITDIST_UNIT_TEST_STTY_SIZE")
      if mockSttySize:
        sttySize = mockSttySize
      else:
        sttySize = os.popen("stty size", "r").read()
      rows, columns = sttySize.split()
    except:
      shrink = False
  if shrink:
    terminalWidth = int(columns)
    numDividers = len(tableData) + 1
    numSpaces = 2 * len(tableData)
    fullTableWidth = sum(tableFieldWidth) + numDividers + numSpaces
    if fullTableWidth > terminalWidth:
      widthToShrink = sum(tableFieldWidth[1:4])
      availableWidth = (terminalWidth
                        - tableFieldWidth[0]
                        - sum(tableFieldWidth[4:])
                        - numDividers
                        - numSpaces)
      newWidth = {}
      remainingWidth = availableWidth
      for i in range(1, 3):
        ratio = float(tableFieldWidth[i]) / widthToShrink
        newWidth[i] = int((ratio*availableWidth) // 1)
        remainingWidth = remainingWidth - newWidth[i]
      newWidth[3] = remainingWidth
      for i in range(1, 4):
        if newWidth[i] < len(tableData[i]["label"]):
          shrink = False
          break
      if shrink:
        for i in range(1, 4):
          tableFieldWidth[i] = newWidth[i]
          for j, field in enumerate(tableData[i]["fields"]):
            tableData[i]["fields"][j] = shrinkString(field, tableFieldWidth[i])
        fullTableWidth = terminalWidth

  # c) Write the header of the table (always left-align the column labels).
  table = "┌" if utf8 else "-"
  for index, width in enumerate(tableFieldWidth):
    table += (("─" if utf8 else "-")*(width+2))
    if index != len(tableFieldWidth)-1:
      table += "┬" if utf8 else "-"
    else:
      table += "┐" if utf8 else "-"
  table += "\n"+("│" if utf8 else "|")
  fieldIdx = 0
  for fieldDict in tableData:
    table += " "
    table += getTableField(fieldDict["label"], tableFieldWidth[fieldIdx], "L")
    table += " "+("│" if utf8 else "|")
    fieldIdx += 1
  table += "\n"+("┝" if utf8 else "|")
  for field_i in range(numFields):
    table += (("━" if utf8 else "-")*(tableFieldWidth[field_i]+2))
    if field_i != numFields-1:
      table += "┿" if utf8 else "|"
    else:
      table += "┥" if utf8 else "|"
  table += "\n"

  # d) Write each row of the table
  for row_i in range(numRows):
    table += "│" if utf8 else "|"
    field_i = 0
    for fieldDict in tableData:
      table += " "+getTableField(fieldDict["fields"][row_i],
        tableFieldWidth[field_i], fieldDict["align"] )+" "
      table += "│" if utf8 else "|"
      field_i += 1
    table += "\n"
  table += "└" if utf8 else "-"
  for index, width in enumerate(tableFieldWidth):
    table += (("─" if utf8 else "-")*(width+2))
    if index != len(tableFieldWidth)-1:
      table += "┴" if utf8 else "-"
    else:
      table += "┘" if utf8 else "-"
  table += "\n"

  return table


# Format a Markdown table from a set of fields.
#
# The format of the tableData input is:
#
#   [ { "label":"<label0>:, "align":"<align0>, "fields":[<fld00>, ... ]},
#     { "label":"<label1>:, "align":"<align1>, "fields":[<fld10>, ... ]},
#     ...
#     ]
#
# The "align" field is either "R" for right, "C" for center, or "L" for left.
#
def createMarkdownTable(tableData):

  # Table size
  numFields = len(tableData)
  numRows = len(tableData[0]["fields"])

  # a) Get the max field width for each column.
  tableFieldWidth = []
  for fieldDict in tableData:
    label = fieldDict["label"]
    maxFieldWidth = len(label)
    if len(fieldDict["fields"]) != numRows:
      raise Exception("Error: column '"+label+"' numfields = " + \
        str(len(fieldDict["fields"])) + " != numRows = "+str(numRows)+"\n" )
    for field in fieldDict["fields"]:
      fieldWidth = len(field)
      if fieldWidth > maxFieldWidth:
        maxFieldWidth = fieldWidth
    tableFieldWidth.append(maxFieldWidth)

  # b) Write the header of the table.
  table = "|"
  fieldIdx = 0
  for fieldDict in tableData:
    table += " "
    table += getTableField(fieldDict["label"], tableFieldWidth[fieldIdx],
      fieldDict["align"])
    table += " |"
    fieldIdx += 1
  table += "\n|"
  for i, fieldDict in enumerate(tableData):
    if ((fieldDict["align"] == "L") or (fieldDict["align"] == "C")):
      table += ":"
    else:
      table += " "
    table += "-"*tableFieldWidth[i]
    if ((fieldDict["align"] == "C") or (fieldDict["align"] == "R")):
      table += ":"
    else:
      table += " "
    table += "|"

  # c) Write each row of the table
  for row_i in range(numRows):
    table += "\n|"
    field_i = 0
    for fieldDict in tableData:
      table += " "+getTableField(fieldDict["fields"][row_i],
        tableFieldWidth[field_i], fieldDict["align"] )+" |"
      field_i += 1

  return table


#
# Helper functions for gitdist
#

import sys
import os
import subprocess
import re

from optparse import OptionParser

def addOptionParserChoiceOption(
  optionName,
  optionDest,
  choiceOptions,
  defaultChoiceIndex,
  helpStr,
  optionParser
  ):
  """ Add a general choice option to a optparse.OptionParser object"""
  defaultOptionValue = choiceOptions[defaultChoiceIndex]
  optionParser.add_option(
    optionName,
    dest=optionDest,
    type="choice",
    choices=choiceOptions,
    default=defaultOptionValue,
    help='%s Choices = (\'%s\').  [default = \'%s\']'
    % (helpStr, '\', \''.join(choiceOptions), defaultOptionValue)
    )


def getDistHelpTopicStr(helpTopicVal):
  helpTopicStr = ""
  if helpTopicVal == "":
    return "" # Don't add any text
  elif helpTopicVal == "all":
    for helpTopic in helpTopics:
      helpTopicStr += helpTopicsDict.get(helpTopic)
  else:
    helpTopicHelpStr = helpTopicsDict.get(helpTopicVal, None)
    if helpTopicHelpStr:
      helpTopicStr += helpTopicHelpStr
    else:
      # Invalid help topic so return nonthing and help error handler deal!
      return ""
  return helpTopicStr


def getUsageHelpStr(helpTopicArg):
  usageHelpStr = helpUsageHeader
  if helpTopicArg == "":
    # No help topic option so just use the standard help header
    None
  else:
    helpTopicArgArray = helpTopicArg.split("=")
    if len(helpTopicArgArray) == 1:
      # Option not formatted correctly, set let error handler get it."
      return ""
    (helpTopicArgName, helpTopicVal) = helpTopicArg.split("=")
    usageHelpStr += getDistHelpTopicStr(helpTopicVal)
  return usageHelpStr


def filterWarningsGen(lines): 
  for line in lines:
    if not line.startswith(s('warning')) and not line.startswith(s('error')): yield line


# Filter warning and error lines from output
def filterWarnings(lines): 
  g = filterWarningsGen(lines)
  if g is not None: 
    return list(g)
  return g


# Get output from command
def getCmndOutput(cmnd, rtnCode=False):
  child = subprocess.Popen(cmnd, shell=True, stdout=subprocess.PIPE,
    stderr = subprocess.STDOUT)
  output = child.stdout.read()
  child.wait()
  if rtnCode:
    return (s(output), child.returncode)
  return s(output)


# Run a command and synchronize the output
def runCmnd(options, cmnd):
  if options.debug:
    print("*** Running command: %s" % cmnd)
  if options.noOpt:
    print(cmnd)
  else:
    subprocess.Popen(cmnd, stdout=sys.stdout, stderr=sys.stderr).communicate()
    print("")


# Determine if a command exists:
def commandExists(cmnd):
  whichCmnd = getCmndOutput("which "+cmnd).strip()
  if os.path.exists(whichCmnd):
    return True
  return False


# Get the terminal colors
txtbld=getCmndOutput(r"tput bold")       # Bold
txtblu=getCmndOutput(r"tput setaf 4")    # Blue
txtred=getCmndOutput(r"tput setaf 1")    # Red
txtrst=getCmndOutput(r"tput sgr0")       # Text reset


# Add color to the repo dirs printed out
def addColorToRepoDir(useColor, strIn):
  if useColor:
    return txtbld+txtblu+strIn+txtrst
  return strIn


# Add color to the error messages printed out
def addColorToErrorMsg(useColor, strIn):
  if useColor:
    return txtred+strIn+txtrst
  return strIn


# Get the paths to all the repos gitdist will work on, along with any optional
# default branches.
def parseGitdistFile(gitdistfile):
  reposFullList = []
  defaultBranchDict = {}
  with open(gitdistfile, 'r') as file:
    for line in file:
      line = line.strip()
      if line == "": continue  # ignore blank lines!
      entries = line.split()
      reposFullList.append(entries[0])
      if len(entries) > 1:
        defaultBranchDict[entries[0]] = entries[1]
      else:
        defaultBranchDict[entries[0]] = "master"
  return (reposFullList, defaultBranchDict)


# Get the commandline options
def getCommandlineOps():

  #
  # A) Define the native gitdist command-line arguments
  #

  distHelpArgName = "--dist-help" # Must match --dist-help before --help!
  helpArgName = "--help"
  withGitArgName = "--dist-use-git"
  reposArgName = "--dist-repos"
  notReposArgName = "--dist-not-repos"
  versionFileName = "--dist-version-file"
  versionFile2Name = "--dist-version-file2"
  noColorArgName = "--dist-no-color"
  debugArgName = "--dist-debug"
  noOptName = "--dist-no-opt"
  modifiedOnlyName = "--dist-mod-only"
  legendName = "--dist-legend"
  shortName = "--dist-short"

  nativeArgNames = [ distHelpArgName, helpArgName, withGitArgName, \
    reposArgName, notReposArgName, \
    versionFileName, versionFile2Name, noColorArgName, debugArgName, noOptName, \
    modifiedOnlyName, legendName, shortName ]
  if sys.version_info > (3,):
    utf8Name = "--dist-utf8-output"
    nativeArgNames.append(utf8Name)

  distRepoStatus = "dist-repo-status"
  distRepoVersionTable = "dist-repo-versions-table"
  nativeCmndNames = [ distRepoStatus, distRepoVersionTable ]

  # Select a version of git (see above help documentation)
  defaultGit = "git" # Try system git
  if not commandExists(defaultGit):
    defaultGit = "" # Give up and make the user specify

  #
  # B) Pull the native commandline arguments out of the commandline
  #

  argv = sys.argv[1:]
  nativeArgs = []
  nativeCmnds = []
  otherArgs = []
  helpTopicArg = "" 

  for arg in argv:
    matchedNativeArg = False
    for nativeArgName in nativeArgNames:
      currentArgName = arg[0:len(nativeArgName)]
      if currentArgName == nativeArgName:
        nativeArgs.append(arg)
        matchedNativeArg = True
        if currentArgName == distHelpArgName:
          helpTopicArg = arg
        break
    matchedNativeCmnd = False
    for nativeCmndName in nativeCmndNames:
      if arg == nativeCmndName:
        nativeCmnds.append(nativeCmndName)
        matchedNativeCmnd = True
        break
    if not (matchedNativeArg or matchedNativeCmnd):
      otherArgs.append(arg)

  if len(nativeCmnds) == 0:
    nativeCmnd = None
  elif len(nativeCmnds) == 1:
    nativeCmnd = nativeCmnds[0]
  elif len(nativeCmnds) > 1:
    raise Exception("Error: Can't have more than one dist-xxx command "+\
      " but was passed in "+str(nativeCmnds))

  #
  # C) Set up the commandline parser and parse the native args
  #

  usageHelp = getUsageHelpStr(helpTopicArg)

  clp = OptionParser(usage=usageHelp)

  addOptionParserChoiceOption(
    distHelpArgName, "helpTopic", [""]+helpTopics+["all"], 0,
    "Print a gitdist help topic.  Using" \
    +" --dist-help=all prints all help topics.  If" \
    +" --help is also specified, then the help usage header and" \
    +" command-line 'options' are also printed." ,
    clp )

  clp.add_option(
    withGitArgName, dest="useGit", type="string",
    default=defaultGit,
    help="Path to the git executable to use for each git repo command."
    +"  By default, gitdist will use 'git' in the environment.  If it can't find"
    +" 'git' in the environment, then it will require setting"
    +" --dist-use-git=<path-to-git>.  (Typically only  used in automated"
    +" testing.) (default='"+defaultGit+"')"
    )

  clp.add_option(
    reposArgName, dest="repos", type="string", default="",
    help="Comma-separated list of repo relative paths '<repo0>,<repo1>,...'."
    +" The base repo is specified with '.' and should usually be listed first."
    +" If left empty '', then the list of repos to process is taken from"
    +" the file ./.gitdist (which lists the relative path of each git repo"
    +" separated by newlines).  If the file"
    +" ./.gitdist does not exist, then the repos listed in the file"
    +" ./.gitdist.default are processed.  If the file"
    +" the file ./.gitdist.default is missing, then no extra repos are"
    +" processed and it is assumed that the base repo will be processed."
    +" Also, any git repos listed that don't exist are ignored."
    +" See --dist-help=repo-selection-and-setup."
    +" (default='')"
    )

  clp.add_option(
    notReposArgName, dest="notRepos", type="string", default="",
    help="Comma-separated list of extra repo relative paths" \
    +" '<repoX>,<repoY>,...' to *not* process. (default='')"
    )

  clp.add_option(
    modifiedOnlyName, dest="modifiedOnly", action="store_true",
    help="If set, then only git repos that have changes w.r.t." \
      " their tracking branches will be processed.  That is, only repos" \
      " that have modified or untracked files or where" \
      " 'git diff --name-only ^<tracking-branch>' returns non-empty output" \
      " will be processed (where <tracking-branch> is returned" \
      " from 'rev-parse --abbrev-ref --symbolic-full-name @{u})'." \
      "  If a local repo does not have a tracking branch, then the repo will" \
      " be skipped as well.  Therefore, be careful to first run 'gitdist-status'" \
      " (see --dist-help=dist-repo-status) to see the" \
      " status of each local git repo to know which repos don't have tracking branches.",
    default=False )

  clp.add_option(
    legendName, dest="printLegend", action="store_true",
    help="If set, then a legend will be printed below the repo summary table"\
      " for the special dist-repo-status command.  Only applicable with" \
      " dist-repo-status (see --dist-help=dist-repo-status).",
    default=False )

  if sys.version_info > (3,):
    clp.add_option(
      utf8Name, dest="utf8", action="store_true",
      help="If set, use UTF-8 box drawing characters instead of ASCII ones" \
        " when creating the repo summary table.",
      default=False )

  clp.add_option(
    versionFileName, dest="versionFile", type="string",
    default="",
    help="Path to a file which contains a list of extra repo relative directories"
    +" and git versions (replaces _VERSION_)." \
    +" (See --dist-help=repo-versions.) (default='')"
    )

  clp.add_option(
    versionFile2Name, dest="versionFile2", type="string",
    default="",
    help="Path to a second file contains a list of extra repo relative"
    +" directories and git versions (replaces _VERSION2_)."
    +" (See --dist-help=repo-versions.) (default='')"
    )

  clp.add_option(
    noColorArgName, dest="useColor", action="store_false",
    help="If set, don't use color in the output for gitdist (better for output to a file).",
    default=True )

  clp.add_option(
    debugArgName, dest="debug", action="store_true",
    help="If set, then debugging info is printed.",
    default=False )

  clp.add_option(
    noOptName, dest="noOpt", action="store_true",
    help="If set, then no git commands will be run but instead will just be printed.",
    default=False )

  clp.add_option(
    shortName, dest="short", action="store_true",
    help="If set, then the repo versions table will only include the Repo " \
      "Dir and SHA1 columns; Commit Date, Author, and Summary will be " \
      "omitted.",
    default=False )

  (options, args) = clp.parse_args(nativeArgs)

  debugFromEnv = os.environ.get("GITDIST_DEBUG_OVERRIDE")
  if debugFromEnv:
    options.debug = True

  #
  # D) Print --dist-topic=<topic-name>, check for valid usage
  #

  if options.helpTopic:
    print(getDistHelpTopicStr(options.helpTopic))
    sys.exit(0)

  if not nativeCmnd and len(otherArgs) == 0:
    print(addColorToErrorMsg(options.useColor,
                             "Must specify git command. See 'git --help' for "
                             "options."))
    sys.exit(1)

  if not options.useGit:
    print(addColorToErrorMsg(options.useColor,
                             "Can't find git, please set --dist-use-git"))
    sys.exit(1)

  #
  # E) Change to top-level git directory (in case of nested git repos)
  #

  moveToBaseDir = os.environ.get("GITDIST_MOVE_TO_BASE_DIR")
  if (moveToBaseDir == None) or (moveToBaseDir == ""):
    # Run gitdist in the current directory
    None
  elif moveToBaseDir == "EXTREME_BASE":
    # Run gitdist in the most base dir where .gitdist[.default] exists
    drive, currentPath = os.path.splitdrive(os.getcwd())
    pathList = []
    while 1:
      currentPath, currentDir = os.path.split(currentPath)
      if currentDir != "":
        pathList.append(currentDir)
      else:
        if currentPath != "":
          pathList.append(currentPath)
        break
    pathList.reverse()
    newPath = drive+pathList[0]
    for directory in pathList:
      newPath = os.path.join(newPath, directory)
      if ((os.path.isfile(os.path.join(newPath, ".gitdist"))) or
        (os.path.isfile(os.path.join(newPath, ".gitdist.default")))):
        break
    os.chdir(newPath)
  elif moveToBaseDir == "IMMEDIATE_BASE":
    # Run gitdist in the immediate base dir where .gitdist[.default] exists
    currentPath = os.getcwd()
    foundIt = False
    while 1:
      if ((os.path.isfile(os.path.join(currentPath, ".gitdist"))) or
        (os.path.isfile(os.path.join(currentPath, ".gitdist.default")))):
        foundIt = True
        break
      currentPath, currentDir = os.path.split(currentPath)
      if currentDir == "":
        break
    if foundIt:
      os.chdir(currentPath)
  else:
    print(
      "Error, env var GITDIST_MOVE_TO_BASE_DIR='"+moveToBaseDir+"' is invalid!"
      + "  Valid choices include empty '', IMMEDIATE_BASE, and EXTREME_BASE.")
    sys.exit(1)

  #
  # F) Get the list of extra repos
  #

  if options.repos:
    reposFullList = options.repos.split(",")
    defaultBranchDict = {}
    for repo in reposFullList:
      defaultBranchDict[repo] = "master"
  else:
    if os.path.exists(".gitdist"):
      gitdistfile = ".gitdist"
    elif os.path.exists(".gitdist.default"):
      gitdistfile = ".gitdist.default"
    else:
      gitdistfile = None
    if gitdistfile:
      (reposFullList, defaultBranchDict) = parseGitdistFile(gitdistfile)
    else:
      reposFullList = ["."] # The default is the base repo
      defaultBranchDict = {".": "master"}

  # Get list of not extra repos

  if options.notRepos:
    notReposFullList = options.notRepos.split(",")
  else:
    notReposFullList = []

  #
  # G) Return
  #

  return (options, nativeCmnd, otherArgs, reposFullList, defaultBranchDict,
    notReposFullList)


# Requote commandline arguments into an array
def requoteCmndLineArgsIntoArray(inArgs):
  argsArray = []
  for arg in inArgs:
    splitArg = arg.split("=")
    newArg = None
    if len(splitArg) == 1:
      newArg = arg
    else:
      newArg = splitArg[0]+"="+'='.join(splitArg[1:])
    argsArray.append(newArg)
  return argsArray


# Get a data-structure for a set of repos from a string
def getRepoVersionDictFromRepoVersionFileString(repoVersionFileStr):
  repoVersionFileStrList = repoVersionFileStr.splitlines()
  repoVersionDict = {}
  len_repoVersionFileStrList = len(repoVersionFileStrList)
  i = 0
  while i < len_repoVersionFileStrList:
    repoDirLine = repoVersionFileStrList[i]
    if repoDirLine[0:3] == "***":
      repoDir = repoDirLine.split(":")[1].strip()
      repoVersionLine = repoVersionFileStrList[i+1]
      repoSha1 = repoVersionLine.split(" ")[0].strip()
      repoDirToEnter = ("." if repoDir == baseRepoName else repoDir)
      repoVersionDict.update({repoDirToEnter : repoSha1})
    else:
      break
    nextRepoNoSummary_i = i+2
    if nextRepoNoSummary_i >= len_repoVersionFileStrList:
      break
    if repoVersionFileStrList[nextRepoNoSummary_i][0:3] == "***":
      # Has no summary line
      i = i + 2
    else:
      # Has a summary line
      i = i + 3
  return repoVersionDict


# Get a data-structure for a set of repos from a file
def getRepoVersionDictFromRepoVersionFile(repoVersionFileName):
  if repoVersionFileName:
    repoVersionFileStr = open(repoVersionFileName, 'r').read()
    return getRepoVersionDictFromRepoVersionFileString(repoVersionFileStr)
  else:
    None


def assertAndGetRepoVersionFromDict(repoDirName, repoVersionDict):
  if repoVersionDict:
    
    repoSha1 = repoVersionDict.get(repoDirName, "")
    if not repoSha1:
      print(addColorToErrorMsg(options.useColor,
                               "Repo '" + repoDirName + "' is not in the "
                               + "list of repos " +
                               str(sorted(repoVersionDict.keys())) + " read in from"
                               + " the version file."))
      sys.exit(3)
    return repoSha1
  else:
    return ""


def replaceRepoVersionInCmndLineArg(cmndLineArg, verToken, repoDirName, repoSha1):
  if repoSha1:
    newCmndLineArg = re.sub(verToken, repoSha1, cmndLineArg)
    return newCmndLineArg
  return cmndLineArg


def replaceRepoVersionInCmndLineArgs(cmndLineArgsArray, repoDirName, \
  repoVersionDict, repoVersionDict2 \
  ):
  repoSha1 = assertAndGetRepoVersionFromDict(repoDirName, repoVersionDict)
  repoSha1_2 = assertAndGetRepoVersionFromDict(repoDirName, repoVersionDict2)
  cmndLineArgsArrayRepo = []
  for cmndLineArg in cmndLineArgsArray:
    newCmndLineArg = replaceRepoVersionInCmndLineArg(cmndLineArg, \
      "_VERSION_", repoDirName, repoSha1)
    newCmndLineArg = replaceRepoVersionInCmndLineArg(newCmndLineArg, \
      "_VERSION2_", repoDirName, repoSha1_2)
    cmndLineArgsArrayRepo.append(newCmndLineArg)
  return cmndLineArgsArrayRepo


# Replace _DEFAULT_BRANCH_ in the command line arguments with the appropriate
# default branch name.
def replaceDefaultBranchInCmndLineArgs(cmndLineArgsArray, repoDirName, \
  defaultBranchDict \
  ):
  cmndLineArgsArrayDefaultBranch = []
  for cmndLineArg in cmndLineArgsArray:
    newCmndLineArg = re.sub("_DEFAULT_BRANCH_", \
      defaultBranchDict[repoDirName], cmndLineArg)
    cmndLineArgsArrayDefaultBranch.append(newCmndLineArg)
  return cmndLineArgsArrayDefaultBranch


# Generate the command line arguments
def runRepoCmnd(options, cmndLineArgsArray, repoDirName, baseDir, \
  repoVersionDict, repoVersionDict2, defaultBranchDict \
  ):
  cmndLineArgsArrayRepo = replaceRepoVersionInCmndLineArgs(cmndLineArgsArray, \
    repoDirName, repoVersionDict, repoVersionDict2)
  cmndLineArgsArrayDefaultBranch = replaceDefaultBranchInCmndLineArgs( \
    cmndLineArgsArrayRepo, repoDirName, defaultBranchDict)
  egCmndArray = [ options.useGit ] + cmndLineArgsArrayDefaultBranch
  runCmnd(options, egCmndArray)


# Get the name of the base directory
def getBaseDirNameFromPath(dirPath):
  dirPathArray = dirPath.split("/")
  return dirPathArray[-1]


# Get the name of the base repo to insert into the table
def getBaseRepoTblName(baseRepoName):
  return baseRepoName+" (Base)"


# Determine if the extra repo should be processed or not
def repoExistsAndNotExcluded(options, extraRepo, notReposList):
  if not os.path.isdir(extraRepo): return False
  if extraRepo in notReposList: return False
  return True


# Get the tracking branch for a repo
def getLocalBranch(options, getCmndOutputFunc):
  (resp, rtnCode) = getCmndOutputFunc(
    options.useGit + " rev-parse --abbrev-ref HEAD",
    rtnCode=True )
  if rtnCode == 0:
    filteredLines = filterWarnings(resp.strip().splitlines())
    if filteredLines and len(filteredLines) > 0:
      localBranch = filteredLines[0].strip()
    else:
      localBranch = "<AMBIGUOUS-HEAD>"
    return s(localBranch)
  return ""


# Get the tracking branch for a repo
def getTrackingBranch(options, getCmndOutputFunc):
  (trackingBranch, rtnCode) = getCmndOutputFunc(
    options.useGit + " rev-parse --abbrev-ref --symbolic-full-name @{u}",
    rtnCode=True )
  if rtnCode == 0:
    return s(trackingBranch.strip())
  return ""
  # Above, if the command failed, there is likely no tracking branch.
  # However, this could fail for other reasons so it is a little dangerous to
  # just fail and return "" but I don't know of another way to do this.


# Get number of commits as a str wr.t.t tracking branch
def getNumCommitsWrtTrackingBranch(options, trackingBranch, getCmndOutputFunc):
  if trackingBranch == "":
    return ""
  (summaryLines, rtnCode) = \
    getCmndOutputFunc(options.useGit + " shortlog -s HEAD ^" +
                      trackingBranch, rtnCode=True)
  if rtnCode != 0:
    raise Exception(summaryLines)
  numCommits = 0
  summaryLines = summaryLines.strip()
  if summaryLines:
    for summaryLine in filterWarnings(summaryLines.splitlines()):
      numAuthorCommits = int(summaryLine.strip().split()[0].strip())
      numCommits += numAuthorCommits
  return str(numCommits)
  # NOTE: Above, we would like to use 'git ref-list --count' but that is not
  # supported in older versions of git (at least not in 1.7.0.4).  Using 'git
  # shortlog -s' will return just one line per author so this is not likely to
  # return a lot of data and the cost of the python code to process this
  # should be insignificant compared to the process execution command.


def matchFieldOneOrTwo(findIdx):
  if findIdx == 0 or findIdx == 1:
    return True
  return False


# Get the number of modified
def getNumModifiedAndUntracked(options, getCmndOutputFunc):
  (rawStatusOutput, rtnCode) = getCmndOutputFunc(
    options.useGit + " status --porcelain", rtnCode=True )
  if rtnCode == 0:
    numModified = 0
    numUntracked = 0
    for line in rawStatusOutput.splitlines():
      if matchFieldOneOrTwo(line.find(s("M"))):
        numModified += 1
      elif matchFieldOneOrTwo(line.find(s("A"))):
        numModified += 1
      elif matchFieldOneOrTwo(line.find(s("D"))):
        numModified += 1
      elif matchFieldOneOrTwo(line.find(s("T"))):
        numModified += 1
      elif matchFieldOneOrTwo(line.find(s("U"))):
        numModified += 1
      elif matchFieldOneOrTwo(line.find(s("R"))):
        numModified += 1
      elif line.find(s("??")) == 0:
        numUntracked += 1
    return (str(numModified), str(numUntracked))
  return ("", "")


#
# Get the repo statistics
#

class RepoStatsStruct:

  def __init__(self, branch, trackingBranch, numCommits, numModified, numUntracked):
    self.branch = branch
    self.trackingBranch = trackingBranch
    self.numCommits = numCommits
    self.numModified = numModified
    self.numUntracked = numUntracked

  def __str__(self):
    return "{" \
     "branch='" + self.branch + "'," \
     " trackingBranch='" + self.trackingBranch + "'," \
     " numCommits='" + self.numCommits + "'," \
     " numModified='" + self.numModified + "'," \
     " numUntracked='" + self.numUntracked + "'" \
     "}"

  def numCommitsInt(self):
    if self.numCommits == '': return 0
    return int(self.numCommits)

  def numModifiedInt(self):
    if self.numModified == '': return 0
    return int(self.numModified)

  def numUntrackedInt(self):
    if self.numUntracked == '': return 0
    return int(self.numUntracked)

  def hasLocalChanges(self):
    if self.numCommitsInt() + self.numModifiedInt() + self.numUntrackedInt() > 0:
      return True
    return False


def getRepoStats(options, getCmndOutputFunc=None):
  if not getCmndOutputFunc:
    getCmndOutputFunc = getCmndOutput
  branch         = getLocalBranch(options, getCmndOutputFunc)
  trackingBranch = getTrackingBranch(options, getCmndOutputFunc)
  numCommits     = getNumCommitsWrtTrackingBranch(options,
                                                  trackingBranch,
                                                  getCmndOutputFunc)
  (numModified, numUntracked) = getNumModifiedAndUntracked(options,
                                                           getCmndOutputFunc)
  return RepoStatsStruct(branch,
                         trackingBranch,
                         numCommits,
                         numModified,
                         numUntracked)


class RepoVersionStruct:

  def __init__(self, sha1, commitDate, author, summary):
    self.sha1 = sha1
    self.commitDate = commitDate
    self.author = author
    self.summary = summary

  def __str__(self):
    return "{" \
     "sha1='" + self.sha1 + "'," \
     " commitDate='" + self.commitDate + "'," \
     " author='" + self.author + "'," \
     " summary='" + self.summary + "'" \
     "}"


# Get the SHA1 for the current commit.
def getCommitSha1(options, getCmndOutputFunc):
  (resp, rtnCode) = getCmndOutputFunc(
    options.useGit + " rev-parse --short HEAD",
    rtnCode=True )
  if rtnCode == 0:
    return s(resp.strip())
  return ""

# Get the commit date for the current commit.
def getCommitDate(options, getCmndOutputFunc):
  (resp, rtnCode) = getCmndOutputFunc(
    options.useGit + " log -1 --pretty=format:\"%cd\" --date=format:\"%G-%m-%d %H:%M:%S\"",
    rtnCode=True )
  if rtnCode == 0:
    return s(resp.strip())
  return ""

# Get the author of the current commit.
def getCommitAuthor(options, getCmndOutputFunc):
  (resp, rtnCode) = getCmndOutputFunc(
    options.useGit + " log -1 --pretty=format:\"%ae\"",
    rtnCode=True )
  if rtnCode == 0:
    return s(resp.strip())
  return ""

# Get the first line of the current commit message.
def getCommitSummary(options, getCmndOutputFunc):
  (resp, rtnCode) = getCmndOutputFunc(
    options.useGit + " log -1 --pretty=format:\"%s\"",
    rtnCode=True )
  if rtnCode == 0:
    return s(resp.strip())
  return ""


def getRepoVersions(options, getCmndOutputFunc=None):
  if not getCmndOutputFunc:
    getCmndOutputFunc = getCmndOutput
  sha1 = getCommitSha1(options, getCmndOutputFunc)
  commitDate = getCommitDate(options, getCmndOutputFunc)
  author = getCommitAuthor(options, getCmndOutputFunc)
  summary = getCommitSummary(options, getCmndOutputFunc)
  return RepoVersionStruct(sha1, commitDate, author, summary)


def convertZeroStrToEmpty(strIn):
  if strIn == "0":
    return ""
  return strIn


class RepoStatTable:
  
  def __init__(self):
    self.tableData = [
      { "label" : "ID", "align" : "R", "fields" : [] },
      { "label" : "Repo Dir", "align" : "L", "fields" : [] },
      { "label" : "Branch", "align":"L", "fields" : [] },
      { "label" : "Tracking Branch", "align":"L", "fields" : [] },
      { "label" : "C", "align":"R", "fields" : [] },
      { "label" : "M", "align":"R", "fields" : [] },
      { "label" : "?", "align":"R", "fields" : [] },
      ]

  def insertRepoStat(self, repoDir, repoStat, repoID):
    self.tableData[0]["fields"].append(str(repoID))
    self.tableData[1]["fields"].append(repoDir)
    self.tableData[2]["fields"].append(repoStat.branch)
    self.tableData[3]["fields"].append(repoStat.trackingBranch)
    self.tableData[4]["fields"].append(convertZeroStrToEmpty(repoStat.numCommits))
    self.tableData[5]["fields"].append(convertZeroStrToEmpty(repoStat.numModified))
    self.tableData[6]["fields"].append(convertZeroStrToEmpty(repoStat.numUntracked))

  def getTableData(self):
    return self.tableData

  
class RepoVersionTable:

  def __init__(self):
    self.tableData = [
      { "label" : "Repository", "align" : "L", "fields" : [] },
      { "label" : "SHA1", "align" : "C", "fields" : [] },
      { "label" : "Commit Date", "align" : "L", "fields" : [] },
      { "label" : "Author", "align" : "L", "fields" : [] },
      { "label" : "Summary", "align" : "L", "fields" : [] },
      ]

  def insertRepoVersion(self, repoDir, repoVersion):
    self.tableData[0]["fields"].append(repoDir)
    self.tableData[1]["fields"].append(repoVersion.sha1)
    self.tableData[2]["fields"].append(repoVersion.commitDate)
    self.tableData[3]["fields"].append(repoVersion.author)
    self.tableData[4]["fields"].append(repoVersion.summary)

  def getTableData(self):
    return self.tableData


def getRepoName(repoDir, baseRepoName):
  if repoDir == ".":
    return baseRepoName
  return repoDir
  
#
# Run the script
#

global baseRepoName
baseRepoName = None

if __name__ == '__main__':

  (options, nativeCmnd, otherArgs, reposFullList, defaultBranchDict, \
    notReposList) = getCommandlineOps()

  if nativeCmnd == "dist-repo-status":
    distRepoStatus = True
    if len(otherArgs) > 0:
      print("Error, passing in extra git commands/args ='" + " ".join(otherArgs)
            + "' with special command 'dist-repo-status' is not allowed!")
      sys.exit(1)
  else:
    distRepoStatus = False
  if nativeCmnd == "dist-repo-versions-table":
    distRepoVersionTable = True
    if len(otherArgs) > 0:
      print("Error, passing in extra git commands/args ='" + " ".join(otherArgs)
            + "' with special command 'dist-repo-versions-table' is not allowed!")
      sys.exit(1)
  else:
    distRepoVersionTable = False

  # Get the reference base directory
  baseDir = os.getcwd()

  # Get the name of the base repo
  baseRepoName = getBaseDirNameFromPath(baseDir)

  # Get the repo version files
  repoVersionDict = getRepoVersionDictFromRepoVersionFile(options.versionFile)
  repoVersionDict2 = getRepoVersionDictFromRepoVersionFile(options.versionFile2)

  # Reform the commandline arguments correctly
  cmndLineArgsArray = requoteCmndLineArgsIntoArray(otherArgs)

  if options.debug:
    print("*** Using git: " + str(options.useGit))

  repoStatTable = RepoStatTable()
  repoVersionTable = RepoVersionTable()

  repoID = 0

  for repo in reposFullList:

    # Determine if we should process this repo
    processThisExtraRepo = True
    if not repoExistsAndNotExcluded(options, repo, notReposList):
      processThisExtraRepo = False
    if processThisExtraRepo:
      repoDoesExistsAndNotExcluded = True
      # cd into extrarepo dir
      if options.debug:
        print("\n*** Changing to directory " + repo)
      os.chdir(repo)
      # Get repo stats
      repoStats = None
      if options.modifiedOnly or distRepoStatus:
        repoStats = getRepoStats(options)
      repoVersions = None
      if distRepoVersionTable:
        repoVersions = getRepoVersions(options)
      # See if we should process based on --dist-mod-only
      if options.modifiedOnly and not repoStats.hasLocalChanges():
         processThisExtraRepo = False
    else:
      repoDoesExistsAndNotExcluded = False

    # Process this repo
    if processThisExtraRepo:
      repoName = getRepoName(repo, baseRepoName)
      repoNameInTpl = repoName + (" (Base)" if repo=="." else "") 
      if distRepoStatus:
        repoStatTable.insertRepoStat(repoNameInTpl, repoStats, repoID)
        processThisExtraRepo = False
      elif distRepoVersionTable:
        repoVersionTable.insertRepoVersion(repoName.split("/")[-1],
          repoVersions)
        processThisExtraRepo = False
      else:
        print("")
        print(
          "*** " + ("Base " if repo=="." else "") + "Git Repo: "
          + addColorToRepoDir(options.useColor,repoName) )
        sys.stdout.flush()
        if options.debug:
          print("*** Tracking branch for git repo '" + repoName + "' = '" +
                repoStats.trackingBranch + "'")
        runRepoCmnd(options, cmndLineArgsArray, repo, baseDir, \
          repoVersionDict, repoVersionDict2, defaultBranchDict)
        if options.debug:
          print("*** Changing to directory " + baseDir)

    if repoDoesExistsAndNotExcluded:
      repoID += 1

    os.chdir(baseDir)

  if distRepoStatus:
    if sys.version_info < (3,):
      print(createTable(repoStatTable.getTableData()))
    else:
      print(createTable(repoStatTable.getTableData(), options.utf8))
    if options.printLegend:
      print(distRepoStatusLegend)
    else:
      print("(tip: to see a legend, pass in --dist-legend.)")
  elif distRepoVersionTable:
    if (options.short):
      print(createMarkdownTable(repoVersionTable.getTableData()[0:2]))
    else:
      print(createMarkdownTable(repoVersionTable.getTableData()))
  else:
    print("")

  sys.stdout.flush()
