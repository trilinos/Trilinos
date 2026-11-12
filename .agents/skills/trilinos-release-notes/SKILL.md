---
name: trilinos-release-notes
description: Generates a draft of Trilinos RELEASE_NOTES entries from PRs merged since the previous release and writes them as a patch file. Run from a full (non-partial, non-shallow) clone with the trilinos-release-M-m-branch being documented checked out. Use when asked to write, update, or draft Trilinos release notes or the RELEASE_NOTES file.
---

# Trilinos Release Notes

Drafts `RELEASE_NOTES` entries for a Trilinos release by surveying PRs merged
since the previous release, filtering for ones worth mentioning, confirming
the list with the user, and writing the result as a patch file. This skill
never edits `RELEASE_NOTES` directly and never commits anything.

The skill runs on the release branch being documented
(`trilinos-release-<M>-<m>-branch`), not on `develop`. Release notes are
finalized on the release branch; `develop`'s copy may be stale and its
`Version.cmake` already points at the next development version.

## Notation

These placeholders are used throughout this skill. Replace them with real
values worked out at runtime -- never hardcode a specific Trilinos version,
commit, date, or PR number in this file or in your reasoning about it.

- `<TARGET>` -- the version being documented, `M.m.p`. A *minor* release has
  `p` absent or `0`; a *patch* release has `p > 0`.
- `<BASE>` -- the release the notes are measured from:
  - minor `<TARGET>`: the most recent earlier minor release (`M'.m'` /
    `M'.m'.0`). The notes cover everything merged since `<BASE>`, even if
    patch releases of `<BASE>` happened in between.
  - patch `<TARGET>` `M.m.p`: the previous release on the same branch,
    `M.m.(p-1)` (`M.m.0` when `p = 1`). The notes cover only what is in
    `<TARGET>` and not in `<BASE>`.
- `<BASE_COMMIT>` / `<BASE_DATE>` -- the commit and date that correspond to
  `<BASE>`.
- `<N>` -- a GitHub PR number.
- *Draft section* -- a `RELEASE_NOTES` section whose header date still
  contains a placeholder such as `XX`. *Released section* -- one whose header
  carries an actual date.

Compare versions after normalizing: `M.m` and `M.m.0` are the same version.

## Step 0: Preconditions

Stop and tell the user how to fix it if any of these fail. Do not try to
work around them.

1. **Full clone.** The skill reads old file contents and diffs, which a
   partial or shallow clone does not have locally.
   ```
   git config --get remote.origin.partialclonefilter   # must print nothing
   git rev-parse --is-shallow-repository               # must print false
   ```
   To fix, either make a fresh clone without `--filter`/`--depth`, or convert
   in place (needs network). Unset the filter *before* refetching, or the
   refetch re-applies it:
   ```
   git config --unset remote.origin.partialclonefilter
   git config --unset remote.origin.promisor
   git fetch --refetch origin && git repack -a -d
   ```
   (or `git fetch --unshallow` for a shallow clone).
2. **On a release branch.** `git branch --show-current` must match
   `trilinos-release-<M>-<m>-branch`. If it is `develop` or anything else,
   ask the user to check out the release branch being documented.
3. **Up to date.** Try `git fetch origin`, then compare the branch with its
   upstream (`git status -sb`). If it is behind, warn the user and suggest
   `git pull --ff-only` before continuing; a stale branch silently produces
   wrong `<TARGET>` / `<BASE>` and missing PRs. If the fetch fails (no
   network), don't stop: tell the user the result is only "up to date as of
   the last fetch" and continue. `origin/develop` is used in Step 3, so it
   is subject to the same caveat.
4. **No missing objects in range.** After Step 2 resolves `<BASE_COMMIT>`:
   ```
   git rev-list --objects --missing=print <BASE_COMMIT>..HEAD | grep -c '^?'
   ```
   must print `0`. (Checking only this range is deliberate: an otherwise
   full clone can still lack objects on old, unrelated branches.)

## Step 1: Determine `<TARGET>` and `<BASE>`

1. Read `Trilinos_VERSION` from `Version.cmake` (`V`) and the newest section
   header in `RELEASE_NOTES` (format: `Trilinos Release <version> Release
   Notes <date>`; `H` is its version).
2. `Version.cmake` is authoritative for `<TARGET>`; `RELEASE_NOTES` decides
   whether to extend an existing section or start a new one:

   | Case | Action |
   |---|---|
   | `V == H`, draft section | `<TARGET> = V`; extend that section. |
   | `V > H`, released section | `<TARGET> = V`; start a new section above it. |
   | `V == H`, released section | Ask: add to the already-released notes, or is a new release coming that `Version.cmake` doesn't reflect yet? |
   | `H > V`; `Trilinos_VERSION_STRING` ends in `-dev`; `V`'s `M.m` doesn't match the branch name; or the header doesn't parse | Ask the user; don't guess. |

   On release branches `Version.cmake` is normally bumped before the notes
   are written, so `V > H` with a released section is the usual starting
   point for a new release.
3. `<BASE>`:
   - minor `<TARGET>`: the newest *released* minor-release header in
     `RELEASE_NOTES` below `<TARGET>`.
   - patch `<TARGET>` `M.m.p`: `M.m.(p-1)`, which should be the newest
     released header with the same `M.m`.
4. Let the user override `<TARGET>`, `<BASE>`, or the release date for the
   header at any point. If they do, skip the corresponding detection step.

## Step 2: Resolve `<BASE_COMMIT>` and `<BASE_DATE>`

Release tags may not be present in the clone; don't depend on them.

- **Minor `<TARGET>`:** use the fork point from the previous release branch:
  ```
  git merge-base HEAD origin/trilinos-release-<BASE M>-<BASE m>-branch
  ```
- **Patch `<TARGET>`:** try, in order:
  1. A release tag for `<BASE>`, if one exists (`git tag -l '*<BASE>*'`).
  2. The first parent of the first-parent commit that changed
     `Trilinos_VERSION` in `Version.cmake` to `<TARGET>`:
     ```
     git log --first-parent --format='%H %ci %s' -- Version.cmake
     ```
     Walk back from HEAD to the oldest commit whose `Version.cmake` has
     `<TARGET>` while its first parent's does not; `<BASE_COMMIT>` is that
     first parent. Sanity-check that `RELEASE_NOTES` at `<BASE_COMMIT>` has
     `<BASE>` as its newest released header.

     This is an inference: fixes can land before the version bump, so it
     may be too late. Always let the user correct it.
- If neither works, ask the user to supply `<BASE_COMMIT>` or `<BASE_DATE>`
  directly.

`<BASE_DATE>` is the committer date of `<BASE_COMMIT>`. Use it only for the
optional GitHub query in Step 4 -- never decide release membership by date.
Commit dates don't track when a change reached a branch; use ancestry.

Report the resolved `<TARGET>`, `<BASE>`, `<BASE_COMMIT>`, and `<BASE_DATE>`
to the user before continuing, so they can catch a wrong guess early. Then
run the missing-objects check from Step 0.

## Step 3: Collect candidate PRs from git

Everything in this step uses only local git history. It needs a full clone
(Step 0) but no network access.

Do **not** use `--first-parent`. Release branches receive most changes
through `master_merge_*` PRs (`develop` -> `master` -> release branch), so
first-parent history shows only that plumbing and misses nearly every real
PR.

### Minor `<TARGET>`

- **Merge commits** (created by GitHub's "Create a merge commit" button):
  ```
  git log --merges <BASE_COMMIT>..HEAD --format='%H%x09%s%x09%b'
  ```
  Match subjects like `Merge pull request #<N> from <user>/<branch>`. The PR
  title is the first line of the commit body.

  Keep only merges that lie on the first-parent line of `origin/develop` or
  of HEAD (`git rev-list --first-parent origin/develop HEAD ^<BASE_COMMIT>`).
  Other "Merge pull request" commits are merges inside a contributor's fork
  (typically with small, unrelated PR numbers) and are not Trilinos PRs.

- **Squash merges** (no merge commit, PR number folded into the subject):
  ```
  git log --no-merges <BASE_COMMIT>..HEAD --format='%H%x09%s'
  ```
  Match a trailing `(#<N>)` in the subject. Both merge styles are present in
  Trilinos history, so check for both, then de-duplicate by `<N>`.

- **Drop branch plumbing**: PRs whose source branch is `master_merge_*`,
  PRs that change only `Version.cmake` and/or `RELEASE_NOTES` (version bumps,
  earlier release-notes PRs), and merges of `master` or release branches.
  Keep the PRs they carry -- those are already in the list because Step 3
  doesn't use `--first-parent` to walk history.

- **Expand cherry-pick carrier PRs.** After dropping plumbing (which also
  has commits that aren't on `develop`), a PR merged directly into the
  release branch (often titled like "Cherry picks to M.m release branch")
  may bundle fixes copied from `develop`. Detect it by checking whether its
  commits are absent from `develop`:
  ```
  git rev-list --no-merges <merge>^1..<merge>^2 ^origin/develop
  ```
  If that prints anything, drop the carrier PR itself and map each listed
  commit to its `develop` PR as described in "Mapping a cherry-pick to its
  PR" below.
  Commits that don't map stay attributed to the carrier PR and are flagged
  to the user.

### Patch `<TARGET>`

Patch fixes typically arrive as cherry-picks with no PR merge of their own
on the release branch, so look at individual commits:

1. List `git log --no-merges <BASE_COMMIT>..HEAD --format='%H%x09%s'`.
2. Skip version-bump and release-notes commits (commits touching only
   `Version.cmake` or only `RELEASE_NOTES`).
3. For each remaining commit, find its PR, in order:
   - a trailing `(#<N>)` in the subject;
   - a `Merge pull request #<N>` on the release branch that contains the
     commit, unless that PR is itself plumbing or a cherry-pick carrier;
   - its `develop` PR, as described in "Mapping a cherry-pick to its PR".
4. List any commit that can't be mapped to a PR for the user; don't drop it
   silently.

### Mapping a cherry-pick to its PR

Used for both minor and patch targets.

1. Build the patch-id table for `develop` once, over the whole range since
   `<BASE_COMMIT>` (don't narrow it by date -- commit dates are unreliable):
   ```
   git log -p --no-merges <BASE_COMMIT>..origin/develop | git patch-id --stable
   ```
   Each output line is `<patch-id> <develop-sha>`.
2. For each cherry-picked commit, compute `git show <sha> | git patch-id
   --stable` and look it up in the table. If there's no match (the pick was
   edited, e.g. to resolve conflicts), fall back to a `develop` commit with
   an identical subject in the same range.
3. Find the PR that brought the `develop` commit in: the *first* merge on
   the ancestry path from that commit to `origin/develop` that is also on
   `develop`'s first-parent line:
   ```
   git rev-list --first-parent origin/develop ^<BASE_COMMIT>      # FP set
   git rev-list --merges --ancestry-path --reverse <develop-sha>..origin/develop
   ```
   Walk the second list in order and take the first commit that is in the
   FP set; its subject is `Merge pull request #<N> ...`. Don't combine
   `--first-parent` with `--ancestry-path` in one command: when a PR's
   branch reached `develop` through a later merge, that combination returns
   nothing.

### Both cases

- **Package attribution**: prefer the `Package:` prefix in the PR title (the
  PR template asks contributors to title PRs this way). If that's missing or
  ambiguous, fall back to the top-level `packages/<pkg>` directories touched,
  via `git diff --stat <sha>^1 <sha>` (or `<sha>^ <sha>` for a squash commit).

- **Already-documented PRs**: drop any candidate whose `/pull/<N>` URL
  already appears in `RELEASE_NOTES`, anywhere between the `<BASE>` section
  and the top of the file. For a minor `<TARGET>`, this covers PRs already
  recorded in a patch release that shipped between `<BASE>` and `<TARGET>`,
  without having to name that release.

## Step 4: Optional GitHub cross-check (read-only)

The agent may have been granted read-only access to
`github.com/trilinos/Trilinos`. If so, use it to catch PRs that are hard to
judge from the commit message alone -- never as the sole source of candidates,
since this step is optional and may be unavailable.

1. Probe access once, e.g. `gh api repos/trilinos/Trilinos` or a plain fetch
   of `api.github.com`. If it fails, say so in one line and continue using
   git-only candidates. Do not treat this as an error.
2. If access works:
   - List PRs labeled `Release Note` merged on or after `<BASE_DATE>`:
     ```
     gh pr list -R trilinos/Trilinos --state merged --label "Release Note" \
       --search "merged:>=<BASE_DATE>" --json number,title,body,labels
     ```
     (or the equivalent `api.github.com/search/issues` query.)
   - Treat every PR with this label as must-include, regardless of the
     judgment call in Step 5.
   - For candidates from Step 3 whose inclusion is unclear from the subject
     alone, fetch the PR body for more context.
3. Use read-only operations only. Never comment, label, or otherwise write
   to the upstream repository.

## Step 5: Decide what's worth mentioning

**Include:**
- New features, options, or capabilities.
- API or user-visible behavior changes.
- Deprecations and removals.
- Notable performance improvements.
- New backend or GPU support.
- Version or snapshot bumps for bundled/vendored packages (e.g. Kokkos,
  KokkosKernels, STK, Krino, Compadre, SEACAS) -- these are usually easy to
  spot from a title like "Inclusion of version X" or "update snapshot to vX".
- Anything carrying the `Release Note` label (see Step 4).

**Exclude:**
- Warning fixes, compiler-flag tweaks with no behavior change.
- CI, autotester (AT2), or GitHub workflow-only changes.
- Test-only changes, doc typos, comment fixes.
- Dependabot and other dependency-bot PRs.
- Small, purely internal build fixes.
- Revert commits and the PRs they revert (drop both, unless the net effect
  after the revert is itself worth noting).

When the title and body aren't enough to decide, look at the diff (`git show
--stat <sha>` or a scoped diff) before excluding or including.

## Step 6: Confirm with the user before writing anything

Present the findings and wait for confirmation, edits, or trims:

- A table grouped by package: PR number, title, a proposed one-line summary,
  and why it was included (`label` or `judgment`).
- A collapsed/short list of excluded PRs with a one-phrase reason each, so
  the user can rescue anything mis-filtered.

Do not proceed to Step 7 until the user has signed off on the list.

## Step 7: Match the house style

Look at the newest *complete* section already in `RELEASE_NOTES` (not the
one you're drafting) and mirror its formatting exactly -- this keeps the
draft in sync with however the style has evolved, rather than copying a
style frozen at the time this skill was written.

General shape, using placeholders only:

```
###############################################################################
#                                                                             #
# Trilinos Release <TARGET> Release Notes                       <Month DD, YYYY> #
#                                                                             #
###############################################################################


<Package>:

  - <One- to few-sentence summary, wrapped at roughly 80 columns, with
    continuation lines indented four spaces.>
    https://github.com/trilinos/Trilinos/pull/<N>
    https://github.com/trilinos/Trilinos/pull/<N>


<NextPackage>:
  ...
```

Rules:
- The `#`-box header line must render as a fixed total width matching the
  surrounding boxes in the file (count columns in the reference section you
  copied, don't assume the number above is exact) -- right-align the date so
  the right-hand `#` lines up with the boxes above and below it.
- Packages are ordered alphabetically. `Kokkos` and `KokkosKernels` are
  combined into one `Kokkos, KokkosKernels:` heading when both are touched
  by the same set of PRs.
- Bullets start with `  - `; one bullet can cite multiple PR URLs when the
  PRs are part of the same piece of work.
- Two blank lines separate packages; follow the reference section's spacing
  before the next release's box.
- A snapshot/version bump may also link to the upstream changelog, when one
  exists, as an extra line under the PR link(s).
- A new patch-release section opens with the line `This is a patch release
  with no new features.` (followed by blank lines matching the reference),
  unless the most recent patch section in the file does something different.

## Step 8: Write the patch, not the file

1. Copy `RELEASE_NOTES` to a scratch file and edit the copy only.
   - If a `<TARGET>` draft section already exists (Step 1), add the new
     bullets into it, preserving existing content and ordering.
   - Otherwise, insert a brand-new section at the top of the file.
2. Diff the original against the scratch copy and rewrite the paths so the
   patch applies at the repo root:
   ```
   git diff --no-index RELEASE_NOTES <scratch-copy>
   ```
   then rewrite the `---`/`+++` lines to `a/RELEASE_NOTES` / `b/RELEASE_NOTES`.
3. Save the result as `release-notes-<TARGET>.patch` in the repository root
   (sibling to `RELEASE_NOTES`).
4. Sanity-check it with `git apply --check release-notes-<TARGET>.patch`.
5. Never modify the real `RELEASE_NOTES`, and never `git add`/`commit`/push
   anything. Tell the user the patch is ready and how to apply it:
   `git apply release-notes-<TARGET>.patch`.

## Pitfalls

- Backports/cherry-picks can appear twice (once on a release branch, once on
  develop) -- Step 3's "already-documented" check only catches the ones
  already written into `RELEASE_NOTES`; watch for duplicate PR numbers
  within your own candidate list too.
- A PR touching several packages may deserve a bullet under more than one
  package heading if the change is significant in each; otherwise list it
  once under the primary package.
- `master_merge_*` PRs (and merges of `master` or other release branches)
  are plumbing, not changes -- skip them, but do scan the PRs they carry.
  This is why Step 3 must not use `--first-parent`.
- Cherry-picked fixes on a release branch usually have no PR merge and no
  `(#<N>)` of their own; a single "release branch" PR may carry several of
  them, in minor as well as patch releases. Map each commit back to its
  `develop` PR (Step 3, "Mapping a cherry-pick to its PR") rather than
  citing the carrier PR.
- Decide what's in a release by ancestry (`git merge-base --is-ancestor`),
  never by commit date -- a commit authored before a release can still have
  landed on the branch after it.
- A stale local branch silently gives the wrong `<TARGET>`, `<BASE>`, and
  PR list. Check Step 0's up-to-date condition before trusting results.
- Large snapshot-update PRs (e.g. a Kokkos bump) typically touch huge diffs;
  summarize as "Inclusion of version X", don't enumerate the diff.
- Don't invent detail that isn't supported by the commit message, diff, or
  PR body/label -- when unsure, ask the user rather than guessing.
