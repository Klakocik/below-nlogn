# Publishing to GitHub

The repository is already a git repository with one clean commit on `main`.
Nothing has been pushed anywhere yet.

## Option A — web UI (works right now)

1. Create an empty repository at <https://github.com/new>:
   * owner: `Klakocik`
   * name: `below-nlogn`
   * **do not** initialise with README, .gitignore or license.
2. Then, from the project root:

```console
git remote add origin https://github.com/Klakocik/below-nlogn.git
git push -u origin main
```

## Option B — GitHub CLI

`gh` is not installed on this machine. Install it, authenticate once, then:

```console
winget install --id GitHub.cli -e
gh auth login            # browser or device-code flow
gh repo create below-nlogn --public --source . --remote origin --push
```

## After pushing

* The CI badge in `README.md` already points at
  `Klakocik/below-nlogn/actions/workflows/ci.yml`; the workflow runs on the
  next push.
* If the repository name or owner differs, update the badge URL and the
  `CI` badge link in `README.md`.
* Suggested repository topics: `bigint`, `integer-multiplication`, `ntt`,
  `fft`, `complexity-theory`, `rust`, `research`.
* Suggested description: *"Experimental Rust toolbox for OpenAI's 'Integer
  multiplication below n log n' — plus an exact multiplier that makes no
  impossible promises."*
