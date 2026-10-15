# Creating a Release

This document describes how to cut a release. The `release.yml` GitHub Action
automates sdist/wheel building, PyInstaller binaries, the changelog, the GitHub
Release and PyPI publication.

## Release process

1. Bump the version in `src/pyfortool/__init__.py`:

   ```python
   __version__ = "0.2.21"
   ```

2. Commit the bump and push to `main`:

   ```bash
   git add src/pyfortool/__init__.py
   git commit -m "Version 0.2.21"
   git push origin main
   ```

3. Create and push a tag **exactly matching** the package version (bare version
   numbers, no `v` prefix):

   ```bash
   git tag -a 0.2.21
   git push origin 0.2.21
   ```

   Pushing the tag triggers the `Release` workflow.

## What the workflow does

The `Release` workflow (`on: push: tags`) runs five jobs:

| Job | Runner | Produces |
|-----|--------|----------|
| `tests` | reusable `tests-suite.yml` | Full test suite (pytest matrix, regression, flake8, pylint, docs). Gate for publication/release |
| `build` | ubuntu-latest | sdist `pyfortool-<ver>.tar.gz`, wheel, `CHANGELOG.md` |
| `pyinstaller` (matrix) | ubuntu-22.04 / ubuntu-22.04-arm | `pyfortool-x86_64`, `pyfortool-aarch64` |
| `publish-pypi` | ubuntu-latest | wheel published to PyPI |
| `release` | ubuntu-latest | GitHub Release with all assets attached |

`publish-pypi` and `release` only run if the `tests` job passes (i.e. every job
in the reusable `tests-suite.yml` workflow succeeds); the release is not created
and nothing is published otherwise.

The test checks themselves are defined once in `.github/workflows/tests-suite.yml`
and reused by both `test.yml` (on every push/PR) and `release.yml` (on tag push),
so there is a single source of truth for the test suite.

### Changelog

The changelog is generated from commits between the previous version tag and the
current one. Each listed commit shows the commit title, author and date.

You can preview it locally:

```bash
bash .github/scripts/generate_changelog.sh 0.2.21
```

### PyInstaller binaries

Both binaries follow the same recipe
(`--collect-all pyfxtran --onefile $(which pyfortool) src/pyfortool/scripting.py`):

- `pyfortool-x86_64`: built on `ubuntu-latest` with Python 3.9 (glibc 2.26,
  matches the NUWA build).
- `pyfortool-aarch64`: built on the GitHub-hosted ARM runner
  `ubuntu-24.04-arm` with Python 3.11 (matches the TURPAN build).

## Version / tag consistency

The workflow refuses to proceed if the pushed tag does not equal the version
reported by `pyfortool.__version__`. This prevents publishing a release whose
artifacts report a different version than the release name.

## PyPI

The wheel is published to PyPI using the `PYPI_API_TOKEN` repository
secret set by the maintainers. Publication is gated behind the `publish`
environment in the repository settings.
