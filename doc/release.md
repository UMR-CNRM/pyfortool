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
| `pyinstaller` (matrix) | ubuntu-latest / ubuntu-26.04-arm (Docker containers) | `pyfortool-x86_64`, `pyfortool-aarch64` |
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
(`--collect-all pyfxtran --onefile $(which pyfortool) src/pyfortool/scripting.py`),
but are built **inside Docker containers** that pin the target glibc, so the
binaries run on older systems than the runner itself:

- `pyfortool-x86_64`: built in `quay.io/pypa/manylinux2014_x86_64`
  (CentOS 7.9, **glibc 2.17**) on `ubuntu-latest`.
- `pyfortool-aarch64`: built in `quay.io/pypa/manylinux_2_28_aarch64`
  (RHEL/AlmaLinux 8, **glibc 2.28**) on the ARM runner `ubuntu-26.04-arm`.

Because CentOS 7 (glibc 2.17) is too old to run GitHub Actions' Node-based
actions, checkout and artifact upload run on the host runner and **only the
PyInstaller build runs inside the manylinux container**. The workflow bind-mounts
the repo into the container and runs the dedicated script
`.github/scripts/build_pyinstaller.sh` via `docker run` with the checked-out
repo mounted into `/src`.

The `manylinux` images ship static (wheel-only) CPythons under `/opt/python/`,
which PyInstaller can't use (`--onefile` needs a shared `libpython`), so the
script installs the **distribution's shared CPython** instead
(`dnf install python39` on the AlmaLinux 8 / glibc-2.28 ARM image, and the
SCL `rh-python38` on the CentOS 7 / glibc-2.17 x86 image). Because CentOS 7 is
EOL, the script first disables the leftover SCLo repositories and provides a
fresh repository pointing at the archived `vault.centos.org` mirrors before
installing. It then installs `pyfxtran` (`pyfortool`) and PyInstaller from
that interpreter, uses the C build toolchain (`gcc`/`make`/`git`) already
present, and `chown`s the build output back to the runner user so it can be
picked up on the host. Each job also verifies the binary's version and logs
the highest required GLIBC symbol (an `objdump`-based check) to confirm the
glibc bound.

## Version / tag consistency

The workflow refuses to proceed if the pushed tag does not equal the version
reported by `pyfortool.__version__`. This prevents publishing a release whose
artifacts report a different version than the release name.

## PyPI

The wheel is published to PyPI using the `PYPI_API_TOKEN` repository
secret set by the maintainers. Publication is gated behind the `publish`
environment in the repository settings.
