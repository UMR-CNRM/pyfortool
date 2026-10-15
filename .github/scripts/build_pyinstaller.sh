#!/usr/bin/env bash
#
# Build a one-file PyInstaller binary of pyfortool inside a manylinux Docker
# image so the resulting binary is pinned to the target glibc of that image.
#
# Usage:
#   build_pyinstaller.sh <runs-in-container> <arch>
#
# The script is designed to run INSIDE the manylinux container with the
# repository bind-mounted at /src. It installs a distribution (shared) CPython
# and uses PyInstaller to build the binary, then chowns the build output back
# to the host runner user.
#
# Arguments:
#   $1 : the host user id (to chown build outputs back to)
#   $2 : the host group id
#   $3 : architecture suffix for the output binary (e.g. x86_64, aarch64)

set -e

HOST_UID="$1"
HOST_GID="$2"
ARCH="$3"

cd /src

# Install a distribution CPython built with --enable-shared. The manylinux
# images ship static (wheel-only) interpreters under /opt/python, which
# PyInstaller cannot use for --onefile builds (it needs libpython shared).
if [ -x /usr/bin/dnf ]; then
    # AlmaLinux 8 based (glibc 2.28)
    DNF_CMD="/usr/bin/dnf"
    if command -v dnf >/dev/null 2>&1; then
        DNF_CMD="$(command -v dnf)"
    fi
    "$DNF_CMD" -y install python39 python39-pip
    PYEXE=/usr/bin/python3.9
else
    # CentOS 7 based (glibc 2.17). CentOS 7 is EOL, so the default SCLo
    # repositories still point at the dead mirrorlist.centos.org. Disable
    # them and provide a vault-based repository for rh-python38.
    yum -y install centos-release-scl 2>/dev/null || true
    sed -i "s/^enabled=1/enabled=0/" /etc/yum.repos.d/CentOS-SCLo-*.repo 2>/dev/null || true
    printf "%s\n" \
        "[vault-sclo-rh]" \
        "name=Vault CentOS-7 SCLo rh (glibc 2.17)" \
        "baseurl=http://vault.centos.org/7.9.2009/sclo/\$basearch/rh/" \
        "enabled=1" \
        "gpgcheck=0" \
        > /etc/yum.repos.d/vault-sclo-rh.repo
    yum -y install rh-python38 rh-python38-python-pip
    # shellcheck disable=SC1091
    . /opt/rh/rh-python38/enable
    PYEXE=python
fi

"$PYEXE" --version
"$PYEXE" -m pip --version
"$PYEXE" -m pip install --upgrade pip
"$PYEXE" -m pip install -e . pyinstaller
"$PYEXE" -m PyInstaller --collect-all pyfxtran --onefile "$(command -v pyfortool)" src/pyfortool/scripting.py

# The container runs as root; hand the build output back to the host runner
# user so it can be renamed/uploaded by the workflow.
chown -R "${HOST_UID}:${HOST_GID}" /src/dist /src/*.spec 2>/dev/null || true
