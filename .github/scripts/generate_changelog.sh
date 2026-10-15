#!/usr/bin/env bash
#
# Generate a Markdown changelog for the given release tag.
#
# The changelog lists, for every commit between the previous version tag and
# the current one, the commit title, author and date.
#
# Usage:
#   generate_changelog.sh <tag> > CHANGELOG.md
#
set -euo pipefail

TAG="${1:?usage: generate_changelog.sh <tag>}"
DATE="$(date -u +%Y-%m-%d)"

# Find the nearest previous version tag that is an ancestor of the target tag.
# Tags are bare version numbers like 0.2.21.
TAG_COMMIT="$(git rev-list -n 1 "${TAG}")"
PREV_TAG="$(git tag --sort=-creatordate | while read -r t; do
    [ "${t}" = "${TAG}" ] && continue
    if git merge-base --is-ancestor "${t}" "${TAG_COMMIT}" 2>/dev/null; then
        printf '%s' "${t}"
        break
    fi
done)"

RANGE="${PREV_TAG:+${PREV_TAG}..}${TAG}"

echo "# pyfortool ${TAG}"
echo
echo "**Release date:** ${DATE}"
echo

# List each commit as: - title — author (date)
while IFS='|' read -r subject author commit_date; do
    [ -z "${subject}" ] && continue
    printf -- '- %s — %s (%s)\n' "${subject}" "${author}" "${commit_date}"
done < <(git log --no-merges --pretty=format:'%s|%an|%ad' --date=short \
    "${RANGE}" 2>/dev/null)

echo
