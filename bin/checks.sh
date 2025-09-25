#!/bin/bash

#set -x
set -e
set -o pipefail #abort if left command on a pipe fails

function usage {
  echo "Usage: $0 [-h] [-s] [test1 [test2 [...]]]"
  echo "  -h   Prints this help message"
  echo "  -s   Makes this script silent (and only return a status)"
  echo "  -v   Print output of tests"
  echo ""
  echo "If the script is called without specifying tests, all tests are performed."
}

silence=0
tests=""
verbose=0
while [ -n "$1" ]; do
  case "$1" in
    '-h') usage; exit;;
    '-s') silence=1;;
    '-v') verbose=1;;
    *) tests="${tests} $1";;
  esac
  shift
done

PYFTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/.."
cd ${PYFTDIR}

global_ret=0

##### Check 'version' consistency
if [ "${tests}" == "" -o "$(echo ${tests} | grep -w "version")" != "" ]; then
  set +e
  vinit=$(python3 -c "$(grep __version__ src/pyfortool/__init__.py); print(__version__)")
  retval1=$?
  vtag=$(git describe --abbrev=0)
  retval2=$?
  set -e
  if [ $retval1 -ne 0 -o $retval2 -ne 0 ]; then
    ret=1
    result="\e[31mERROR\e[0m"
    [ ${silence} == 1 ] && echo "ERROR during version check"
  else
    if [ "${vinit}" == "${vtag}" ]; then
      ret=0
      result="\e[32mOK\e[0m"
    else
      ret=1
      result="\e[31mproblem\e[0m"
    fi
  fi

  global_ret=$((${global_ret} + ${ret}))
  [ ${silence} == 0 ] && echo -e "  'version' consistency: ${result}"
fi

##### Check pylint score
if [ "${tests}" == "" -o "$(echo ${tests} | grep -w "pylint")" != "" ]; then
  set +e
  output=$(pylint -d R0912,C0209,R0915,R1702,C0302,R0913,R0914,W1202,R0904,R0902 \
                  --persistent=n -f parseable src/pyfortool/)
  set -e
  score=$(echo "${output}" | grep 'Your code has been rated at' | cut -d\  -f 7 | cut -d/ -f 1)
  [ ${verbose} == 1 ] && echo "${output}"
  if [ $(python3 -c "print(0 if ${score} >= 9.8 else 1)") -ne 0 ]; then
    ret=1
    result="\e[31mproblem\e[0m"
  else
    ret=0
    result="\e[32mOK\e[0m"
  fi

  global_ret=$((${global_ret} + ${ret}))
  [ ${silence} == 0 ] && echo -e "  pylint score: ${result}"
fi

##### Check flake8 score
if [ "${tests}" == "" -o "$(echo ${tests} | grep -w "flake8")" != "" ]; then
  set +e
  output=$(flake8 src/pyfortool/)
  retval=$?
  set -e
  score=$(echo ${output} | wc -l)
  [ ${verbose} == 1 ] && echo "${output}"
  if [ $retval -ne 0 ]; then
    ret=1
    if [ ${score} -ne 0 ]; then
      result="\e[31mproblem\e[0m"
    else
      result="\e[31mERROR\e[0m"
      [ ${silence} == 1 ] && echo "ERROR during flake8 execution"
    fi
  else
    ret=0
    result="\e[32mOK\e[0m"
  fi

  global_ret=$((${global_ret} + ${ret}))
  [ ${silence} == 0 ] && echo -e "  flake8 score: ${result}"
fi

##### Check examples
if [ "${tests}" == "" -o "$(echo ${tests} | grep -w "examples")" != "" ]; then
  cd examples
  set +e
  output=$(./tests.sh 2>&1)
  retval=$?
  set -e
  [ ${verbose} == 1 ] && echo "${output}"
  if [ ${retval} -ne 0 ]; then
    ret=1
    result="\e[31mproblem\e[0m"
  else
    ret=0
    result="\e[32mOK\e[0m"
  fi

  global_ret=$((${global_ret} + ${ret}))
  [ ${silence} == 0 ] && echo -e "  test examples: ${result}"
fi

exit ${global_ret}
