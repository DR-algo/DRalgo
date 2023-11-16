#!/bin/sh
#
# @file update-copyright.sh
#
# Updates the copyright notices, like
#
#     Copyright (C) 2021-2023 Andreas Ekstedt
#     Copyright (C) 2021-2023 Philipp Schicho
#     Copyright (C) 2021-2023 Tuomas V.I. Tenkanen
#
# in all files in the current directory and its subdirectories.
#
# NOTE: GNU grep/sed required.
#
# Usage:
#   update-copyright.sh
#
set -eu

year=$(date +%Y)

grep -l -r 'Copyright *(C).*Ekstedt' . | while IFS= read -r f; do
  case $f in
    *update-copyright.sh)
      ;;
    *)
      sed -i "s/Copyright *(C).*Ekstedt/Copyright (C) 2021-$year Andreas Ekstedt/g" "$f"
      ;;
  esac
done
grep -l -r 'Copyright *(C).*Schicho' . | while IFS= read -r f; do
  case $f in
    *update-copyright.sh)
      ;;
    *)
      sed -i "s/Copyright *(C).*Schicho/Copyright (C) 2021-$year Philipp Schicho/g" "$f"
      ;;
  esac
done
grep -l -r 'Copyright *(C).*Tenkanen' . | while IFS= read -r f; do
  case $f in
    *update-copyright.sh)
      ;;
    *)
      sed -i "s/Copyright *(C).*Tenkanen/Copyright (C) 2021-$year Tuomas V.I. Tenkanen/g" "$f"
      ;;
  esac
done