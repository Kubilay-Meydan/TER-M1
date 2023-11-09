#!/usr/bin/env bash
set -o errexit   # abort on nonzero exitstatus
set -o nounset   # abort on unbound variable
set -o pipefail  # don't hide errors within pipes

variables(){
   in_salmon="${1}"
   out_sym="${2}"
}

main(){
    variables $@
    ln --force --relative --symbolic $in_salmon $out_sym
}

main "$@"
