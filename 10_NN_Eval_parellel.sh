#!/bin/bash

waitall() { # PID...
  ## Wait for children to exit and indicate whether all exited with 0 status.
  local errors=0
  while :; do
    #debug "Processes remaining: $*"
    for pid in "$@"; do
      shift
      if kill -0 "$pid" 2>/dev/null; then
        #debug "$pid is still alive."
        set -- "$@" "$pid"
      elif wait "$pid"; then
        debug "$pid exited with zero exit status."
      else
        debug "$pid exited with non-zero exit status."
        ((++errors))
      fi
    done
    (("$#" > 0)) || break
    # TODO: how to interrupt this sleep when a child terminates?
    sleep ${WAITALL_DELAY:-1}
   done
  ((errors == 0))
}

debug() { echo "DEBUG: $*" >&2; }

numbOfProc=3

tissues=('zero'
BRAIN
BRAIN
BRAIN
FB
FB
FB
GUT
GUT
GUT
MT
MT
MT
BRAIN
BRAIN
BRAIN
FB
FB
FB
GUT
GUT
GUT
MT
MT
MT
)

outs=('zero'
Without
Without
Without
Without
Without
Without
Without
Without
Without
Without
Without
Without
With
With
With
With
With
With
With
With
With
With
With
With
)

schemes=('zero'
n12_6_3
n12_6
n12
n12_6_3
n12_6
n12
n12_6_3
n12_6
n12
n12_6_3
n12_6
n12
n12_6_3
n12_6
n12
n12_6_3
n12_6
n12
n12_6_3
n12_6
n12
n12_6_3
n12_6
n12
)

echo "Here"

sampleCount=${#tissues[@]}
for (( i=1; i<$sampleCount; i+=$numbOfProc )); do
   # bash arrays are 0-indexed
   pids=""
   endInd=$(($numbOfProc - 1))

   for t in $(seq 0 $endInd); do
      position=$(( $i + $t ))

      # Get current sample
      currTiss=${tissues[$position]}
      currOut=${outs[$position]}
      currScheme=${schemes[$position]}

      echo $currTiss $currOut $currScheme
      Rscript --vanilla 10_NN_Eval_parellel.R $currTiss $currOut $currScheme 1>$currTiss "_"$currOut"_"$currScheme".out" 2>$currTiss"_"$currOut"_"$currScheme".err" &

      pids="$pids $!"
   done
   waitall $pids
done

