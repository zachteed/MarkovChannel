#!/bin/bash

export OMP_NUM_THREADS=1
mkdir -p experiments/log

for i in `seq 1 1`; do
  ./MarkovChannel experiments/config/solver_expm.prototxt > \
    experiments/log/expm_$i.txt &
done

for i in `seq 1 1`; do
  ./MarkovChannel experiments/config/solver_ode.prototxt > \
    experiments/log/ode_$i.txt &
done


