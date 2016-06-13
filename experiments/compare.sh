#!/bin/bash

export OMP_NUM_THREADS=1
mkdir -p experiments/log

for i in `seq 1 10`; do
  ./MarkovChannel experiments/config/solver1.prototxt > \
    experiments/log/1.$i.txt &
  sleep 1
done

for i in `seq 1 10`; do
  ./MarkovChannel experiments/config/solver2.prototxt > \
    experiments/log/2.$i.txt &
  sleep 1
done

for i in `seq 1 10`; do
  ./MarkovChannel experiments/config/solver3.prototxt > \
    experiments/log/3.$i.txt &
  sleep 1
done



