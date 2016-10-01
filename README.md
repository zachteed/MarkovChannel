# MarkovChannel
Ion channel Markov model parameter optimization

### Introduction

MarkovChannel is a tool for optimizing Markov models of ion channels, using the matrix exponential to greatly reduce fitting times.  It provides a flexable framework for encoding a range of experimentally collected protject into the model fitting.

### License

MarkovChannel is released under the MIT Licensse (refer to the LICENSE file for details).

### Contents
0. [Requirements: software](#requirements-software)
0. [Compilation](#compilation)
0. [Running the Demos](#running-the-demos)
0. [Understanding the Results](#understanding-to-results)
0. [Training on Other Data](#training-on-other-data)
0. [Troubleshooting](#troubleshooting)


### Requirements: software

0. 'Protobuf'
    - if you are using ubuntu run 'sudo apt-get install libprotobuf-dev protobuf-compiler'
0. 'mkl' Intel Math Kernel Library
    - [mkl](https://software.intel.com/en-us/intel-mkl)

### Compilation

0. To compile simply run 'make'

### Running the Demos

MarkovChannel comes with two demos; one for a Na<sup>+</sup> channel and one for a K<sup>+</sup> channel.

Before executing either of the examples, you must first set the mkl enviornment variables
    + run 'source /opt/intel/mkl/bin/mkvars.sh intel64'

To execute the Na^+ optimization
    + run './MarkovChannel solver.prototxt k-protocols.txt'

To exectute the K^+ optimization
    + run './MarkovChannel solver.prototxt na-protocols.txt'

When running these commands, optimization progress will be periodically displayed.  More detailed information and fitted models will be written to the snapshot directory.

### Understanding the Results

Results of the demos are written to the snapshot directory in the "K+\_Demo" and "Na+\_Demo" respectivly.  In each of these subdirectories, you will find iter\_(%d).txt and iter\_(%d).model.  The .txt files contain the model fits at that iteration for each of the protocols.  The .model file provides the model structure and rate parameters that determine the behavior of the model.

Included in MarkovChannel are some MATLAB scripts that can interpret the .model files.


### Training on Other Data

To fit models on other data, you must encode the experimental protocols in the .prototxt format.

For example, the Na+ inactivation protocol is represented as:
```
name: "inac"
source: "inac.dat"
v0: -120.0
normalize: true

step {
  dt: 200.0
  stype: NONE
}

step {
  dt: 2.5
  vm: -20
  stype: PEAK
  stepsizze: 0.05
}
```




