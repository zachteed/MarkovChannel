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
with inac.dat
```
-120.0000    1.0000
-110.0000    0.9796
-100.0000    0.9181
 -90.0000    0.7607
 -80.0000    0.4866
 -70.0000    0.2226
 -60.0000    0.0800
 -40.0000    0.0080
 -20.0000    0.0007
        0    0.0001
  20.0000         0
```

The header contains to following fields
* name - the name of the protocol
* source - the data used by the protocol
* v0 - the initial voltage
* normalize - whether to normalize output to [0, 1]

The protocol can then consist of any number of steps.  Each step contains the following fields:
* dt - the duration of the step
* vm - step voltage
* stype - the type of step this can be one of the following forms
	+ NONE - simply perform the step, produce no output
	+ PEAK - record the peak channel conductance
	+ MIN - record the minimum channel conductance
	+ TAU - record rate constants (parameterized by extra_args)
	+ TRACE - record channel conductance after each stepsize
* stepsize - optimal parameter for size of ODE/EXPM stepping
* extra_args - any number of additional step arguments, used for TAU stype

If either 'dt' or 'vm' is missing from the step, then the program treats this value as a variable and searches the .dat file for the missing paramters.  For example, in the inactivation protocol, the first step
```
step {
  dt: 200.0
  stype: NONE
}
```
is missing the 'vm' argument.  The .dat file is then searched and fills in 'vm' with the values of the first row in the .dat file.  So 'vm' becomes {-120, -110, -100, ..., 20}.

More examples of protocol encodings can be found in the demos folder.







