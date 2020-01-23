#!/bin/bash

./readFitness.py ${BASH_ARGV[0]}
./plotDevice.py ${BASH_ARGV[0]}
./convergence.py ${BASH_ARGV[0]}