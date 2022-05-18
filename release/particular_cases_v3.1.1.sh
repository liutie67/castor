#!/bin/bash

version="3.1.1"

# The analytical projection tool will be released later
rm castor_v${version}/castor-proj.cc
rm -r castor_v${version}/include/analytic_simulator castor_v${version}/src/analytic_simulator

# Tools that will be released later
rm castor_v${version}/toolkits/castor-mergeBedPositions.cc
rm castor_v${version}/toolkits/castor-makeReplicates.cc
rm castor_v${version}/toolkits/castor-sumImages.cc
