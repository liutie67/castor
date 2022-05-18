#!/bin/bash

version="2.1"

# The analytical projection tool will be released later
rm castor_v${version}/castor-proj.cc
rm -r castor_v${version}/include/analytic_simulator castor_v${version}/src/analytic_simulator

# Parts of the documentation that will be released later
rm castor_v${version}/docs/CASToR__image_deformation.pdf
rm castor_v${version}/docs/CASToR__dynamic_model.pdf

# Tools that will be released later
rm castor_v${version}/toolkits/castor-mergeBedPositions.cc
rm castor_v${version}/toolkits/castor-makeReplicates.cc
rm castor_v${version}/toolkits/castor-sumImages.cc
rm castor_v${version}/toolkits/castor-ImageBasedDynamicModel.cc

# Dynamic classes will be released later
rm castor_v${version}/src/dynamic/iLinearPatlakModel.cc
rm castor_v${version}/include/dynamic/iLinearPatlakModel.hh
rm castor_v${version}/src/image/iDeformationRigid.cc
rm castor_v${version}/include/image/iDeformationRigid.hh

