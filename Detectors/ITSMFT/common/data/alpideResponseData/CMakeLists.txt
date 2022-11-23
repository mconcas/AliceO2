# Copyright 2019-2020 CERN and copyright holders of ALICE O2.
# See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
# All rights not expressly granted are reserved.
#
# This software is distributed under the terms of the GNU General Public
# License v3 (GPL Version 3), copied verbatim in the file "COPYING".
#
# In applying this license CERN does not waive the privileges and immunities
# granted to it by virtue of its status as an Intergovernmental Organization
# or submit itself to any jurisdiction.

message(STATUS "sourcedir is: " ${CMAKE_CURRENT_SOURCE_DIR})
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/AlpideResponseData.root)
# add_custom_command(
#     OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/AlpideResponseData.root
#     COMMAND $ENV{ROOTSYS}/bin/root -b -q  ${CMAKE_CURRENT_SOURCE_DIR}/alpideResponse.C++
#     DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/alpideResponse.C O2::ITSMFTSimulation
# )
add_custom_target(NAME ALL
                  COMMAND /usr/bin/env
                          #$ENV{ROOTSYS}/bin/root -b -q  ${CMAKE_CURRENT_SOURCE_DIR}/alpideResponse.C++
                  DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/alpideResponse.C
                          O2::ITSMFTSimulation
                  COMMENT "Generating ROOT file for Alpide response"
)

o2_add_test_root_macro(alpideResponse.C LABELS "its;mft")