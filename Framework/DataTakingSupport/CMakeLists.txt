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
o2_add_library(FrameworkDataTakingSupport
               SOURCES src/Plugin.cxx
               PRIVATE_INCLUDE_DIRECTORIES ${CMAKE_CURRENT_LIST_DIR}/src
               PUBLIC_LINK_LIBRARIES O2::Framework AliceO2::InfoLogger
               PRIVATE_LINK_LIBRARIES $<$<TARGET_EXISTS:AliceO2::Tracing>:AliceO2::Tracing>
               TARGETVARNAME targetName)

# Propagate tracing support so Plugin.cxx can inject traceid/spanid into InfoLogger messages.
if(TARGET AliceO2::Tracing)
  target_compile_definitions(${targetName} PRIVATE O2_WITH_DPL_TRACING)
endif()
