#
# Copyright (c) 2011 CNRS
# Authors: Florent Lamiraux
#
#
# This file is part of roboptim-core-plugin-cminpack
# roboptim-core-plugin-cminpack is free software: you can redistribute it
# and/or modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation, either version
# 3 of the License, or (at your option) any later version.
#
# roboptim-core-plugin-cminpack is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied warranty
# of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Lesser Public License for more details.  You should have
# received a copy of the GNU Lesser General Public License along with
# roboptim-core-plugin-cminpack  If not, see
# <http://www.gnu.org/licenses/>.

CONFIG_FILES(shared-tests/local-libdir.hh)

# Make sure local headers are found
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}/tests)
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/tests)

# Add Boost path to include directories.
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

# Make Boost.Test generates the main function in test cases.
ADD_DEFINITIONS(-DBOOST_TEST_DYN_LINK -DBOOST_TEST_MAIN)

SET(COMMON_HEADERS
  ${HEADERS}
  shared-tests/common.hh
  shared-tests/distance-to-sphere.hh
  )

# ROBOPTIM_CORE_CMINPACK_PLUGIN_TEST(NAME)
# ----------------------------------------
#
# Define a test named `NAME'.
#
# This macro will create a binary from `NAME.cc', link it against
# Boost and add it to the test suite.
#
MACRO(ROBOPTIM_CORE_PLUGIN_CMINPACK_TEST NAME)
  ADD_EXECUTABLE(${NAME}
    ${NAME}.cc
    ${COMMON_HEADERS})
  ADD_TEST(${NAME} ${RUNTIME_OUTPUT_DIRECTORY}/${NAME})

  PKG_CONFIG_USE_DEPENDENCY(${NAME} roboptim-core)
  PKG_CONFIG_USE_DEPENDENCY(${NAME} cminpack)

  # Link against Boost.
  TARGET_LINK_LIBRARIES(${NAME} ${Boost_LIBRARIES})
ENDMACRO(ROBOPTIM_CORE_PLUGIN_CMINPACK_TEST)


# simple
# Link the dummy solver with this test to bypass the plug-in system.
ADD_EXECUTABLE(simple
  simple.cc
  ${COMMON_HEADERS}
  ${CMAKE_SOURCE_DIR}/src/cminpack.cc
  )
PKG_CONFIG_USE_DEPENDENCY(simple roboptim-core)
PKG_CONFIG_USE_DEPENDENCY(simple cminpack)
TARGET_LINK_LIBRARIES(simple ${Boost_LIBRARIES})
ADD_TEST(simple ${RUNTIME_OUTPUT_DIRECTORY}/simple)

ROBOPTIM_CORE_PLUGIN_CMINPACK_TEST(plugin)
