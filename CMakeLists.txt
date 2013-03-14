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

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

#SET(CXX_DISABLE_WERROR TRUE)
INCLUDE(cmake/base.cmake)
INCLUDE(cmake/boost.cmake)
INCLUDE(cmake/cpack.cmake)

SET(PROJECT_NAME roboptim-core-plugin-cminpack)
SET(PROJECT_DESCRIPTION
  "Non linear least-square optimization using Levenberg Marquardt algorithm."
  )
SET(PROJECT_URL "")
SET(CUSTOM_HEADER_DIR "roboptim/core/plugin")

# Where to compile shared objects
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

SETUP_PROJECT()

# Declare headers
SET(${PROJECT_NAME}_HEADERS
  include/roboptim/core/plugin/cminpack.hh
  )

# Search for dependencies.
SEARCH_FOR_BOOST()

ADD_REQUIRED_DEPENDENCY("cminpack >= 1.0")
ADD_REQUIRED_DEPENDENCY("roboptim-core >= 0.5")

# Add dependency toward pg-slider library in pkg-config file.
PKG_CONFIG_APPEND_LIBS("roboptim-core-plugin-cminpack")

ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(tests)

SETUP_PROJECT_FINALIZE()
