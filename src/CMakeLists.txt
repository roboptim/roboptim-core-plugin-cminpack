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

SET(PLUGIN_NAME ${PROJECT_NAME})

# Define the directory where plug-ins will be installed.
IF(ROBOPTIM_CORE_FOUND)
  GET_FILENAME_COMPONENT(RELPLUGINDIR ${ROBOPTIM_CORE_PLUGINDIR} NAME)
ELSE()
  GET_FILENAME_COMPONENT(RELPLUGINDIR ${ROBOPTIM_CORE_DEBUG_PLUGINDIR} NAME)
ENDIF(ROBOPTIM_CORE_FOUND)
SET(PLUGINDIR ${CMAKE_INSTALL_LIBDIR}/${RELPLUGINDIR})

ADD_LIBRARY(${PLUGIN_NAME} MODULE cminpack.cc)
SET_TARGET_PROPERTIES(${PLUGIN_NAME} PROPERTIES
  COMPILE_FLAGS "-Droboptim_core_plugin_cminpack_EXPORTS"
  PREFIX "")
IF(NOT APPLE)
  SET_TARGET_PROPERTIES(${PLUGIN_NAME}
    PROPERTIES VERSION 3.2.0 SOVERSION 3)
ENDIF()

PKG_CONFIG_USE_DEPENDENCY(${PLUGIN_NAME} roboptim-core)
PKG_CONFIG_USE_DEPENDENCY(${PLUGIN_NAME} cminpack)

INSTALL(TARGETS ${PLUGIN_NAME} DESTINATION ${PLUGINDIR})
