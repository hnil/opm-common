# - Find OPM Flow Diagnostics Library
#
# Defines the following variables:
#   opm-flowdiagnostics-applications_INCLUDE_DIRS    Directory of header files
#   opm-flowdiagnostics-applications_LIBRARIES       Directory of shared object files
#   opm-flowdiagnostics-applications_DEFINITIONS     Defines that must be set to compile
#   opm-flowdiagnostics-applications_CONFIG_VARS     List of defines that should be in config.h
#   HAVE_OPM_FLOWDIAGNOSTICS_APPLICATIONS            Binary value to use in config.h

# Copyright (C) 2012 Uni Research AS
# This code is licensed under The GNU General Public License v3.0

# use the generic find routine
include (opm-flowdiagnostics-applications-prereqs)

include (OpmPackage)

find_opm_package (
  # module name
  "opm-flowdiagnostics-applications"

  # dependencies
  "${opm-flowdiagnostics-applications_DEPS}"
  
  # header to search for
  "opm/utility/ECLUnitHandling.hpp"

  # library to search for
  "opmflowdiagnostics-applications"

  # defines to be added to compilations
  ""

  # test program
"#include <opm/utility/ECLUnitHandling.hpp>

int main()
{
    auto M = Opm::ECLUnits::createUnitSystem(1); // METRIC
}
"

  # config variables
  "${opm-flowdiagnostics-applications_CONFIG_VAR}"
  )
