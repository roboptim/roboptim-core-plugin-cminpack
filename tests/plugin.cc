
// Copyright (c) 2011 CNRS
// Authors: Florent Lamiraux


// This file is part of roboptim-core-plugin-cminpack
// roboptim-core-plugin-cminpack is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.

// roboptim-core-plugin-cminpack is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// roboptim-core-plugin-cminpack  If not, see
// <http://www.gnu.org/licenses/>.

#include <roboptim/core/solver.hh>
#include <roboptim/core/solver-error.hh>
#include <roboptim/core/solver-factory.hh>

#include "roboptim/core/plugin/cminpack.hh"
#include "shared-tests/common.hh"
#include "shared-tests/distance-to-sphere.hh"

using roboptim::SumOfC1Squares;
using roboptim::Result;
using roboptim::GenericSolver;
using roboptim::Solver;
using roboptim::SolverError;
using roboptim::SolverFactory;
using roboptim::cminpack::SolverWithJacobian;

typedef Solver <SumOfC1Squares, boost::mpl::vector <> > solver_t;

BOOST_AUTO_TEST_CASE (plugin)
{
  lt_dlinit();
  BOOST_REQUIRE_EQUAL (lt_dlsetsearchpath (PLUGIN_PATH), 0);

  boost::shared_ptr<boost::test_tools::output_test_stream>
    output = retrievePattern ("plugin");

  boost::shared_ptr <roboptim::cminpack::F> f (new roboptim::cminpack::F());
  SumOfC1Squares soq (f, "");
  solver_t::problem_t pb(soq);
  roboptim::cminpack::initialize_problem <SolverWithJacobian::problem_t>(pb);

  SolverFactory <solver_t> factory ("cminpack", pb);
  solver_t& solver = factory ();

  // Solve the problem and get result
  SolverWithJacobian::result_t res = solver.minimum ();
  // Display solver information.
  (*output) << solver << std::endl;

  if ((res.which () != GenericSolver::SOLVER_VALUE) &&
      (res.which () != GenericSolver::SOLVER_VALUE_WARNINGS))
    {
      (*output) << "A solution should have been found. Failing..."
		<< std::endl
		<< boost::get<SolverError> (res).what ()
		<< std::endl;
      std::cout << output->str () << std::endl;
      BOOST_FAIL ("optimization failed");
    }

  // Get the result.
  Result& result = boost::get<Result> (res);

  // Display the result.
  (*output) << "A solution has been found: " << std::endl;
  (*output) << result << std::endl;

  std::cout << output->str () << std::endl;
  BOOST_CHECK (output->match_pattern ());

  if (lt_dlexit ())
    std::cerr << "lt_dlexit failed" << std::endl;
}
