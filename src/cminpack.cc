
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

#include <cstring>

#include <cminpack.h>

#include "roboptim/core/function.hh"
#include "roboptim/core/problem.hh"
#include "roboptim/core/solver-error.hh"
#include "roboptim/core/plugin/cminpack.hh"


namespace roboptim {
  namespace cminpack {
    // Copy row iRow of matrix to array of double
    void
    matrix_to_array (Function::value_type* dst, const Function::matrix_t& src,
		     int iRow);
  } // namespace cminpack
} // end of namespace roboptim

using roboptim::detail::array_to_vector;
using roboptim::detail::vector_to_array;
using roboptim::cminpack::matrix_to_array;

extern "C" {
  int  roboptim_plugin_cminpack_fcn
  (void *p, int m, int n, const double *x, double *fvec, double *fjrow,
   int iflag)
  {
    // Get instance of solver.
    roboptim::cminpack::SolverWithJacobian* solver =
      static_cast <roboptim::cminpack::SolverWithJacobian*> (p);
    assert (n == solver->n ());
    assert (m == solver->m ());
    // Get parameter as a roboptim vector
    array_to_vector (solver->parameter (), x);
    // Get cost value and convert to array
    vector_to_array (fvec, solver->value());
    if (iflag < 2) return 0;
    // Get cost jacobian and convert to array
    ::roboptim::Function::size_type row =
	static_cast< ::roboptim::Function::size_type> (iflag) - 2;
    vector_to_array (fjrow, solver->jacobianRow (row));
    return 0;
  }
}

namespace roboptim
{
  namespace cminpack
  {
    SolverWithJacobian::SolverWithJacobian (const problem_t& problem) :
      Solver <SumOfC1Squares, boost::mpl::vector<> >
      (problem),
      n_ (problem.function ().baseFunction ()->inputSize ()),
      m_ (problem.function ().baseFunction ()->outputSize ()),
      x_ (new double [n_]),
      fvec_ (new double [m_]),
      fjac_ (new double [n_*n_]),
      ipvt_ (new int [n_]),
      lwa_ (static_cast<int> (5 * n_ + m_)),
      wa_ (new double [lwa_]),
      parameter_ (n_),
      value_ (m_),
      jacobianRow_ (n_),
      cost_ (problem.function ().baseFunction ())
    {
      std::size_t n = static_cast<std::size_t> (n_);
      std::size_t m = static_cast<std::size_t> (m_);
      std::size_t lwa = static_cast<std::size_t> (lwa_);

      // Initialize memory
      memset (x_, 0, n * sizeof (double));
      memset (fvec_, 0, m * sizeof (double));
      memset (fjac_, 0, n * n * sizeof (double));
      memset (ipvt_, 0, n * sizeof (int));
      memset (wa_, 0, lwa * sizeof (double));

      // Initialize this class parameters
      parameter_.setZero ();
      value_.setZero ();
      jacobianRow_.setZero ();
    }

    SolverWithJacobian::~SolverWithJacobian () throw ()
    {
      delete[] x_;
      delete[] fvec_;
      delete[] fjac_;
      delete[] ipvt_;
      delete[] wa_;
    }

    void SolverWithJacobian::solve () throw ()
    {
      int ldfjac = static_cast<int> (n_);
      double tol = 1e-6;
      // Set initial guess
      if (problem().startingPoint()) {
	vector_to_array (x_, *(problem().startingPoint()));
      }
      int info = lmstr1(roboptim_plugin_cminpack_fcn,
			(void*)this,
			static_cast<int> (m_), static_cast<int> (n_),
			x_, fvec_, fjac_, ldfjac,
			tol, ipvt_, wa_, lwa_);
      switch (info) {
      case 0:
	result_ = SolverError ("improper input parameters");
	break;
      case 1:
      case 2:
      case 3:
	{
	  Result result (n_, 1);
	  array_to_vector(result.x, x_);
	  result.value = problem().function()(result.x);
	  result_ = result;
	}
	break;
      case 4:
	{
	  ResultWithWarnings result (n_, 1);
	  array_to_vector(result.x, x_);
	  result.value = problem().function()(result.x);
	  result.warnings.push_back(SolverWarning
				    ("fvec is orthogonal to the columns of"
				     " the jacobian to machine precision"));
	  result_ = result;
	}
	break;
      case 5:
	{
	  ResultWithWarnings result (n_, 1);
	  array_to_vector(result.x, x_);
	  result.value = problem().function()(result.x);
	  result.warnings.push_back(SolverWarning
				    ("number of calls to fcn with iflag = 1 "
				     "has reached 100*(n+1)"));
	  result_ = result;
	}
	break;
      case 6:
	{
	  ResultWithWarnings result (n_, 1);
	  array_to_vector(result.x, x_);
	  result.value = problem().function()(result.x);
	  result.warnings.push_back(SolverWarning
				    ("tol is too small. no further reduction"
				     " in the sum of squares is possible"));
	  result_ = result;
	}
	break;
      case 7:
	{
	  ResultWithWarnings result (n_, 1);
	  array_to_vector(result.x, x_);
	  result.value = problem().function()(result.x);
	  result.warnings.push_back(SolverWarning
				    ("tol is too small. no further"
				     " improvement in the approximate"
				     " solution x is possible"));
	  result_ = result;
	}
	break;
      default:
	result_ = SolverError ("Return value not documented");
      }
    }

    void
    matrix_to_array (Function::value_type* dst, const Function::matrix_t& src,
		     int iRow)
    {
      if (src.cols () == 0)
	return;

      std::size_t size =
	static_cast<std::size_t> (src.cols ()) * sizeof (Function::value_type);
      memcpy (dst, &src(iRow,0), size);

      // NaN != NaN, handle this case.
      for (Eigen::VectorXd::Index j = 0; j < src.cols (); ++j)
	if (src(iRow,j) != src(iRow,j))
	  assert (dst[j] != dst[j]);
	else
	  assert (dst[j] == src(iRow,j));
    }
  } // namespace cminpack
} // end of namespace roboptim

extern "C"
{
  using namespace roboptim::cminpack;
  typedef SolverWithJacobian::parent_t solver_t;

  ROBOPTIM_DLLEXPORT unsigned getSizeOfProblem ();
  ROBOPTIM_DLLEXPORT solver_t* create (const SolverWithJacobian::problem_t& pb);
  ROBOPTIM_DLLEXPORT void destroy (solver_t* p);

  ROBOPTIM_DLLEXPORT unsigned getSizeOfProblem ()
  {
    return sizeof (solver_t::problem_t);
  }

  ROBOPTIM_DLLEXPORT solver_t* create (const SolverWithJacobian::problem_t& pb)
  {
    return new SolverWithJacobian (pb);
  }

  ROBOPTIM_DLLEXPORT void destroy (solver_t* p)
  {
    delete p;
  }
}
