
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
#include "roboptim/core/plugin/cminpack/config.hh"

using roboptim::detail::array_to_vector;
using roboptim::detail::vector_to_array;

extern "C" {
  int  roboptim_plugin_cminpack_fcn
  (void *p, int ROBOPTIM_DEBUG_ONLY(m), int ROBOPTIM_DEBUG_ONLY(n),
   const double *x, double *fvec, double *fjrow, int iflag)
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
    SolverWithJacobian::SolverWithJacobian (const problem_t& pb) :
      parent_t (pb),
      n_ (),
      m_ (),
      x_ (),
      fvec_ (),
      fjac_ (),
      ipvt_ (),
      lwa_ (),
      wa_ (),
      parameter_ (),
      value_ (),
      jacobianRow_ (),
      baseCost_ ()
    {
      const SumOfC1Squares* cost
        = dynamic_cast<const SumOfC1Squares*> (&pb.function ());

      if (!cost)
      {
        throw std::runtime_error ("the cminpack plugin expects"
                                  " a SumOfC1Squares cost function");
      }

      baseCost_ = cost->baseFunction ();

      n_ = baseCost_->inputSize ();
      m_ = baseCost_->outputSize ();
      lwa_ = static_cast<int> (5 * n_ + m_);

      // Initialize memory
      x_.resize (n_);
      x_.setZero ();
      fvec_.resize (m_);
      fvec_.setZero ();
      fjac_.resize (n_*n_);
      fjac_.setZero ();
      ipvt_.resize (n_);
      ipvt_.setZero ();
      wa_.resize (lwa_);
      wa_.setZero ();

      // Initialize this class parameters
      parameter_.resize (n_);
      parameter_.setZero ();
      value_.resize (m_);
      value_.setZero ();
      jacobianRow_.resize (n_);
      jacobianRow_.setZero ();
    }

    SolverWithJacobian::~SolverWithJacobian ()
    {
    }

    void SolverWithJacobian::solve ()
    {
      int ldfjac = static_cast<int> (n_);
      double tol = 1e-6;

      // Set initial guess
      if (problem().startingPoint()) {
        x_ = *(problem().startingPoint());
      }

      int info = lmstr1(roboptim_plugin_cminpack_fcn,
			(void*)this,
			static_cast<int> (m_), static_cast<int> (n_),
			x_.data(), fvec_.data(), fjac_.data(), ldfjac,
			tol, ipvt_.data(), wa_.data(), lwa_);
      switch (info) {
      case 0:
	result_ = SolverError ("improper input parameters");
	break;
      case 1:
      case 2:
      case 3:
	{
	  Result result (n_, 1);
	  result.x = x_;
	  result.value = problem().function()(result.x);
	  result_ = result;
	}
	break;
      case 4:
	{
	  Result result (n_, 1);
	  result.x = x_;
	  result.value = problem().function()(result.x);
	  result.warnings.push_back(SolverWarning
				    ("fvec is orthogonal to the columns of"
				     " the jacobian to machine precision"));
	  result_ = result;
	}
	break;
      case 5:
	{
	  Result result (n_, 1);
	  result.x = x_;
	  result.value = problem().function()(result.x);
	  result.warnings.push_back(SolverWarning
				    ("number of calls to fcn with iflag = 1 "
				     "has reached 100*(n+1)"));
	  result_ = result;
	}
	break;
      case 6:
	{
	  Result result (n_, 1);
	  result.x = x_;
	  result.value = problem().function()(result.x);
	  result.warnings.push_back(SolverWarning
				    ("tol is too small. no further reduction"
				     " in the sum of squares is possible"));
	  result_ = result;
	}
	break;
      case 7:
	{
	  Result result (n_, 1);
	  result.x = x_;
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

  } // namespace cminpack
} // end of namespace roboptim

extern "C"
{
  using namespace roboptim::cminpack;
  typedef SolverWithJacobian::parent_t solver_t;

  ROBOPTIM_CORE_PLUGIN_CMINPACK_DLLEXPORT unsigned getSizeOfProblem ();
  ROBOPTIM_CORE_PLUGIN_CMINPACK_DLLEXPORT const char* getTypeIdOfConstraintsList ();
  ROBOPTIM_CORE_PLUGIN_CMINPACK_DLLEXPORT solver_t* create (const SolverWithJacobian::problem_t& pb);
  ROBOPTIM_CORE_PLUGIN_CMINPACK_DLLEXPORT void destroy (solver_t* p);

  ROBOPTIM_CORE_PLUGIN_CMINPACK_DLLEXPORT unsigned getSizeOfProblem ()
  {
    return sizeof (solver_t::problem_t);
  }

  ROBOPTIM_CORE_PLUGIN_CMINPACK_DLLEXPORT const char* getTypeIdOfConstraintsList ()
  {
    return typeid (solver_t::problem_t::constraintsList_t).name ();
  }

  ROBOPTIM_CORE_PLUGIN_CMINPACK_DLLEXPORT solver_t* create (const SolverWithJacobian::problem_t& pb)
  {
    return new SolverWithJacobian (pb);
  }

  ROBOPTIM_CORE_PLUGIN_CMINPACK_DLLEXPORT void destroy (solver_t* p)
  {
    delete p;
  }
}
