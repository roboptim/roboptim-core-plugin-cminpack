
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

#ifndef ROBOPTIM_CORE_PLUGIN_CMINPACK_CMINPACK_HH
# define ROBOPTIM_CORE_PLUGIN_CMINPACK_CMINPACK_HH

# include <boost/mpl/vector.hpp>
# include <roboptim/core/solver.hh>

# include "roboptim/core/sum-of-c1-squares.hh"

namespace roboptim {
  namespace cminpack {
    /// \brief Solver implementing a variant of Levenberg-Marquardt algorithm.
    ///
    /// This solver tries to minimize the euclidean norm of a vector valued
    /// function.
    class SolverWithJacobian :
      public Solver<SumOfC1Squares, boost::mpl::vector<> >
    {
    public:
      /// \brief Parent type
      typedef Solver<SumOfC1Squares, boost::mpl::vector<> > parent_t;

      /// \brief Cost function type
      typedef problem_t::function_t function_t;

      /// \brief Type of result
      typedef function_t::argument_t argument_t;
      typedef function_t::argument_ref argument_ref;
      typedef function_t::const_argument_ref const_argument_ref;

      /// \brief Type of result
      typedef function_t::result_t result_t;

      /// \brief Type of gradient
      typedef function_t::gradient_t gradient_t;
      typedef function_t::const_gradient_ref const_gradient_ref;

      /// \brief Size type
      typedef function_t::size_type size_type;

      /// \brief Constructot by problem
      explicit SolverWithJacobian (const problem_t& problem);
      virtual ~SolverWithJacobian ();
      /// \brief Solve the optimization problem
      virtual void solve ();

      /// Number of variables
      size_type n () const
      {
	return n_;
      }

      /// Number of functions
      size_type m () const
      {
	return m_;
      }

      /// Get parameter
      argument_t& parameter ()
      {
	return parameter_;
      }

      const argument_t& parameter () const
      {
	return parameter_;
      }

      /// Get value
      const argument_t& value () const
      {
	(*cost_)(value_, parameter_);
	return value_;
      }

      /// Get Jacobian
      const gradient_t& jacobianRow (size_type iRow) const
      {
	(*cost_).gradient (jacobianRow_, parameter_, iRow);
	return jacobianRow_;
      }

    private:
      /// Number of variables
      size_type n_;
      /// Dimension of the cost function
      size_type m_;
      /// Array of double to store variable of optimization problem
      double* x_;
      /// Array of double to store value of optimization problem
      double* fvec_;
      /// Array of double to store one line of Jacobian
      double* fjac_;
      /// Array of int used by the optimizer
      int* ipvt_;
      /// Positive integer not less than 5*n_+m_
      int lwa_;
      /// array of double of size lwa_;
      double* wa_;

      /// Parameter of the function
      argument_t parameter_;
      /// Value of the function
      mutable argument_t value_;
      /// Jacobian of the cost function
      mutable gradient_t jacobianRow_;
      /// Reference to cost function
      boost::shared_ptr <const DifferentiableFunction> cost_;
    }; // class Solver
  } // namespace cminpack
} // namespace roboptim
#endif // ROBOPTIM_CORE_PLUGIN_CMINPACK_CMINPACK_HH
