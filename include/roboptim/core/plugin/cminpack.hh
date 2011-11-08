
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

# include "roboptim/core/plugin/sum-of-c1-squares.hh"

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
      typedef Solver<SumOfC1Squares, boost::mpl::vector<> >
      parent_t;
      /// \brief Cost function type
      typedef problem_t::function_t function_t;
      /// \brief type of result
      typedef parent_t::result_t result_t;
      /// \brief type of gradient
      typedef DerivableFunction::gradient_t gradient_t;
      /// \brief Size type
      typedef Function::size_type size_type;

      /// \brief Constructot by problem
      explicit SolverWithJacobian (const problem_t& problem);
      virtual ~SolverWithJacobian () throw ();
      /// \brief Solve the optimization problem
      virtual void solve () throw ();
      /// Number of variables
      int n() const {return n_;}
      /// Number of functions
      int m() const {return m_;}
      /// Get parameter
      Function::argument_t& parameter () {
	return parameter_;
      }
      const Function::argument_t& parameter () const {return parameter_;}
      /// Get value
      const Function::argument_t& value () const {
	(*cost_)(value_, parameter_);
	return value_;
      }
      /// Get Jacobian
      const gradient_t& jacobianRow (size_t iRow) const {
	(*cost_).gradient(jacobianRow_, parameter_, iRow);
	return jacobianRow_;
      }

    private:
      /// Number of variables
      size_t n_;
      /// Dimension of the cost function
      size_t m_;
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
      Function::argument_t parameter_;
      /// Value of the function
      mutable Function::argument_t value_;
      /// Jacobian of the cost function
      mutable gradient_t jacobianRow_;
      /// Reference to cost function
      boost::shared_ptr <const DerivableFunction> cost_;
    }; // class Solver
  } // namespace cminpack
} // namespace roboptim
#endif // ROBOPTIM_CORE_PLUGIN_CMINPACK_CMINPACK_HH
