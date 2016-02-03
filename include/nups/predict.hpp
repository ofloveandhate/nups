// This file is part of Numerical Univaraite Polynomial Solver (NUPS).

// NUPS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// NUPS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with NUPS.  If not, see <http://www.gnu.org/licenses/>.

// Copyright Daniel Brake, 2016
// University of Notre Dame
// Applied and Computational Mathematics and Statistics
// danielthebrake@gmail.com

// provides numerical predictor types for NUPS

#ifndef NUPS_PREDICT_HPP
#define NUPS_PREDICT_HPP



namespace nups {

	namespace predict {

		template<class PolyT>
		struct PredictorBase
		{
			
			template<typename NumT>
			static void Predict(std::vector<NumT> & solution, std::vector<NumT> const& variables, std::vector<NumT> const& rhs, NumT const& t,  NumT const& delta_t)
			{
				return PolyT::DoPredict(solution, variables, rhs, t, delta_t);
			}
		};



		template<typename LinearSolverT>
		struct RK4 : public PredictorBase<RK4<LinearSolverT> >
		{


			/**
			Fourth-order Runge-Kutta predictor, using a linear solver provided as a template type.
			*/
			template<typename NumT>
			static void DoPredict(std::vector<NumT> & solution, std::vector<NumT> const& variables, std::vector<NumT> const& rhs, NumT const& t,  NumT const& delta_t)
			{

			}
		};

	} // re: nups::predict


} // re: namespace nups

#endif



