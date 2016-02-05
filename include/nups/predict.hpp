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
			static void Predict(std::vector<NumT> & delta_x, std::vector<NumT> const& x, std::vector<NumT> const& rhs, typename NumTraits<NumT>::RealType const& delta_t)
			{	
				return PolyT::DoPredict(delta_x, x, rhs, delta_t);
			}
		};


		template<typename LinearSolverT>
		struct Euler : public PredictorBase<Euler<LinearSolverT> >
		{
			typedef LinearSolverT LinSolver;

			/**
			First-order Euler predictor, using a linear solver provided as a template type.
			*/
			template<typename NumT>
			static void DoPredict(std::vector<NumT> & delta_x, std::vector<NumT> const& x, std::vector<NumT> const& rhs, typename NumTraits<NumT>::RealType const& delta_t)
			{

				// the euler step
				LinearSolverT::Solve(delta_x, x, rhs);

				for (unsigned ii=0; ii<LinearSolverT::Degree; ++ii)
				{
					delta_x[ii] *= delta_t;
				}

			}
		};



		template<typename LinearSolverT>
		struct RK4 : public PredictorBase<RK4<LinearSolverT> >
		{
			typedef LinearSolverT LinSolver;

			/**
			Fourth-order Runge-Kutta predictor, using a linear solver provided as a template type.
			*/
			template<typename NumT>
			static void DoPredict(std::vector<NumT> & delta_x, std::vector<NumT> const& x, std::vector<NumT> const& rhs, typename NumTraits<NumT>::RealType const& delta_t)
			{
				std::vector<NumT> k1, k2, k3, k4;

				// get the first prediction, the euler step
				LinearSolverT::Solve(k1, x, rhs);


				// use the k1 to compute k2.  
				// use new space value x_2 = x+delta_t/2*k1
				std::vector<NumT> x_2 = x;
				for (unsigned ii=0; ii<LinearSolverT::Degree; ii++)
					x_2[ii] += delta_t/NumT(2)*k1[ii];

				LinearSolverT::Solve(k2, x_2, rhs);

				// use the k2 to compute k3.  
				// use new space value x_2 = x+delta_t/2*k2
				for (unsigned ii=0; ii<LinearSolverT::Degree; ii++)
					x_2[ii] = x[ii] + delta_t/NumT(2)*k2[ii];

				LinearSolverT::Solve(k3, x_2, rhs);

				// use the k3 to compute k4.  
				// use new space value x_2 = x+delta_t*k3
				for (unsigned ii=0; ii<LinearSolverT::Degree; ii++)
					x_2[ii] = x[ii] + delta_t*k3[ii];

				LinearSolverT::Solve(k4, x_2, rhs);


				// finally, combine the results of the above predictions, using a weighted sum
				delta_x.resize(LinearSolverT::Degree);
				for (unsigned ii=0; ii<LinearSolverT::Degree; ii++)
					delta_x[ii] = delta_t/NumT(6) * (k1[ii] + NumT(2)*k2[ii] + NumT(2)*k3[ii] + k4[ii]);

			}
		};

	} // re: nups::predict


} // re: namespace nups

#endif



