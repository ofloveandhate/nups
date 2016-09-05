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
			
			template<typename NumT, 
								template <typename, typename...> class ContT, 
			          			typename... Alloc>
			static void Predict(ContT<NumT, Alloc...> & delta_x, ContT<NumT, Alloc...> const& x, ContT<NumT, Alloc...> const& rhs, typename NumTraits<NumT>::RealType const& delta_t)
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
			template<typename NumT, 
								template <typename, typename...> class ContT, 
			          			typename... Alloc>
			static void DoPredict(ContT<NumT, Alloc...> & delta_x, ContT<NumT, Alloc...> const& x, ContT<NumT, Alloc...> const& rhs, typename NumTraits<NumT>::RealType const& delta_t)
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
			template<typename NumT, 
								template <typename, typename...> class ContT, 
			          			typename... Alloc>
			static void DoPredict(ContT<NumT, Alloc...> & delta_x, ContT<NumT, Alloc...> const& x, ContT<NumT, Alloc...> const& rhs, typename NumTraits<NumT>::RealType const& delta_t)
			{
				ContT<NumT, Alloc...> k1, k2, k3, k4;

				// get the first prediction, the euler step
				LinearSolverT::Solve(k1, x, rhs);


				// use the k1 to compute k2.  
				// use new space value x_2 = x+delta_t/2*k1
				ContT<NumT, Alloc...> x_2 = x;
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




		template<typename LinearSolverT>
		struct RKKC45 : public PredictorBase<RKKC45<LinearSolverT> >
		{
			typedef LinearSolverT LinSolver;

			/**
			Fifth-order Cash-Karp predictor, using a linear solver provided as a template type.
			*/
			template<typename NumT, 
								template <typename, typename...> class ContT, 
			          			typename... Alloc>
			static void DoPredict(ContT<NumT, Alloc...> & delta_x, ContT<NumT, Alloc...> const& x, ContT<NumT, Alloc...> const& rhs, typename NumTraits<NumT>::RealType const& delta_t)
			{
				ContT<NumT, Alloc...> k1, k2, k3, k4, k5, k6;

				// get the first prediction, the euler step
				LinearSolverT::Solve(k1, x, rhs);


				// use the k1 to compute k2.  
				ContT<NumT, Alloc...> x_2 = x;
				for (unsigned ii=0; ii<LinearSolverT::Degree; ii++)
					x_2[ii] += delta_t/NumT(5)*k1[ii];
				LinearSolverT::Solve(k2, x_2, rhs);

				// use the k2 to compute k3.  
				for (unsigned ii=0; ii<LinearSolverT::Degree; ii++)
					x_2[ii] = x[ii] + delta_t*(NumT(3)/NumT(40)*k1[ii] + NumT(9)/NumT(40)*k2[ii]);
				LinearSolverT::Solve(k3, x_2, rhs);

				// use the k3 to compute k4.  
				for (unsigned ii=0; ii<LinearSolverT::Degree; ii++)
					x_2[ii] = x[ii] + delta_t*(NumT(3)/NumT(10)*k1[ii] + NumT(-9)/NumT(10)*k2[ii] + NumT(6)/NumT(5)*k3[ii]);
				LinearSolverT::Solve(k4, x_2, rhs);

				// use the k4 to compute k5.  
				for (unsigned ii=0; ii<LinearSolverT::Degree; ii++)
					x_2[ii] = x[ii] + delta_t*(NumT(-11)/NumT(54)*k1[ii] + NumT(5)/NumT(2)*k2[ii] + NumT(-70)/NumT(27)*k3[ii] + NumT(35)/NumT(27)*k4[ii]);
				LinearSolverT::Solve(k5, x_2, rhs);

				// use the k5 to compute k6.  
				for (unsigned ii=0; ii<LinearSolverT::Degree; ii++)
					x_2[ii] = x[ii] + delta_t*(NumT(1631)/NumT(55296)*k1[ii] + NumT(175)/NumT(512)*k2[ii] + NumT(575)/NumT(13824)*k3[ii] + NumT(44275)/NumT(110592)*k4[ii] + NumT(253)/NumT(4096)*k5[ii]);
				LinearSolverT::Solve(k6, x_2, rhs);

				// finally, combine the results of the above predictions, using a weighted sum
				delta_x.resize(LinearSolverT::Degree);
				for (unsigned ii=0; ii<LinearSolverT::Degree; ii++)
					delta_x[ii] = delta_t * (NumT(37)/NumT(378)*k1[ii] + /* 0 */ NumT(250)/NumT(621)*k3[ii] + NumT(125)/NumT(594)*k4[ii] + /* 0 */ NumT(512)/NumT(1771)*k6[ii]);

			}
		};
	} // re: nups::predict


} // re: namespace nups

#endif



