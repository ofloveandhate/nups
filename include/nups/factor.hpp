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

// declares factoring algorithms

#ifndef NUPS_FACTOR_HPP
#define NUPS_FACTOR_HPP
#include <iostream>
#include <vector>
#include <sstream>
#include "type_traits.hpp"

namespace nups {

	namespace factor {

		template<int degree, typename PolyT>
		struct FactorizerBase
		{
			template<typename NumT>
			static void Factor(std::vector<NumT> & coefficients_a, std::vector<NumT> & coefficients_b, std::vector<NumT> const& coefficients_start)
			{
				if (coefficients_start.size()!=degree && coefficients_start.size()!=degree+1)
				{
					std::stringstream error_message;
					error_message << "factoring a polynomial of degree " << degree << " must have " << degree << " or " << degree+1 << " coefficients.  yours has " << coefficients_start.size();
					throw std::runtime_error(error_message.str());
				}

				if (coefficients_start.size()==degree+1)
				{
					std::vector<NumT> re_scaled_coefficients(degree);
					for (unsigned ii = 0; ii < degree; ++ii)
						re_scaled_coefficients[ii] = coefficients_start[ii] / coefficients_start[degree];

					return PolyT::DoFactorMonic(coefficients_a, coefficients_b, coefficients_start);
				}

				return PolyT::DoFactorMonic(coefficients_a, coefficients_b, coefficients_start);
			}
		};
		
		template<typename PredictorT>
		struct Octic : public FactorizerBase<8,Octic<PredictorT> >
		{

			// factors a monic octic univariate polynomial into two quartics.
			template<typename NumT>
			static void DoFactorMonic(std::vector<NumT> & coefficients_a, std::vector<NumT> & coefficients_b, std::vector<NumT> const& a)
			{
				unsigned num_steps = 1000;


				// first, we make the start point
				std::vector<NumT> x(8);
				for (typename std::vector<NumT>::iterator iter=x.begin(); iter!=x.end(); iter++)
					*iter = nups::Random<NumT>::Generate();


				const NumT& r3 = x[0];
				const NumT& r2 = x[1];
				const NumT& r1 = x[2];
				const NumT& r0 = x[3];

				const NumT& s3 = x[4];
				const NumT& s2 = x[5];
				const NumT& s1 = x[6];
				const NumT& s0 = x[7];



				// now, we need to feed the start point into the polynomial to generate the random coefficients, a_star  (a^\ast).
				std::vector<NumT> a_star(8);
				//[r3+s3, r3*s3+r2+s2, r2*s3+r3*s2+r1+s1, r1*s3+r2*s2+r3*s1+r0+s0, r0*s3+r1*s2+r2*s1+r3*s0, r0*s2+r1*s1+r2*s0, r0*s1+r1*s0, r0*s0]
				a_star[0] = r3+s3;
				a_star[1] = r3*s3+r2+s2;
				a_star[2] = r2*s3+r3*s2+r1+s1;
				a_star[3] = r1*s3+r2*s2+r3*s1+r0+s0;
				a_star[4] = r0*s3+r1*s2+r2*s1+r3*s0;
				a_star[5] = r0*s2+r1*s1+r2*s0;
				a_star[6] = r0*s1+r1*s0;
				a_star[7] = r0*s0;

				
				
				std::vector<NumT> a_minus_a_star(8);
				for (unsigned ii=0; ii<8; ++ii)
					a_minus_a_star[ii] = a[ii] - a_star[ii];

				std::vector<NumT> a_star_minus_a(8);
				for (unsigned ii=0; ii<8; ++ii)
					a_star_minus_a[ii] = a_star[ii] - a[ii];

				typename TypeTraits<NumT>::RealType t(1);
				typename TypeTraits<NumT>::RealType delta_t(typename TypeTraits<NumT>::RealType(-1)/num_steps);
				std::vector<NumT> delta_x;

				
				print_to_screen_matlab(x,"x_start");
				print_to_screen_matlab(a,"a");
				print_to_screen_matlab(a_star,"a_star");
				std::cout << "delta_t: " << delta_t << std::endl;

				std::vector<NumT> residuals;
				for (unsigned n=0; n<num_steps; ++n)
				{
					if (0)
					{
						std::cout << "\n\ntaking step " << n << ", t: " << t << std::endl;
						print_to_screen_matlab(x,"x");
						
						Evaluate(residuals, x, a, a_star, t);
						print_to_screen_matlab(residuals,"residuals");
						
					}

					PredictorT::Predict(delta_x, x, a_star_minus_a, delta_t);

					if (0)
					{
						print_to_screen_matlab(delta_x,"delta_x");
					}

					for (unsigned ii=0; ii<8; ++ii)
						x[ii] += delta_x[ii];

					t += delta_t;

					for (unsigned kk=0; kk<2; kk++)
						Correct(x,a,a_star,t);
					
				}

				for (unsigned kk=0; kk<6; kk++)
					Correct(x,a,a_star,t);

				if (0)
					print_to_screen_matlab(x,"final_x");

				coefficients_a.resize(4);
				for (unsigned ii = 0; ii < 4; ++ii)
					coefficients_a[ii] = x[ii];

				coefficients_b.resize(4);
				for (unsigned ii = 0; ii < 4; ++ii)
					coefficients_b[ii] = x[ii+4];
			}

		private:

			template<typename NumT>
			static void Correct(std::vector<NumT> & x, std::vector<NumT> const& a, std::vector<NumT> const& a_star, typename TypeTraits<NumT>::RealType const& t)
			{
				std::vector<NumT> residuals;

				Evaluate(residuals, x, a, a_star, t);

				for (unsigned ii=0; ii<8; ++ii)
					residuals[ii] = -residuals[ii];

				std::vector<NumT> delta_x;
				PredictorT::LinSolver::Solve(delta_x, x, residuals);
				for (unsigned ii=0; ii<8; ++ii)
					x[ii] += delta_x[ii];
			}

			template<typename NumT>
			static void Evaluate(std::vector<NumT> & f, std::vector<NumT> const& rs, std::vector<NumT> const& a, std::vector<NumT> const& a_star, typename TypeTraits<NumT>::RealType const& t)
			{	
				const NumT& r3 = rs[0];
				const NumT& r2 = rs[1];
				const NumT& r1 = rs[2];
				const NumT& r0 = rs[3];

				const NumT& s3 = rs[4];
				const NumT& s2 = rs[5];
				const NumT& s1 = rs[6];
				const NumT& s0 = rs[7];

				f.resize(8);
				f[0] = r3+s3 					- (t*a_star[0] + (NumT(1)-t)*a[0]);
				f[1] = r3*s3+r2+s2 				- (t*a_star[1] + (NumT(1)-t)*a[1]);
				f[2] = r2*s3+r3*s2+r1+s1 		- (t*a_star[2] + (NumT(1)-t)*a[2]);
				f[3] = r1*s3+r2*s2+r3*s1+r0+s0 	- (t*a_star[3] + (NumT(1)-t)*a[3]);
				f[4] = r0*s3+r1*s2+r2*s1+r3*s0 	- (t*a_star[4] + (NumT(1)-t)*a[4]);
				f[5] = r0*s2+r1*s1+r2*s0 		- (t*a_star[5] + (NumT(1)-t)*a[5]);
				f[6] = r0*s1+r1*s0 				- (t*a_star[6] + (NumT(1)-t)*a[6]);
				f[7] = r0*s0 					- (t*a_star[7] + (NumT(1)-t)*a[7]);
			}
		};

		//[ r0*s0, r0*s1 + r1*s0, r0*s2 + r1*s1 + r2*s0, r0*s3 + r1*s2 + r2*s1 + r3*s0, r0 + s0 + r1*s3 + r2*s2 + r3*s1, r1 + s1 + r2*s3 + r3*s2, r2 + s2 + r3*s3, r3 + s3, 1]
		struct Decic : public FactorizerBase<10,Decic>
		{

			// factors a monic decic univariate polynomial
			template<typename NumT>
			static void DoFactorMonic(std::vector<NumT> & coefficients_a, std::vector<NumT> & coefficients_b, std::vector<NumT> const& coefficients_start)
			{
			}
		};

	} // re: namespace factor

} // re: namespace nups

#endif
