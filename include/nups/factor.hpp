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

		template<typename PolyT>
		struct FactorizerBase
		{
			/**
			Factors the polynomial represented by p into two polynomials, r and s.
			*/
			template<typename NumT>
			static void Factor(std::vector<typename TypeTraits<NumT>::ComplexType> & r, std::vector<typename TypeTraits<NumT>::ComplexType> & s, std::vector<NumT> const& p)
			{
				if (p.size()!=PolyT::Degree && p.size()!=PolyT::Degree+1)
				{
					std::stringstream error_message;
					error_message << "factoring a polynomial of PolyT::Degree " << PolyT::Degree << " must have " << PolyT::Degree << " or " << PolyT::Degree+1 << " coefficients.  yours has " << p.size();
					throw std::runtime_error(error_message.str());
				}

				std::vector<NumT> re_ordered_coefficients(PolyT::Degree);

				if (p.size()==PolyT::Degree+1)
				{
					for (unsigned ii = 0; ii < PolyT::Degree; ++ii)
						re_ordered_coefficients[PolyT::Degree-1 - ii] = p[ii] / p[PolyT::Degree];
				}
				else
				{
					for (unsigned ii = 0; ii < PolyT::Degree; ++ii)
						re_ordered_coefficients[PolyT::Degree-1 - ii] = p[ii];
				}

				return PolyT::DoFactorMonic(r, s, re_ordered_coefficients);
			}

		private:
			// factors a monic octic univariate polynomial into two quartics.
			template<typename NumT>
			static void DoFactorMonic(std::vector<typename TypeTraits<NumT>::ComplexType> & r, std::vector<typename TypeTraits<NumT>::ComplexType> & s, std::vector<NumT> const& a)
			{

				unsigned num_steps = 25;
				unsigned num_corrects_during = 5;
				unsigned num_corrects_after = 6;





				typename TypeTraits<NumT>::RealType t(1);
				typename TypeTraits<NumT>::RealType delta_t(typename TypeTraits<NumT>::RealType(-1)/num_steps);
				std::vector<NumT> delta_x;

				std::vector<NumT> x(PolyT::Degree);
				std::vector<NumT> a_star(PolyT::Degree);
				std::vector<NumT> a_star_minus_a(PolyT::Degree);
				std::vector<NumT> residuals(PolyT::Degree);



				// make the start point
				GenerateStartPoint(x);

				// now, we need to feed the start point into the polynomial to generate the random coefficients, a_star  (a^\ast).
				PolyT::EvaluateF(a_star, x);
				
				// we use the vector a^\ast - a as the right hand side during linear solves, so we cache the vector.
				ComputeA_Star_Minus_A(a_star_minus_a, a_star, a);

				

				

				if (0)
				{
					print_to_screen_matlab(x,"x_start");
					print_to_screen_matlab(a,"a");
					print_to_screen_matlab(a_star,"a_star");
					std::cout << "delta_t: " << delta_t << std::endl;
				}
				
				for (unsigned n=0; n<num_steps; ++n)
				{
					if (0)
					{
						std::cout << "\n\ntaking step " << n << ", t: " << t << std::endl;
						print_to_screen_matlab(x,"x");
						
						PolyT::EvaluateHomotopy(residuals, x, a, a_star, t);
						print_to_screen_matlab(residuals,"residuals");
						
					}

					PolyT::Predictor::Predict(delta_x, x, a_star_minus_a, delta_t);

					if (0)
					{
						print_to_screen_matlab(delta_x,"delta_x");
					}

					for (unsigned ii=0; ii<PolyT::Degree; ++ii)
						x[ii] += delta_x[ii];

					t += delta_t;

					for (unsigned kk=0; kk<num_corrects_during; kk++)
						Correct(x,a,a_star,t);
					
				}

				for (unsigned kk=0; kk<num_corrects_after; kk++)
					Correct(x,a);

				if (0)
					print_to_screen_matlab(x,"final_x");


				// here we reverse the order of the coefficients back into `subscripts matching`.
				r.resize(PolyT::DegreeFactorR);
				for (unsigned ii = 0; ii < PolyT::DegreeFactorR; ++ii)
					r[PolyT::DegreeFactorR-1-ii] = x[ii];

				s.resize(PolyT::DegreeFactorS);
				for (unsigned ii = 0; ii < PolyT::DegreeFactorS; ++ii)
					s[PolyT::DegreeFactorS-1-ii] = x[ii+PolyT::DegreeFactorR];
			}


			template<typename NumT>
			static void GenerateStartPoint(std::vector<NumT> & rs_start)
			{
				for (typename std::vector<NumT>::iterator iter=rs_start.begin(); iter!=rs_start.end(); iter++)
					*iter = nups::Random<NumT>::Generate();
			}

			template<typename NumT>
			static void ComputeA_Star_Minus_A(std::vector<NumT> & a_star_minus_a, std::vector<NumT> & a_star, std::vector<NumT> const& a)
			{
				for (unsigned ii=0; ii<PolyT::Degree; ++ii)
					a_star_minus_a[ii] = a_star[ii] - a[ii];
			}


			template<typename NumT>
			static void Correct(std::vector<NumT> & x, std::vector<NumT> const& a, std::vector<NumT> const& a_star, typename TypeTraits<NumT>::RealType const& t)
			{
				std::vector<NumT> residuals;

				PolyT::EvaluateHomotopy(residuals, x, a, a_star, t);

				for (unsigned ii=0; ii<PolyT::Degree; ++ii)
					residuals[ii] = -residuals[ii];

				std::vector<NumT> delta_x;
				PolyT::Predictor::LinSolver::Solve(delta_x, x, residuals);
				for (unsigned ii=0; ii<PolyT::Degree; ++ii)
					x[ii] += delta_x[ii];
			}

			template<typename NumT>
			static void Correct(std::vector<NumT> & x, std::vector<NumT> const& a)
			{
				std::vector<NumT> residuals;

				PolyT::EvaluateHomotopy(residuals, x, a);

				for (unsigned ii=0; ii<PolyT::Degree; ++ii)
					residuals[ii] = -residuals[ii];

				std::vector<NumT> delta_x;
				PolyT::Predictor::LinSolver::Solve(delta_x, x, residuals);
				for (unsigned ii=0; ii<PolyT::Degree; ++ii)
					x[ii] += delta_x[ii];
			}
		}; // re: FactorizerBase
		









		template<typename PredictorT>
		struct Octic : public FactorizerBase<Octic<PredictorT> >
		{
			enum {
				  Degree = 8,
				  DegreeFactorR = 4,
				  DegreeFactorS = 4
				 };

			typedef PredictorT Predictor;


			


			template<typename NumT>
			static void EvaluateF(std::vector<NumT> & a_star, std::vector<NumT> const& rs)
			{
				const NumT& r3 = rs[0];
				const NumT& r2 = rs[1];
				const NumT& r1 = rs[2];
				const NumT& r0 = rs[3];

				const NumT& s3 = rs[4];
				const NumT& s2 = rs[5];
				const NumT& s1 = rs[6];
				const NumT& s0 = rs[7];
				
				//[r3+s3, r3*s3+r2+s2, r2*s3+r3*s2+r1+s1, r1*s3+r2*s2+r3*s1+r0+s0, r0*s3+r1*s2+r2*s1+r3*s0, r0*s2+r1*s1+r2*s0, r0*s1+r1*s0, r0*s0]
				a_star[0] = r3+s3;
				a_star[1] = r3*s3+r2+s2;
				a_star[2] = r2*s3+r3*s2+r1+s1;
				a_star[3] = r1*s3+r2*s2+r3*s1+r0+s0;
				a_star[4] = r0*s3+r1*s2+r2*s1+r3*s0;
				a_star[5] = r0*s2+r1*s1+r2*s0;
				a_star[6] = r0*s1+r1*s0;
				a_star[7] = r0*s0;
			}



			


			template<typename NumT>
			static void EvaluateHomotopy(std::vector<NumT> & f, std::vector<NumT> const& rs, std::vector<NumT> const& a, std::vector<NumT> const& a_star, typename TypeTraits<NumT>::RealType const& t)
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

			template<typename NumT>
			static void EvaluateHomotopy(std::vector<NumT> & f, std::vector<NumT> const& rs, std::vector<NumT> const& a)
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
				f[0] = r3+s3 					- a[0];
				f[1] = r3*s3+r2+s2 				- a[1];
				f[2] = r2*s3+r3*s2+r1+s1 		- a[2];
				f[3] = r1*s3+r2*s2+r3*s1+r0+s0 	- a[3];
				f[4] = r0*s3+r1*s2+r2*s1+r3*s0 	- a[4];
				f[5] = r0*s2+r1*s1+r2*s0 		- a[5];
				f[6] = r0*s1+r1*s0 				- a[6];
				f[7] = r0*s0 					- a[7];
			}
		};




		//[ r0*s0, r0*s1 + r1*s0, r0*s2 + r1*s1 + r2*s0, r0*s3 + r1*s2 + r2*s1 + r3*s0, r0 + s0 + r1*s3 + r2*s2 + r3*s1, r1 + s1 + r2*s3 + r3*s2, r2 + s2 + r3*s3, r3 + s3, 1]
		struct Decic : public FactorizerBase<Decic>
		{
			enum {
				  Degree = 10,
				  DegreeFactorR = 8,
				  DegreeFactorS = 2
				 };

		};

	} // re: namespace factor

} // re: namespace nups

#endif
