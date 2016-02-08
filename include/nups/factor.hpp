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
			static void Factor(std::vector<typename NumTraits<NumT>::ComplexType> & r, std::vector<typename NumTraits<NumT>::ComplexType> & s, std::vector<NumT> const& p)
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
			static void DoFactorMonic(std::vector<typename NumTraits<NumT>::ComplexType> & r, std::vector<typename NumTraits<NumT>::ComplexType> & s, std::vector<NumT> const& a)
			{

				unsigned num_steps = 20;
				unsigned num_corrects_during = 5;
				unsigned num_corrects_after = 5;



				typedef typename NumTraits<NumT>::RealType Real;
				typedef typename NumTraits<NumT>::ComplexType Complex;

				Real t(1);
				Real delta_t(Real(-1)/num_steps);

				std::vector<Complex> delta_x(PolyT::Degree);
				std::vector<Complex> x(PolyT::Degree);
				std::vector<Complex> a_star(PolyT::Degree);
				std::vector<Complex> a_star_minus_a(PolyT::Degree);
				std::vector<Complex> residuals(PolyT::Degree);



				// make the start point
				GenerateStartPoint(x);

				// now, we need to feed the start point into the polynomial to generate the random coefficients, a_star  (a^\ast).
				PolyT::EvaluateF(a_star, x);
				
				// we use the vector a^\ast - a as the right hand side during linear solves, so we cache the vector.
				A_Star_Minus_A(a_star_minus_a, a_star, a);

				
				
				for (unsigned n=0; n<num_steps; ++n)
				{

					PolyT::Predictor::Predict(delta_x, x, a_star_minus_a, delta_t);

					for (unsigned ii=0; ii<PolyT::Degree; ++ii)
						x[ii] += delta_x[ii];

					t += delta_t;

					for (unsigned kk=0; kk<num_corrects_during; kk++)
						Correct(x,a,a_star,t);
				}

				for (unsigned kk=0; kk<num_corrects_after; kk++)
					Correct(x,a);

				// here we reverse the order of the coefficients back into `subscripts matching`.
				r.resize(PolyT::DegreeFactorR);
				for (unsigned ii = 0; ii < PolyT::DegreeFactorR; ++ii)
					r[PolyT::DegreeFactorR-1-ii] = x[ii];

				s.resize(PolyT::DegreeFactorS);
				for (unsigned ii = 0; ii < PolyT::DegreeFactorS; ++ii)
					s[PolyT::DegreeFactorS-1-ii] = x[ii+PolyT::DegreeFactorR];
			}


			template<typename NumT, bool UseBinomial = false>
			static void GenerateStartPoint(std::vector<NumT> & rs_start)
			{
				if (UseBinomial)
				{
					unsigned counter(0);
					for (unsigned ii(0); ii<PolyT::DegreeFactorR; ++ii)
						rs_start[ii] = nups::Random<NumT>::Generate()*NChooseK<typename NumTraits<NumT>::RealType>(PolyT::DegreeFactorR,counter++);

					counter = 0;
					for (unsigned ii(0); ii<PolyT::DegreeFactorS; ++ii)
						rs_start[ii+PolyT::DegreeFactorR] = nups::Random<NumT>::Generate()*NChooseK<typename NumTraits<NumT>::RealType>(PolyT::DegreeFactorS,counter++);
				}
				else
				{
					for (unsigned ii(0); ii<PolyT::DegreeFactorR; ++ii)
						rs_start[ii] = nups::Random<NumT>::Generate();

					for (unsigned ii(0); ii<PolyT::DegreeFactorS; ++ii)
						rs_start[ii+PolyT::DegreeFactorR] = nups::Random<NumT>::Generate();
				}
			}

			template<typename AStarNumT, typename ANumT>
			static void A_Star_Minus_A(std::vector<AStarNumT> & a_star_minus_a, std::vector<AStarNumT> & a_star, std::vector<ANumT> const& a)
			{
				for (unsigned ii=0; ii<PolyT::Degree; ++ii)
					a_star_minus_a[ii] = a_star[ii] - a[ii];
			}


			template<typename ANumT, typename AStarNumT>
			static void Correct(std::vector<AStarNumT> & x, std::vector<ANumT> const& a, std::vector<AStarNumT> const& a_star, typename NumTraits<ANumT>::RealType const& t)
			{
				std::vector<AStarNumT> residuals;

				PolyT::EvaluateHomotopy(residuals, x, a, a_star, t);

				for (unsigned ii=0; ii<PolyT::Degree; ++ii)
					residuals[ii] = -residuals[ii];

				std::vector<AStarNumT> delta_x;
				PolyT::Predictor::LinSolver::Solve(delta_x, x, residuals);
				for (unsigned ii=0; ii<PolyT::Degree; ++ii)
					x[ii] += delta_x[ii];
			}

			template<typename ANumT, typename XNumT>
			static void Correct(std::vector<XNumT> & x, std::vector<ANumT> const& a)
			{
				std::vector<XNumT> residuals;

				PolyT::EvaluateHomotopy(residuals, x, a);

				for (unsigned ii=0; ii<PolyT::Degree; ++ii)
					residuals[ii] = -residuals[ii];

				std::vector<XNumT> delta_x;
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
			static void EvaluateF(std::vector<NumT> & f, std::vector<NumT> const& rs)
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
				f[0] = r3+s3;
				f[1] = r3*s3+r2+s2;
				f[2] = r2*s3+r3*s2+r1+s1;
				f[3] = r1*s3+r2*s2+r3*s1+r0+s0;
				f[4] = r0*s3+r1*s2+r2*s1+r3*s0;
				f[5] = r0*s2+r1*s1+r2*s0;
				f[6] = r0*s1+r1*s0;
				f[7] = r0*s0;
			}



			


			template<typename AStarNumT, typename ANumT>
			static void EvaluateHomotopy(std::vector<AStarNumT> & f, std::vector<AStarNumT> const& rs, std::vector<ANumT> const& a, std::vector<AStarNumT> const& a_star, typename NumTraits<ANumT>::RealType const& t)
			{	

				typedef typename NumTraits<ANumT>::RealType Real;
				const AStarNumT& r3 = rs[0];
				const AStarNumT& r2 = rs[1];
				const AStarNumT& r1 = rs[2];
				const AStarNumT& r0 = rs[3];

				const AStarNumT& s3 = rs[4];
				const AStarNumT& s2 = rs[5];
				const AStarNumT& s1 = rs[6];
				const AStarNumT& s0 = rs[7];

				f.resize(8);
				f[0] = r3+s3 					- (t*a_star[0] + (Real(1)-t)*a[0]);
				f[1] = r3*s3+r2+s2 				- (t*a_star[1] + (Real(1)-t)*a[1]);
				f[2] = r2*s3+r3*s2+r1+s1 		- (t*a_star[2] + (Real(1)-t)*a[2]);
				f[3] = r1*s3+r2*s2+r3*s1+r0+s0 	- (t*a_star[3] + (Real(1)-t)*a[3]);
				f[4] = r0*s3+r1*s2+r2*s1+r3*s0 	- (t*a_star[4] + (Real(1)-t)*a[4]);
				f[5] = r0*s2+r1*s1+r2*s0 		- (t*a_star[5] + (Real(1)-t)*a[5]);
				f[6] = r0*s1+r1*s0 				- (t*a_star[6] + (Real(1)-t)*a[6]);
				f[7] = r0*s0 					- (t*a_star[7] + (Real(1)-t)*a[7]);
			}

			template<typename XNumT, typename ANumT>
			static void EvaluateHomotopy(std::vector<XNumT> & f, std::vector<XNumT> const& rs, std::vector<ANumT> const& a)
			{	
				const XNumT& r3 = rs[0];
				const XNumT& r2 = rs[1];
				const XNumT& r1 = rs[2];
				const XNumT& r0 = rs[3];

				const XNumT& s3 = rs[4];
				const XNumT& s2 = rs[5];
				const XNumT& s1 = rs[6];
				const XNumT& s0 = rs[7];

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
