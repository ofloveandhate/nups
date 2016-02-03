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
			static void DoFactorMonic(std::vector<NumT> & coefficients_a, std::vector<NumT> & coefficients_b, std::vector<NumT> const& coefficients_start)
			{
				unsigned num_steps = 4;

				// first, we make the start point
				std::vector<NumT> x(8);
				for (typename std::vector<NumT>::iterator iter=x.begin(); iter!=x.begin(); iter++)
					*iter = nups::Random<NumT>::Generate();

				// now, we need to feed the start point into the polynomial to generate the random coefficients, a_star  (a^\ast).
				//
				//

				const NumT& r3 = x[0];
				const NumT& r2 = x[1];
				const NumT& r1 = x[2];
				const NumT& r0 = x[3];

				const NumT& s3 = x[4];
				const NumT& s2 = x[5];
				const NumT& s1 = x[6];
				const NumT& s0 = x[7];

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

				NumT t(1);
				NumT delta_t(typename TypeTraits<NumT>::RealType(1)/num_steps);
				std::vector<NumT> delta_x;

				std::vector<NumT> rhs(8);
				for (unsigned ii=0; ii<8; ++ii)
					rhs[ii] = a_star[ii]-coefficients_start[ii];

				for (unsigned n=0; n<num_steps; ++n)
				{
					PredictorT::Predict(delta_x, x, rhs, t, delta_t);
				}
			}

		};


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
