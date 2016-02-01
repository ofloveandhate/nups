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

namespace nups {

	namespace factor {

		template<int degree, typename PolyT>
		struct FactorizerBase
		{
			template<typename FT, typename ST>
			static void Factor(std::vector<FT> & coefficients_a, std::vector<FT> & coefficients_b, std::vector<ST> const& coefficients_start)
			{
				if (coefficients_start.size()!=degree && coefficients_start.size()!=degree+1)
				{
					std::stringstream error_message;
					error_message << "factoring a polynomial of degree " << degree << " must have " << degree << " or " << degree+1 << " coefficients.  yours has " << coefficients_start.size();
					throw std::runtime_error(error_message.str());
				}

				if (coefficients_start.size()==degree+1)
				{
					std::vector<ST> re_scaled_coefficients(degree);
					for (unsigned ii = 0; ii < degree; ++ii)
						re_scaled_coefficients[ii] = coefficients_start[ii] / coefficients_start[degree];

					return PolyT::DoFactorMonic(coefficients_a, coefficients_b, coefficients_start);
				}

				return PolyT::DoFactorMonic(coefficients_a, coefficients_b, coefficients_start);
			}
		};
		

		struct Octic : FactorizerBase<8,Octic>
		{

			// factors a monic octic univariate polynomial
			template<typename FT, typename ST>
			static void DoFactorMonic(std::vector<FT> & coefficients_a, std::vector<FT> & coefficients_b, std::vector<ST> const& coefficients_start)
			{
				
			}
		};


		struct Decic : FactorizerBase<10,Decic>
		{

			// factors a monic decic univariate polynomial
			template<typename FT, typename ST>
			static void DoFactorMonic(std::vector<FT> & coefficients_a, std::vector<FT> & coefficients_b, std::vector<ST> const& coefficients_start)
			{
			}
		};

	} // re: namespace factor

} // re: namespace nups

#endif
