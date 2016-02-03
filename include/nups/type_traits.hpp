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

// declares type traits for the nups library of numerical univariate polynomial solvers

#ifndef NUPS_TYPE_TRAITS_HPP
#define NUPS_TYPE_TRAITS_HPP

#include <complex>

namespace nups
{	
	template <typename T>
	struct TypeTraits
	{
	};

	template<>
	struct TypeTraits<double>
	{
		enum {IsComplex = 0};
		typedef double RealType;
		typedef std::complex<double> ComplexType;
	};

	template<>
	struct TypeTraits<std::complex<double> >
	{
		enum {IsComplex = 1};
		typedef double RealType;
		typedef std::complex<double> ComplexType;
	};

	template<typename T>
	struct Random
	{
	};

	template<>
	struct Random<std::complex<double> >
	{
		static std::complex<double> Generate()
		{
			#ifdef HAVE_CPP_11
				std::default_random_engine generator;
				std::uniform_real_distribution<double> distribution(-1.0,1.0);
				std::complex<double> returnme(distribution(generator), distribution(generator));
			#else
				std::complex<double> returnme(2*rand()/RAND_MAX-1,2*rand()/RAND_MAX-1);
			#endif
			return returnme / sqrt( abs(returnme));
		}
	};

}

#endif

