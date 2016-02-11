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

#include <cstdlib>
#include <complex>

namespace nups
{	
	template <typename T>
	struct NumTraits
	{
	};

	template<>
	struct NumTraits<float>
	{
		enum {IsComplex = 0};
		typedef float RealType;
		typedef std::complex<float> ComplexType;

		static float NewtonTerminationThreshold() 
		{
			return 1e-5f;
		}
	};

	template<>
	struct NumTraits<std::complex<float> >
	{
		enum {IsComplex = 1};
		typedef float RealType;
		typedef std::complex<float> ComplexType;

		static float NewtonTerminationThreshold() 
		{
			return 1e-5f;
		}
	};



	template<>
	struct NumTraits<double>
	{
		enum {IsComplex = 0};
		typedef double RealType;
		typedef std::complex<double> ComplexType;

		static double NewtonTerminationThreshold() 
		{
			return 1e-11;
		}
	};

	template<>
	struct NumTraits<std::complex<double> >
	{
		enum {IsComplex = 1};
		typedef double RealType;
		typedef std::complex<double> ComplexType;

		static double NewtonTerminationThreshold() 
		{
			return 1e-11;
		}
	};

	template<typename T>
	struct Random
	{
	};

	template<>
	struct Random<float>
	{
		static float Generate()
		{
			return 2*(float(rand())/RAND_MAX-0.5f);
		}
	};

	template<>
	struct Random<std::complex<float> >
	{
		static std::complex<float> Generate()
		{
			#ifdef HAVE_CPP_11
				std::default_random_engine generator;
				std::uniform_real_distribution<float> distribution(-1.0,1.0);
				std::complex<float> returnme(distribution(generator), distribution(generator));
			#else
				std::complex<float> returnme(2*(float(rand())/RAND_MAX-0.5f),2*(float(rand())/RAND_MAX-0.5f));
			#endif
			return returnme / float(sqrt(abs(returnme)));
		}
	};


	template<typename T>
	T NChooseK(unsigned n, unsigned k)
	{
		T result(1);

		for (unsigned ii = 1; ii <= k; ++ii)
		{
			result *= (n+1-ii);
		}

		for (unsigned ii = 1; ii <= k; ++ii)
		{
			result /= (ii);
		}
		return result;
	}

	template<>
	struct Random<double>
	{
		static double Generate()
		{
			return 2*(double(rand())/RAND_MAX-0.5);
		}
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
				std::complex<double> returnme(2*(double(rand())/RAND_MAX-0.5),2*(double(rand())/RAND_MAX-0.5));
			#endif
			return returnme / sqrt(abs(returnme));
		}
	};




	template<typename T>
		void print_to_screen_matlab(std::vector<std::complex<T> > const& v, std::string const& name)
		{	
			std::cout.precision(16);
			std::cout << name << " = [...\n";
			for (int ii = 0; ii < v.size(); ii++)
			{ // print kth coordinate
				std::cout << real(v[ii]) << "+1i*" << imag(v[ii])<<";\n";
			}
			std::cout << "];\n\n";
		}

		template<typename T>
		void print_to_screen_matlab(std::vector<T> const& v, std::string const& name)
		{	
			std::cout.precision(16);
			std::cout << name << " = [...\n";
			for (int ii = 0; ii < v.size(); ii++)
			{ // print kth coordinate
				std::cout << v[ii] << ";\n";
			}
			std::cout << "];\n\n";
		}

	

}


#endif

