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

// declares solvers

#ifndef NUPS_POLYNOMIAL_SOLVE_HPP
#define NUPS_POLYNOMIAL_SOLVE_HPP

#include <assert.h>
#include "nups/type_traits.hpp"
#include "nups/factor.hpp"
#include "nups/predict.hpp"
#include "nups/linear_solve.hpp"
#include <vector>

namespace nups {
	namespace solver {


		/**
		\brief The base solver class for NUPS.  Uses CRTP.

		\param degree The degree of the polynomial being solved.
		\param PolyT The type of polynomial being solved.  Seems to duplicate degree.

		\todo Eliminate degree as a template type.  Should be coupled with PolyT
		*/
		template <int degree, class PolyT>
		struct SolverBase
		{


			SolverBase(){}

			/**
			\brief Solve a quadratic equation

			\param[in] coefficients The input coefficients.  May be complex or real.
			\param[out] solutions The output, computed solutions.  Must be complex.

			You may put in a monic or non-monic polynomial to this function.
			*/
			template<typename CoeffT>
			void Solve(std::vector<typename NumTraits<CoeffT>::ComplexType>& solutions, std::vector<CoeffT> const& coefficients)
			{
				return SolveWithComplex(solutions, coefficients);
			}


			/**
			\brief Sharpens using Newton's method, interpreting the last entry in coeffs as the leading coefficient.

			\see SharpenMonic
			*/
			template<typename SolnT, typename CoeffT>
			static void SharpenNonMonic(std::vector<SolnT>& solutions, std::vector<CoeffT> const& coeffs, double final_tolerance, unsigned max_iterations)
			{
				
				for (typename std::vector<SolnT>::iterator iter=solutions.begin(); iter!=solutions.end(); ++iter)
				{	
					for (unsigned int ii(0); ii<max_iterations; ++ii)
					{
						SolnT delta_x = EvaluatePolyNonMonic(*iter,coeffs)/EvaluateDerivNonMonic(*iter,coeffs);
						*iter -= delta_x;
						if (abs(delta_x) < final_tolerance)
							break;
					}
				}
			}


			/**
			\brief Sharpens using Newton's method, assuming the leading coefficient is 1, and was omitted.

			\see SharpenNonMonic
			*/
			template<typename SolnT, typename CoeffT>
			static void SharpenMonic(std::vector<SolnT>& solutions, std::vector<CoeffT> const& coeffs, double final_tolerance, unsigned max_iterations)
			{
				for (typename std::vector<SolnT>::iterator iter=solutions.begin(); iter!=solutions.end(); ++iter)
					for (unsigned int ii(0); ii<max_iterations; ++ii)
					{
						SolnT delta_x = EvaluatePolyMonic(*iter,coeffs)/EvaluateDerivMonic(*iter,coeffs);
						*iter -= delta_x;
						if (abs(delta_x) < final_tolerance)
							break;
					}

			}

			/**
			\brief Evaluates a polynomial, interpreting the last entry in the vector of coefficients as the leading coefficient.
			
			\see EvaluatePolyMonic
			*/
			template<typename SolnT, typename CoeffT>
			static SolnT EvaluatePolyNonMonic(SolnT const& x, std::vector<CoeffT> const& coeffs)
			{
				SolnT x_val(1), f(0);
				for (typename std::vector<CoeffT>::const_iterator iter = coeffs.begin(); iter != coeffs.end(); ++iter)
				{
					f += *iter * x_val;
					x_val *= x;
				}
				return f;
			}

			/**
			\brief Evaluates a polynomial, assuming the leading coefficient is 1, and was omitted.

			assumes MONIC POLYNOMIAL, with omitted 1 at the front.

			\see EvaluatePolyNonMonic
			*/
			template<typename SolnT, typename CoeffT>
			static SolnT EvaluatePolyMonic(SolnT const& x, std::vector<CoeffT> const& coeffs)
			{	
				SolnT x_val(1), f(0);
				for (typename std::vector<CoeffT>::const_iterator iter = coeffs.begin(); iter != coeffs.end(); ++iter)
				{
					f += *iter * x_val;
					x_val *= x;
				}
				f += x_val;
				return f;
			}


			template<typename SolnT, typename CoeffT>
			static SolnT EvaluateDerivNonMonic(SolnT & x, std::vector<CoeffT> const& coeffs)
			{
				SolnT x_val(1), f(0);
				typename NumTraits<SolnT>::RealType pre_factor(1);

				for (typename std::vector<CoeffT>::const_iterator iter = coeffs.begin()+1; iter != coeffs.end(); ++iter)
				{
					f += *iter * x_val * pre_factor;
					pre_factor++;
					x_val *= x;
				}
				return f;
			}


			template<typename SolnT, typename CoeffT>
			static SolnT EvaluateDerivMonic(SolnT & x, std::vector<CoeffT> const& coeffs)
			{
				SolnT x_val(1), f(0);
				typename NumTraits<SolnT>::RealType pre_factor(1);

				for (typename std::vector<CoeffT>::const_iterator iter = coeffs.begin()+1; iter != coeffs.end(); ++iter)
				{
					f += *iter * x_val * pre_factor;
					pre_factor++;
					x_val *= x;
				}
				f += x_val*pre_factor;
				return f;
			}



		private:


			template<typename CoeffT>
			void SolveWithComplex(std::vector<typename NumTraits<CoeffT>::ComplexType>& solutions, std::vector<CoeffT> const& coefficients)
			{

				if (coefficients.size()!=degree+1 && coefficients.size()!=degree)
				{
					std::stringstream error_message;
					error_message << "solving a degree " << degree << " polynomial requires " << degree << " or " << degree+1 << " coefficients.";
					throw std::runtime_error(error_message.str());
				}


				//re-scale by the leading coefficient, which is the first coefficient by convention
				if (coefficients.size()==degree+1)
				{
					std::vector<CoeffT> re_scaled_coefficients(degree);
					for (unsigned ii = 0; ii < degree; ++ii)
						re_scaled_coefficients[ii] = coefficients[ii] / coefficients[degree];

					return static_cast< PolyT * >( this )->SolveWithComplexMonic(solutions, re_scaled_coefficients);
				}
				else
					return static_cast< PolyT * >( this )->SolveWithComplexMonic(solutions, coefficients);
			}

		};


		/**
		\brief Solver functor for solving degree two polynomials.
		*/
		struct Quadratic : public SolverBase<2, Quadratic>
		{

			template<typename CoeffT>
			void SolveWithComplexMonic(std::vector<typename NumTraits<CoeffT>::ComplexType>& solutions, std::vector<CoeffT> const& coefficients)
			{	
				typedef typename NumTraits<CoeffT>::ComplexType Real;
				typedef typename NumTraits<CoeffT>::ComplexType Complex;

				solutions.resize(2);

				Complex sqrt_discriminant = sqrt(pow(Complex(coefficients[1]),2) - Real(4)*Complex(coefficients[0]));

				solutions[1] = (-coefficients[1] + sqrt_discriminant)/Real(2);
				solutions[0] = (-coefficients[1] - sqrt_discriminant)/Real(2);
			}
		}; // Quadratic


		/**
		\brief Solver functor for solving degree four polynomials.
		*/
		struct Quartic : public SolverBase<4, Quartic>
		{
			template<typename CoeffT>
			void SolveWithComplexMonic(std::vector<typename NumTraits<CoeffT>::ComplexType>& solutions, std::vector<CoeffT> const& coefficients)
			{	
				typedef typename NumTraits<CoeffT>::ComplexType Real;
				typedef typename NumTraits<CoeffT>::ComplexType Complex;

				if (coefficients.size()!=4)
					throw std::runtime_error("solving a monic degree 4 polynomial requires 4 coefficients.");

				Real one_third = Real(1)/Real(3);

				solutions.resize(4);

				// the leading coefficient, a, is 1 by hypothesis
				const CoeffT& e = coefficients[0];
				const CoeffT& d = coefficients[1];
				const CoeffT& c = coefficients[2];
				const CoeffT& b = coefficients[3];

				CoeffT delta_0 = c*c - Real(3)*b*d + Real(12)*e;
				CoeffT delta_1 = Real(2)*c*c*c - Real(9)*b*c*d + Real(27)*b*b*e + Real(27)*d*d - Real(72)*c*e;

				CoeffT p = (Real(8)*c-Real(3)*b*b)/Real(8);
				CoeffT q = (b*b*b - Real(4)*b*c + Real(8)*d)/Real(8);

				CoeffT Q = pow((delta_1 + sqrt(delta_1*delta_1 - Real(4)*delta_0*delta_0*delta_0))/Real(2),one_third);
				CoeffT S = sqrt(Real(-2)/Real(3)*p + one_third*(Q+delta_0/Q))/Real(2);

				CoeffT y = sqrt(-Real(4)*S*S - Real(2)*p + q/S)/Real(2);

				CoeffT minus_b_over_four = -b/Real(4);

				solutions[0] = minus_b_over_four - S + y;
				solutions[1] = minus_b_over_four - S - y;

				y = sqrt(-Real(4)*S*S - Real(2)*p - q/S)/Real(2);

				solutions[2] = minus_b_over_four + S + y;
				solutions[3] = minus_b_over_four + S - y;
			}
			//  auto delta_0 = c*c - 3*b*d + 12*a*e;
			// 	auto delta_1 = 2*pow(c,3) - 9*b*c*d + 27*a*d*d - 72*a*c*e;

			// 	auto Q = pow((delta_1 + sqrt(delta_1*delta_1 - 4*pow(delta_0,3)))/2,one_third);
			// 	auto S = sqrt(minus_two_thirds*p + one_third/a*(Q+delta_0/Q))/2;

			// 	auto p = (8*a*c-3*b*b)/(8*a*a);
			// 	auto q = (pow(b,3) - 4*a*b*c + 8*a*a*d)/(8*a*a*a);

			// 	auto y = sqrt(-4*S*S - 2*P + q/S)/2;

			// 	auto min_b_over_four_a = -b/(4*a);

			// 	solutions[0] = min_b_over_four_a - S + y;
			// 	solutions[1] = min_b_over_four_a - S - y;

			// 	y = sqrt(-4*S*S - 2*P - q/S)/2;

			// 	solutions[2] = min_b_over_four_a + S + y;
			// 	solutions[3] = min_b_over_four_a + S - y;

		}; // Quartic



		/**
		Class for solving arbitrary univariate octic polynomials.  

		\tparam PredictorT The predictor method to be used during factorization.  Default is classic RK4.
		\tparam StartT The method for generating start points for the factorization.  Default is UnitCoefficient.
		*/
		template< template<class> class PredictorT = nups::predict::RK4,
				  template<int,int> class StartT = nups::factor::UnitCoefficient>
		struct Octic : public SolverBase<8, Octic< PredictorT, StartT > >
		{	
			typedef PredictorT<OcticLinear> Predictor;



			/**
			\brief Constructor for an octic solver

			\param sharpens The number of sharpen operations performed on the final solutions.
			\param steps The number of predict steps taken during factorization.
			\param corrects_during The number of Newton corrections taken after each step of the homotopy for the factorization.
			\param corrects_after The number of Newton corrections taken after factorization, but before solution of the factored quartics.
			*/
			Octic(double final_tolerance, unsigned max_iterations_final = 7, unsigned steps = 10, unsigned corrects_during = 7, unsigned corrects_after = 7) : max_iterations_(max_iterations_final), final_tolerance_(final_tolerance), factorizer_(steps, corrects_during, corrects_after)
			{
			}


			/**
			\brief Solve a monic octic.  

			Intended to be called from the base solver, but you can call it yourself if you need to
			*/
			template<typename CoeffT>
			void SolveWithComplexMonic(std::vector<typename NumTraits<CoeffT>::ComplexType>& solutions, std::vector<CoeffT> const& coefficients)
			{	

				typedef typename NumTraits<CoeffT>::ComplexType Real;
				typedef typename NumTraits<CoeffT>::ComplexType Complex;


				if (coefficients.size()!=8)
					throw std::runtime_error("solving a monic degree 8 polynomial requires 8 coefficients.");

				std::vector<Complex> factor_coeffs_1, factor_coeffs_2;
				factorizer_.Factor(factor_coeffs_1, factor_coeffs_2, coefficients);

				std::vector<Complex> solns_temp_1, solns_temp_2;
				quartic_solver_.Solve(solns_temp_1, factor_coeffs_1);
				quartic_solver_.Solve(solns_temp_2, factor_coeffs_2);


				solutions.resize(8);
				for (unsigned ii = 0; ii < 4; ++ii)
					solutions[ii]   = solns_temp_1[ii];

				for (unsigned ii = 0; ii < 4; ++ii)
					solutions[ii+4] = solns_temp_2[ii];

				this->SharpenMonic(solutions, coefficients, final_tolerance_, max_iterations_);
			}

			factor::Octic<Predictor,StartT> factorizer_;
			solver::Quartic quartic_solver_;
			double final_tolerance_;
			unsigned max_iterations_;
		}; // Octic


		template<template<class> class PredictorT = nups::predict::RK4,
				 template<int,int> class StartT = nups::factor::UnitCoefficient>
		struct Decic : public SolverBase<10, Decic<PredictorT,StartT> >
		{
			typedef PredictorT<DecicLinear> Predictor;

			template<typename CoeffT>
			void SolveWithComplexMonic(std::vector<typename NumTraits<CoeffT>::ComplexType>& solutions, std::vector<CoeffT> const& coefficients)
			{	
				typedef typename NumTraits<CoeffT>::ComplexType Real;
				typedef typename NumTraits<CoeffT>::ComplexType Complex;

				if (coefficients.size()!=10)
					throw std::runtime_error("solving a monic degree 10 polynomial requires 10 coefficients.");

				std::vector<Complex> factor_coeffs_1, factor_coeffs_2;
				factor::Decic<Predictor,StartT>::Factor(factor_coeffs_1, factor_coeffs_2, coefficients);


				std::vector<Complex> solns_temp_1, solns_temp_2;
				octic_solver_.Solve(solns_temp_1, factor_coeffs_1);
				quadratic_solver_.Solve(solns_temp_2, factor_coeffs_2);


				solutions.resize(10);
				for (unsigned ii = 0; ii < 8; ++ii)
					solutions[ii]   = solns_temp_1[ii];

				for (unsigned ii = 0; ii < 2; ++ii)
					solutions[ii+8] = solns_temp_2[ii];
			}

			solver::Octic<PredictorT> octic_solver_;
			solver::Quadratic quadratic_solver_;

		}; // Decic


	} // re: namespace nups::solver
} // re: namespace nups

#endif
