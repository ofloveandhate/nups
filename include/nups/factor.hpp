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
#include <stdexcept>

namespace nups {

	namespace factor {


		template<typename T>
		struct StartPointGenerator
		{
		};


		template<int DegreeFactorR, int DegreeFactorS>
		struct UnitCoefficient : StartPointGenerator<UnitCoefficient<DegreeFactorR, DegreeFactorS> >
		{
			/**
			\brief Generate a start point for the homotopy factoring a desired polynomial.

			\param[out] rs_start The generated start point
			*/
			template<typename NumT>
			static void Generate(std::vector<NumT> & rs_start)
			{
				rs_start.resize(DegreeFactorR+DegreeFactorS);
				for (unsigned ii(0); ii<DegreeFactorR+DegreeFactorS; ++ii)
					rs_start[ii] = nups::Random<NumT>::Generate();
			}
		};


		template<int DegreeFactorR, int DegreeFactorS>
		struct BinomialCoefficient : StartPointGenerator<BinomialCoefficient<DegreeFactorR, DegreeFactorS> >
		{
			/**
			\brief Generate a start point for the homotopy factoring a desired polynomial.

			\param[out] rs_start The generated start point
			*/
			template<typename NumT>
			static void Generate(std::vector<NumT> & rs_start)
			{
				rs_start.resize(DegreeFactorR+DegreeFactorS);
				unsigned counter(0);
				for (unsigned ii(0); ii<DegreeFactorR; ++ii)
					rs_start[ii] = nups::Random<NumT>::Generate()*NChooseK<typename NumTraits<NumT>::RealType>(DegreeFactorR,counter++);

				counter = 0;
				for (unsigned ii(0); ii<DegreeFactorS; ++ii)
					rs_start[ii+DegreeFactorR] = nups::Random<NumT>::Generate()*NChooseK<typename NumTraits<NumT>::RealType>(DegreeFactorS,counter++);

			}
		};









		/**
		\brief Base class for factoring algorithms
		*/
		template<typename PolyT, typename StartT>
		struct FactorizerBase
		{

			FactorizerBase(unsigned num_steps, unsigned num_corrects_during, unsigned num_corrects_after) : num_steps_(num_steps), num_corrects_during_(num_corrects_during), num_corrects_after_(num_corrects_after)
			{
			}


			/**
			Factors the polynomial represented by p into two polynomials, r and s.

			\param[out] r The higher degree of the two output monics.
			\param[out] s The lower degree of the two output monics.
			\param[int] p The input monic polynomial, which you want to factor.
			*/
			template<typename NumT>
			void Factor(std::vector<typename NumTraits<NumT>::ComplexType> & r, std::vector<typename NumTraits<NumT>::ComplexType> & s, std::vector<NumT> const& p)
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
				// push off the factoring to the private function which assumes monic.
				return static_cast< PolyT * >( this )->DoFactorMonic(r, s, re_ordered_coefficients);
			}

		private:

			unsigned num_steps_;
			unsigned num_corrects_during_;
			unsigned num_corrects_after_;

			/**
			\brief Factors a monic univariate polynomial into two lower-degree monic polynomials.
			*/
			template<typename NumT>
			void DoFactorMonic(std::vector<typename NumTraits<NumT>::ComplexType> & r, std::vector<typename NumTraits<NumT>::ComplexType> & s, std::vector<NumT> const& a)
			{
				typedef typename NumTraits<NumT>::RealType Real;
				typedef typename NumTraits<NumT>::ComplexType Complex;

				Real t(1);
				Real delta_t(Real(-1)/num_steps_);

				std::vector<Complex> delta_x(PolyT::Degree);
				std::vector<Complex> x(PolyT::Degree);
				std::vector<Complex> a_star(PolyT::Degree);
				std::vector<Complex> a_star_minus_a(PolyT::Degree);
				std::vector<Complex> residuals(PolyT::Degree);



				// make the start point
				StartT::Generate(x);

				// now, we need to feed the start point into the polynomial to generate the random coefficients, a_star  (a^\ast).
				PolyT::EvaluateF(a_star, x);
				
				// we use the vector a^\ast - a as the right hand side during linear solves, so we cache the vector.
				A_Star_Minus_A(a_star_minus_a, a_star, a);

				
				
				for (unsigned n=0; n<num_steps_; ++n)
				{

					PolyT::Predictor::Predict(delta_x, x, a_star_minus_a, delta_t);

					for (unsigned ii=0; ii<PolyT::Degree; ++ii)
						x[ii] += delta_x[ii];

					t += delta_t;

					for (unsigned kk=0; kk<num_corrects_during_; kk++)
						if (Correct(x,a,a_star,t))
							break;
				}

				for (unsigned kk=0; kk<num_corrects_after_; kk++)
					if (Correct(x,a))
						break;

				// here we reverse the order of the coefficients back into `subscripts matching`.
				r.resize(PolyT::DegreeFactorR);
				for (unsigned ii = 0; ii < PolyT::DegreeFactorR; ++ii)
					r[PolyT::DegreeFactorR-1-ii] = x[ii];

				s.resize(PolyT::DegreeFactorS);
				for (unsigned ii = 0; ii < PolyT::DegreeFactorS; ++ii)
					s[PolyT::DegreeFactorS-1-ii] = x[ii+PolyT::DegreeFactorR];
			}


			template<typename NumT>
			static void A_Star_Minus_A(std::vector<typename NumTraits<NumT>::ComplexType> & a_star_minus_a, std::vector<typename NumTraits<NumT>::ComplexType> & a_star, std::vector<NumT> const& a)
			{	
				a_star_minus_a.resize(PolyT::Degree);
				for (unsigned ii=0; ii<PolyT::Degree; ++ii)
					a_star_minus_a[ii] = a_star[ii] - a[ii];
			}

			/**
			\return true if correction loop should terminate, false otherwise.
			*/
			template<typename NumT>
			static bool Correct(std::vector<typename NumTraits<NumT>::ComplexType> & x, 
			                    std::vector<NumT> const& a, 
			                    std::vector<typename NumTraits<NumT>::ComplexType> const& a_star, 
			                    typename NumTraits<NumT>::RealType const& t)
			{
				std::vector<typename NumTraits<NumT>::ComplexType> residuals(PolyT::Degree);

				PolyT::EvaluateHomotopy(residuals, x, a, a_star, t);


				std::vector<typename NumTraits<NumT>::ComplexType> delta_x;
				PolyT::Predictor::LinSolver::Solve(delta_x, x, residuals);

				for (unsigned ii=0; ii<PolyT::Degree; ++ii)
					x[ii] -= delta_x[ii];

				for (typename std::vector<typename NumTraits<NumT>::ComplexType>::const_iterator iter = delta_x.begin(); iter!=delta_x.end(); ++iter)
				{
					if (abs(*iter) > NumTraits<NumT>::NewtonTerminationThreshold())
						return false;
				}

				return true;
				// if (max_norm < 1e-15)
					// std::cout << "reached numerical zero " << ++bla << " during\n";
			}

			/**
			\return true if correction loop should terminate, false otherwise.
			*/
			template<typename ANumT, typename XNumT>
			static bool Correct(std::vector<XNumT> & x, std::vector<ANumT> const& a)
			{
				std::vector<XNumT> residuals(PolyT::Degree);

				PolyT::EvaluateHomotopy(residuals, x, a);

				// for (unsigned ii=0; ii<PolyT::Degree; ++ii)
				// 	residuals[ii] = -residuals[ii];

				std::vector<XNumT> delta_x;
				PolyT::Predictor::LinSolver::Solve(delta_x, x, residuals);
				
				for (unsigned ii=0; ii<PolyT::Degree; ++ii)
					x[ii] -= delta_x[ii];
				
				for (typename std::vector<typename NumTraits<XNumT>::ComplexType>::const_iterator iter = delta_x.begin(); iter!=delta_x.end(); ++iter)
				{
					if (abs(*iter) > NumTraits<XNumT>::NewtonTerminationThreshold())
						return false;
				}

				return true;

			}


		}; // re: FactorizerBase
		









		template<typename PredictorT, template<int,int> class StartT>
		struct Octic : public FactorizerBase<Octic<PredictorT, StartT >, StartT<4,4> >
		{
			enum {
				  Degree = 8,
				  DegreeFactorR = 4,
				  DegreeFactorS = 4
				 };

			typedef FactorizerBase<Octic<PredictorT, StartT >, StartT<4,4> > Factorizer;
			typedef PredictorT Predictor;

			
			Octic(unsigned num_steps, unsigned num_corrects_during, unsigned num_corrects_after) : Factorizer(num_steps, num_corrects_during, num_corrects_after)
			{
			}

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

				f.resize(Degree);
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

				f.resize(Degree);
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
		template<typename PredictorT, template<int,int> class StartT>
		struct Decic : public FactorizerBase<Decic<PredictorT, StartT >, StartT<8,2> >
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
