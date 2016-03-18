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

// declares special linear algebra solvers

#ifndef NUPS_LINEAR_SOLVE_HPP
#define NUPS_LINEAR_SOLVE_HPP

namespace nups {
	namespace solver {


		template<typename PolyT>
		struct LinearSolverBase
		{
			template<typename NumT, 
								template <typename, typename...> class ContT, 
			          			typename... Alloc>
			static void Solve(ContT<NumT, Alloc...> & solution, ContT<NumT, Alloc...> const& variables, ContT<NumT, Alloc...> const& rhs)
			{
				return PolyT::DoSolve(solution, variables, rhs);
			}
		};

		struct QuarticDenseLinear : public LinearSolverBase<QuarticDenseLinear>
		{
			enum {Degree = 4};
			// instructions on how to solve a dense 4x4 linear system
			// implemented using Sarrus' rule.
			//
			// M is input as a 16-length vector, but interpreted as a 4x4 matrix.  The first four entries are the top row, etc.
			template<typename NumT, 
								template <typename, typename...> class ContT, 
			          			typename... Alloc>
			static void DoSolve(ContT<NumT, Alloc...> & x, ContT<NumT, Alloc...> const& M, ContT<NumT, Alloc...> const& B)
			{
				#ifndef NUPS_DISABLE_ASSERTS
				assert(M.size()==16 && "must have 16 coefficients for the solution of a dense 4x4 system");
				assert(B.size()==4 && "must have 4 entries in the right hand side for the solution of a dense 4x4 system");
				#endif

				x.resize(4);

				const NumT& a = M[0];
				const NumT& b = M[1];
				const NumT& c = M[2];
				const NumT& d = M[3];
				const NumT& e = M[4];
				const NumT& f = M[5];
				const NumT& g = M[6];
				const NumT& h = M[7];
				const NumT& i = M[8];
				const NumT& j = M[9];
				const NumT& k = M[10];
				const NumT& l = M[11];
				const NumT& m = M[12];
				const NumT& n = M[13];
				const NumT& o = M[14];
				const NumT& p = M[15];

				const NumT& b1 = B[0];
				const NumT& b2 = B[1];
				const NumT& b3 = B[2];
				const NumT& b4 = B[3];

				NumT bg = b*g;
				NumT bh = b*h;
				NumT cf = c*f;
				NumT ch = c*h;
				NumT df = d*f;
				NumT dg = d*g;

				NumT kp = k*p;
				NumT lo = l*o;
				NumT jp = j*p;
				NumT ln = l*n;
				NumT jo = j*o;
				NumT kn = k*n;

				NumT D_234 = f*(kp-lo) + g*(ln-jp) + h*(jo-kn);
				NumT D_134 = b*(kp-lo) + c*(ln-jp) + d*(jo-kn);
				NumT D_124 = (bg-cf)*p + (df-bh)*o + (ch-dg)*n;
				NumT D_123 = (bg-cf)*l + (df-bh)*k + (ch-dg)*j;

				NumT det_M = a*D_234 - e*D_134 + i*D_124 - m*D_123;


				NumT sub134 = b*(kp-lo) + c*(ln-jp) + d*(jo-kn);
				NumT sub124 = (bg-cf)*p + (df-bh)*o + (ch-dg)*n;
				NumT sub123 = (bg-cf)*l + (df-bh)*k + (ch-dg)*j;
				NumT D1 = b1*D_234 - b2*sub134 + b3*sub124 - b4*sub123;


				NumT ag = a*g;
				NumT ah = a*h;
				NumT ce = c*e;
				NumT de = d*e;
				NumT km = k*m;
				NumT io = i*o;
				NumT lm = l*m;
				NumT ip = i*p;
				NumT 
				sub234 = e*(kp-lo) + g*(lm-ip) + h*(io-km);
				sub134 = a*(kp-lo) + c*(lm-ip) + d*(io-km);
				sub124 = (ag-ce)*p + (de-ah)*o + (ch-dg)*m;
				sub123 = (ag-ce)*l + (de-ah)*k + (ch-dg)*i;
				NumT D2 = -b1*sub234 + b2*sub134 - b3*sub124 + b4*sub123;

				NumT in = i*n;
				NumT jm = j*m;
				NumT be = b*e;
				NumT af = a*f;

				sub234 = e*(jp-ln) + f*(lm-ip) + h*(in-jm);
				sub134 = a*(jp-ln) + b*(lm-ip) + d*(in-jm);
				sub124 = (af-be)*p + (de-ah)*n + (bh-df)*m;
				sub123 = (af-be)*l + (de-ah)*j + (bh-df)*i;
				NumT D3 =  b1*sub234 - b2*sub134 + b3*sub124 - b4*sub123;

				sub234 = e*(jo-kn) + f*(km-io) + g*(in-jm);
				sub134 = a*(jo-kn) + b*(km-io) + c*(in-jm);
				sub124 = (af-be)*o + (ce-ag)*n + (bg-cf)*m;
				sub123 = (af-be)*k + (ce-ag)*j + (bg-cf)*i;
				NumT D4 = -b1*sub234 + b2*sub134 - b3*sub124 + b4*sub123;

				x[0] = D1 / det_M;
				x[1] = D2 / det_M;
				x[2] = D3 / det_M;
				x[3] = D4 / det_M;
			}

		};



		struct OcticLinear : public LinearSolverBase<OcticLinear>
		{
			enum {Degree = 8};
			/**
			\brief Solves the linear problem Ax = b, for a special matrix A.

			\param[out] solution The returned value of the solution.
			\param[in] r The values of the \f$r\f$ variables.
			\param[in] s The values of the \f$s\f$ variables.
			\param[in] rhs The right-hand-side of the linear problem being solved.

			In this particular predictor, the matrix A is assumed to have a block structure.  The structure comes from the product \f$r \cdot s\f$, where \f$r\f$ and \f$s\f$ are monic quartic polynomials.  That is, the code here has specific instructions on how to solve the special linear system for the factoring of an octic monic univariate into two monic quartics
			*/
			template<typename XNumT, typename RHSNumT, 
								template <typename, typename...> class ContT, 
			          			typename... Alloc1, typename... Alloc2>
			static void DoSolve(ContT<XNumT, Alloc1...> & solution, ContT<XNumT, Alloc1...> const& variables, ContT<RHSNumT, Alloc2...> const& rhs)
			{
				if (variables.size()!=8)
					throw std::runtime_error("there must be 8 variables in the input for octic linear solver");
				if (rhs.size()!=8)
					throw std::runtime_error("there must be 8 variables in the right hand side for octic linear solver");

				const XNumT& r0 = variables[0];
				const XNumT& r1 = variables[1];
				const XNumT& r2 = variables[2];
				const XNumT& r3 = variables[3];

				const XNumT& s0 = variables[4];
				const XNumT& s1 = variables[5];
				const XNumT& s2 = variables[6];
				const XNumT& s3 = variables[7];
				
				const RHSNumT& b1 = rhs[0];
				const RHSNumT& b2 = rhs[1];
				const RHSNumT& b3 = rhs[2];
				const RHSNumT& b4 = rhs[3];

				const RHSNumT& c1 = rhs[4];
				const RHSNumT& c2 = rhs[5];
				const RHSNumT& c3 = rhs[6];
				const RHSNumT& c4 = rhs[7];

				// first we will solve the problem in the top of the block matrix.

				//[[C D]  [x] = [b]
				//  A B]] [y]   [c]
				//
				//
				//  since A and B invertible, Ax+By=b is solvable by inverting A.
				//  
				//     Ax+By = c
				//        Ax = c-By
				//         x = A_inv(c-By)
				//  hence, 
				//      x = A^{-1}(c - B y)
				//  then substitute into the C,D,c equation

				solution.resize(8);



				//set up some temporaries			
				const XNumT z31 =  s3*s3-s2;
				const XNumT z41 = -s3*s3*s3+XNumT(2)*s2*s3-s1;

				//                        Cx+Dy = b
				//        C(A^{-1}(c - B y))+Dy = b
				// C A_inv c - C A_inv B y + Dy = b
				//           (-C A_inv B + D) y = b - C A_inv c
				// 
				// so we need to write the matrix
				//     (-C A_inv B + D)
				// into a matrix, and pass it into the dense 4x4 solver.
				// also, write 
				//      b - C A_inv c
				// into a vector, and pass into the same.

				

				// compute $C A^{-1}$

				// C A^{-1}
				ContT<XNumT, Alloc1...> C_Ainv(16);
				// C_Ainv = C A^{-1}
				C_Ainv[15] = s0-s1*s3 + s2*z31 + s3*z41;
				C_Ainv[14] = s1-s2*s3 + s3*z31;
				C_Ainv[13] = -z31;
				C_Ainv[12] = s3;
				//---------------------------------------
				C_Ainv[11] = -s0*s3 + s1*z31 + s2*z41;
				C_Ainv[10] = s0-s1*s3 + s2*z31;
				C_Ainv[9] = -s2*s3 + s1;
				C_Ainv[8] = s2;
				//---------------------------------------
				C_Ainv[7] = s0*z31 + s1*z41;
				C_Ainv[6] = -s0*s3 + s1*z31;
				C_Ainv[5] = -s1*s3 + s0;
				C_Ainv[4] = s1;
				//---------------------------------------
				C_Ainv[3] = s0*z41;
				C_Ainv[2] = s0*z31;
				C_Ainv[1] = -s0*s3;
				C_Ainv[0] = s0;
				//---------------------------------------

				const ContT<XNumT, Alloc1...>& w = C_Ainv;

				// M = -C A^{-1} B + D
				ContT<XNumT, Alloc1...> M(16);	

				M[15]  = -(w[14]*r3+w[13]*r2+w[12]*r1+w[15])	+ r0;
				M[14]  = -(w[13]*r3+w[12]*r2+w[14])				+ r1;
				M[13]  = -(w[12]*r3+w[13])						+ r2;
				M[12]  = -(w[12])								+ r3;
				//------------------------------------------
				M[11]  = -(w[10]*r3+w[9]*r2+w[8]*r1+w[11]);
				M[10]  = -(w[9]*r3+w[8]*r2+w[10])				+ r0;
				M[9]  = -(w[8]*r3+w[9])						+ r1;
				M[8]  = -(w[8])								+ r2;
				//------------------------------------------
				M[7]  = -(w[5]*r2+w[4]*r1+w[6]*r3+w[7]);
				M[6]  = -(w[5]*r3+w[4]*r2+w[6]);
				M[5] = -(w[4]*r3+w[5])					+ r0;
				M[4] = -(w[4])							+ r1;
				//------------------------------------------
				M[3] = -(w[2]*r3+w[1]*r2+w[0]*r1+w[3]);
				M[2] = -(w[1]*r3+w[0]*r2+w[2]);
				M[1] = -(w[0]*r3+w[1]);
				M[0] = -(w[0])							+ r0;
				//------------------------------------------

				// copy in the right hand side to a vector for the trailing 4x4 dense solve
				//

				// essentially b - C A_inv c
				ContT<XNumT, Alloc1...> rhs_44(4);
				rhs_44[0] = b1 -  (w[0]*c1 +  w[1]*c2 +  w[2]*c3 +  w[3]*c4);
				rhs_44[1] = b2 -  (w[4]*c1 +  w[5]*c2 +  w[6]*c3 +  w[7]*c4);
				rhs_44[2] = b3 -  (w[8]*c1 +  w[9]*c2 + w[10]*c3 + w[11]*c4);
				rhs_44[3] = b4 - (w[12]*c1 + w[13]*c2 + w[14]*c3 + w[15]*c4);
				
				ContT<XNumT, Alloc1...> temp(4);
				QuarticDenseLinear::Solve(temp, M,rhs_44);

				solution[4] = temp[0];
				solution[5] = temp[1];
				solution[6] = temp[2];
				solution[7] = temp[3];

				const XNumT& y1 = temp[0];
				const XNumT& y2 = temp[1];
				const XNumT& y3 = temp[2];
				const XNumT& y4 = temp[3];

				// then finally compute the first four entries in the solution.  they're easy at this point

				// we will re-use the following temporaries			
				// const NumT z31 =  pow(s3,2)-s2;
				// const NumT z41 = -pow(s3,3)+NumT(2)*s2*s3-s1;

				const XNumT t1 = -y4*r3 + c3-y3;
				const XNumT t2 = -r2*y4 - r3*y3 + c2 - y2;

				//      x = A^{-1}(c - B y)
				solution[3] = c4-y4;

				solution[0] = -r1*y4-r2*y3-r3*y2+c1-y1 - s3*t2 + z31*t1 + z41*solution[3];
				solution[1] = t2 - s3*t1 + z31*solution[3];
				solution[2] = t1 - s3*solution[3];
				
				
			}

		};


		struct DecicLinear : public LinearSolverBase<DecicLinear>
		{
			enum {Degree = 10};
			/**
			\brief Solves the linear problem Ax = b, for a special matrix A.

			\param[out] solution The returned value of the solution.
			\param[in] r The values of the \f$r\f$ variables.
			\param[in] s The values of the \f$s\f$ variables.
			\param[in] rhs The right-hand-side of the linear problem being solved.

			In this particular predictor, the matrix A is assumed to have a block structure.  The structure comes from the product \f$r \cdot s\f$, where \f$r\f$ and \f$s\f$ are monic quartic polynomials.  That is, the code here has specific instructions on how to solve the special linear system for the factoring of an octic monic univariate into two monic quartics
			*/
			template<typename NumT, 
								template <typename, typename...> class ContT, 
			          			typename... Alloc>
			static void DoSolve(ContT<NumT, Alloc...> & solution, ContT<NumT, Alloc...> const& variables, ContT<NumT, Alloc...> const& rhs)
			{
				
			}

		};


	} // re: solve
} // re: nups



#endif

