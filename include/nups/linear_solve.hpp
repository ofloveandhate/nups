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
			template<typename NumT>
			static void Solve(std::vector<NumT> & solution, std::vector<NumT> const& variables, std::vector<NumT> const& rhs)
			{
				return PolyT::DoSolve(solution, variables, rhs);
			}
		};

		struct QuarticDenseLinear : public LinearSolverBase<QuarticDenseLinear>
		{

			// instructions on how to solve a dense 4x4 linear system
			// implemented using Sarrus' rule.
			//
			// M is input as a 16-length vector, but interpreted as a 4x4 matrix.  The first four entries are the top row, etc.
			template<typename NumT>
			static void DoSolve(std::vector<NumT> & x, std::vector<NumT> const& M, std::vector<NumT> const& B)
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

			/**
			\brief Solves the linear problem Ax = b, for a special matrix A.

			\param[out] solution The returned value of the solution.
			\param[in] r The values of the \f$r\f$ variables.
			\param[in] s The values of the \f$s\f$ variables.
			\param[in] rhs The right-hand-side of the linear problem being solved.

			In this particular predictor, the matrix A is assumed to have a block structure.  The structure comes from the product \f$r \cdot s\f$, where \f$r\f$ and \f$s\f$ are monic quartic polynomials.  That is, the code here has specific instructions on how to solve the special linear system for the factoring of an octic monic univariate into two monic quartics
			*/
			template<typename NumT>
			static void DoSolve(std::vector<NumT> & solution, std::vector<NumT> const& variables, std::vector<NumT> const& rhs)
			{
				if (variables.size()!=8)
					throw std::runtime_error("there must be 8 variables in the input for octic linear solver");
				if (rhs.size()!=8)
					throw std::runtime_error("there must be 8 variables in the right hand side for octic linear solver");

				const NumT& r3 = variables[0];
				const NumT& r2 = variables[1];
				const NumT& r1 = variables[2];
				const NumT& r0 = variables[3];

				const NumT& s3 = variables[4];
				const NumT& s2 = variables[5];
				const NumT& s1 = variables[6];
				const NumT& s0 = variables[7];

				const NumT& b1 = rhs[0];
				const NumT& b2 = rhs[1];
				const NumT& b3 = rhs[2];
				const NumT& b4 = rhs[3];

				const NumT& c1 = rhs[4];
				const NumT& c2 = rhs[5];
				const NumT& c3 = rhs[6];
				const NumT& c4 = rhs[7];

				// first we will solve the problem in the top of the block matrix.

				//[[A B]  [x] = [b]
				//  C D]] [y]   [c]
				//
				//
				//  since A and B invertible, Ax+By=b is solvable by inverting A.
				//  
				//     Ax+By = b
				//        Ax = b-By
				//         x = A_inv(b-By)
				//  hence, 
				//      x = A^{-1}(b - B y)
				//  then substitute into the C,D,c equation

				solution.resize(8);

				// aq_inv:=aq^(-1):
				// x:=aq_inv.(B-bq.C);
				//                Matrix(%id = 18446744078454429446)
				// lprint(x)
				// Matrix(4, 1, {
				//(1, 1) = b1-c1, 
				//(2, 1) = -s3*(b1-c1)-r3*c1+b2-c2, 
				//(3, 1) = (s3^2-s2)*(b1-c1)-s3*(-c1*r3+b2-c2)-r2*c1-r3*c2+b3-c3, 
				//(4, 1) = (-s3^3+2*s2*s3-s1)*(b1-c1)+(s3^2-s2)*(-c1*r3+b2-c2)-s3*(-c1*r2-c2*r3+b3-c3)-r1*c1-r2*c2-r3*c3+b4-c4}, datatype = anything, storage = rectangular, order = Fortran_order, shape = [])


				//set up some temporaries			
				const NumT z31 =  pow(s3,2)-s2;
				const NumT z41 = -pow(s3,3)+NumT(2)*s2*s3-s1;

				//                        Cx+Dy = c
				//        C(A^{-1}(b - B y))+Dy = c
				// C A_inv b - C A_inv B y + Dy = c
				//           (-C A_inv B + D) y = c - C A_inv b
				// 
				// so we need to write the matrix
				//     (-C A_inv B + D)
				// into a matrix, and pass it into the dense 4x4 solver.
				// also, write 
				//      c - C A_inv b
				// into a vector, and pass into the same.

				

				// compute $C A^{-1}$

				// C A^{-1}
				// Matrix(4, 4, {
				//(1, 1) = s0-s1*s3+s2*(s3^2-s2)+s3*(-s3^3+2*s2*s3-s1), 
				//(1, 2) = s1-s2*s3+s3*(s3^2-s2), 
				//(1, 3) = -s3^2+s2, 
				//(1, 4) = s3, 
				//
				//(2, 1) = -s0*s3+s1*(s3^2-s2)+s2*(-s3^3+2*s2*s3-s1), 
				//(2, 2) = s0-s1*s3+s2*(s3^2-s2), 
				//(2, 3) = -s2*s3+s1, 
				//(2, 4) = s2, 
				//
				//(3, 1) = s0*(s3^2-s2)+s1*(-s3^3+2*s2*s3-s1), 
				//(3, 2) = -s0*s3+s1*(s3^2-s2), 
				//(3, 3) = -s1*s3+s0, 
				//(3, 4) = s1, 
				//
				//(4, 1) = s0*(-s3^3+2*s2*s3-s1), 
				//(4, 2) = s0*(s3^2-s2), 
				//(4, 3) = -s0*s3, 
				//(4, 4) = s0}, datatype = anything, storage = rectangular, order = Fortran_order, shape = [])

				std::vector<NumT> w(16);
				// w = A^{-1}
				w[0] = s0-s1*s3 + s2*z31 + s3*z41;
				w[1] = s1-s2*s3 + s3*z31;
				w[2] = -pow(s3,2) + s2;
				w[3] = s3;
				//---------------------------------------
				w[4] = -s0*s3 + s1*z31 + s2*z41;
				w[5] = s0-s1*s3 + s2*z31;
				w[6] = -s2*s3 + s1;
				w[7] = s2;
				//---------------------------------------
				w[8] = s0*z31 + s1*z41;
				w[9] = -s0*s3 + s1*z31;
				w[10] = -s1*s3 + s0;
				w[11] = s1;
				//---------------------------------------
				w[12] = s0*z41;
				w[13] = s0*z31;
				w[14] = -s0*s3;
				w[15] = s0;
				//---------------------------------------



				// the following is a dense matrix times B, with B written in terms of the coefficients r_i and s_i.
				//
				//
				// Matrix(4, 4, {(1, 1) = m1*r3+m2*r2+m3*r1+m0, (1, 2) = m2*r3+m3*r2+m1, (1, 3) = m3*r3+m2, (1, 4) = m3, (2, 1) = m5*r3+m6*r2+m7*r1+m4, (2, 2) = m6*r3+m7*r2+m5, (2, 3) = m7*r3+m6, (2, 4) = m7, (3, 1) = m10*r2+m11*r1+m9*r3+m8, (3, 2) = m10*r3+m11*r2+m9, (3, 3) = m11*r3+m10, (3, 4) = m11, (4, 1) = m13*r3+m14*r2+m15*r1+m12, (4, 2) = m14*r3+m15*r2+m13, (4, 3) = m15*r3+m14, (4, 4) = m15}, datatype = anything, storage = rectangular, order = Fortran_order, shape = []
				//
				//

				std::vector<NumT> M(16);	

				M[0]  = -(w[1]*r3+w[2]*r2+w[3]*r1+w[0])		+ r0;
				M[1]  = -(w[2]*r3+w[3]*r2+w[1])				+ r1;
				M[2]  = -(w[3]*r3+w[2])						+ r2;
				M[3]  = -(w[3])								+ r3;
				//------------------------------------------
				M[4]  = -(w[5]*r3+w[6]*r2+w[7]*r1+w[4]);
				M[5]  = -(w[6]*r3+w[7]*r2+w[5])				+ r0;
				M[6]  = -(w[7]*r3+w[6])						+ r1;
				M[7]  = -(w[7])								+ r2;
				//------------------------------------------
				M[8]  = -(w[10]*r2+w[11]*r1+w[9]*r3+w[8]);
				M[9]  = -(w[10]*r3+w[11]*r2+w[9]);
				M[10] = -(w[11]*r3+w[10])					+ r0;
				M[11] = -(w[11])							+ r1;
				//------------------------------------------
				M[12] = -(w[13]*r3+w[14]*r2+w[15]*r1+w[12]);
				M[13] = -(w[14]*r3+w[15]*r2+w[13]);
				M[14] = -(w[15]*r3+w[14]);
				M[15] = -(w[15])							+ r0;
				//------------------------------------------

				// copy in the right hand side to a vector for the trailing 4x4 dense solve
				//

				// essentially c - C A_inv b
				std::vector<NumT> y(4);
				y[0] = c1 -  (w[0]*b1 +  w[1]*b2 +  w[2]*b3 +  w[3]*b4);
				y[1] = c2 -  (w[4]*b1 +  w[5]*b2 +  w[6]*b3 +  w[7]*b4);
				y[2] = c3 -  (w[8]*b1 +  w[9]*b2 + w[10]*b3 + w[11]*b4);
				y[3] = c4 - (w[12]*b1 + w[13]*b2 + w[14]*b3 + w[15]*b4);
				
				std::vector<NumT> temp(4);
				QuarticDenseLinear::Solve(temp, M,y);

				solution[4] = temp[0];
				solution[5] = temp[1];
				solution[6] = temp[2];
				solution[7] = temp[3];

				const NumT& y1 = temp[0];
				const NumT& y2 = temp[1];
				const NumT& y3 = temp[2];
				const NumT& y4 = temp[3];

				// then finally compute the first four entries in the solution.  they're easy at this point
				// Matrix(4, 1, {
				//(1, 1) = b1-y1, 
				//(2, 1) = -s3*(b1-y1)-r3*y1+b2-y2, 
				//(3, 1) = (s3^2-s2)*(b1-y1)-s3*(-r3*y1+b2-y2)-r2*y1-r3*y2+b3-y3, 
				//(4, 1) = (-s3^3+2*s2*s3-s1)*(b1-y1)+(s3^2-s2)*(-r3*y1+b2-y2)-s3*(-r2*y1-r3*y2+b3-y3)-r1*y1-r2*y2-r3*y3+b4-y4}, datatype = anything, storage = rectangular, order = Fortran_order, shape = [])

				// we will re-use the following temporaries			
				// const NumT z31 =  pow(s3,2)-s2;
				// const NumT z41 = -pow(s3,3)+NumT(2)*s2*s3-s1;

				const NumT t1 = -y1*r3 + b2-y2;
				const NumT t2 = -y1*r2 - y2*r3 + b3-y3;


				solution[0] =  b1-y1;
				solution[1] = -s3*solution[0] - r3*y1+b2-y2;
				solution[2] = z31*solution[0] - s3*t1 + t2;
				solution[3] = z41*solution[0] + z31*t1 - s3*t2 - r1*y1 - r2*y2 - r3*y3 + b4-y4;
			}

		};


		struct DecicLinear : public LinearSolverBase<DecicLinear>
		{

			/**
			\brief Solves the linear problem Ax = b, for a special matrix A.

			\param[out] solution The returned value of the solution.
			\param[in] r The values of the \f$r\f$ variables.
			\param[in] s The values of the \f$s\f$ variables.
			\param[in] rhs The right-hand-side of the linear problem being solved.

			In this particular predictor, the matrix A is assumed to have a block structure.  The structure comes from the product \f$r \cdot s\f$, where \f$r\f$ and \f$s\f$ are monic quartic polynomials.  That is, the code here has specific instructions on how to solve the special linear system for the factoring of an octic monic univariate into two monic quartics
			*/
			template<typename NumT>
			static void DoSolve(std::vector<NumT> & solution, std::vector<NumT> const& variables, std::vector<NumT> const& rhs)
			{
				
			}

		};


	} // re: solve
} // re: nups



#endif

