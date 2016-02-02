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

		struct QuarticDenseLinearSolver : public LinearSolverBase<QuarticDenseLinearSolver>
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

				NumT& a = M[0];
				NumT& b = M[1];
				NumT& c = M[2];
				NumT& d = M[3];
				NumT& e = M[4];
				NumT& f = M[5];
				NumT& g = M[6];
				NumT& h = M[7];
				NumT& i = M[8];
				NumT& j = M[9];
				NumT& k = M[10];
				NumT& l = M[11];
				NumT& m = M[12];
				NumT& n = M[13];
				NumT& o = M[14];
				NumT& p = M[15];

				NumT& b1 = B[0];
				NumT& b2 = B[1];
				NumT& b3 = B[2];
				NumT& b4 = B[3];

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

				NumT D_234 = f*(kp-lo) - g*(jp+ln) + h*(jo-kn);
				NumT D_134 = b*(kp-lo) - c*(jp+ln) + d*(jo-kn);

				NumT D_124 = (bg-cf)*p - (bh+df)*o + (ch-dg)*n;
				NumT D_123 = (bg-cf)*l - (bh+df)*k + (ch-dg)*j;

				NumT det_M = a*D_234 - e*D_134 + i*D_124 - m*D_123;


				NumT sub134 = b*(kp-lo) - c*(jp+ln) + d*(jo-kn);
				NumT sub124 = (bg-cf)*p - (bh+df)*o - (dg+ch)*n;
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
				sub234 = e*(kp-lo) - g*(ip+lm) + h*(io-km);
				sub134 = a*(kp-lo) - c*(ip+lm) + d*(io-km);
				sub124 = (ag-ce)*p + (de-ah)*o + (ch-dg)*m;
				sub123 = (ag-ce)*l + (de-ah)*k + (ch-dg)*i;
				NumT D2 = -b1*sub234 + b2*sub134 - b3*sub124 + b4*sub123;

				NumT in = i*n;
				NumT jm = j*m;
				NumT be = b*e;
				NumT af = a*f;

				sub234 = e*(jp-ln) - f*(ip+lm) + h*(in-jm);
				sub134 = a*(jp-ln) - b*(ip+lm) + d*(in-jm);
				sub124 = (af-be)*p + (de-ah)*n + (bh-df)*m;
				sub123 = (af-be)*l + (de-ah)*j + (bh-df)*i;
				NumT D3 =  b1*sub234 - b2*sub134 + b3*sub124 - b4*sub123;

				sub234 = e*(jo-kn) - f*(io+km) + g*(in-jm);
				sub134 = a*(jo-kn) - b*(io+km) + c*(in-jm);
				sub124 = (af-be)*o + (ce-ag)*n + (bg-cf)*m;
				sub123 = (af-be)*k + (ce-ag)*j + (bg-cf)*i;
				NumT D4 = -b1*sub234 + b2*sub134 - b3*sub124 + b4*sub123;

				x[0] = D1 / det_M;
				x[1] = D2 / det_M;
				x[2] = D3 / det_M;
				x[3] = D4 / det_M;
			}

		};



		struct OcticLinearSolver : public LinearSolverBase<OcticLinearSolver>
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

