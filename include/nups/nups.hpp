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

// this is the primary header file for NUPS

#ifndef NUPS_HPP
#define NUPS_HPP

#include "nups/factor.hpp"
#include "nups/polynomial_solve.hpp"
#include "nups/linear_solve.hpp"
#include "nups/predict.hpp"
#include <nups_exports.h>


extern "C"
{
	/**
	\brief Check whether NUPS is available.  

	This function is intended for use with the Autotools, particularly the AC_SEARCH_LIBS command.

	\return The character `y`.
	*/
	char HaveNUPS();
}

#endif
