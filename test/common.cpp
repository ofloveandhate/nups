// This file is part of Numerical Univariate Polynomial Solver (NUPS).

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


// this source file contains the test suite for NUPS
#define BOOST_TEST_DYN_LINK 1
#include <boost/test/unit_test.hpp>


#include "nups/nups.hpp"

BOOST_AUTO_TEST_SUITE(common_functionality)


using namespace nups;

BOOST_AUTO_TEST_CASE(NChoosek)
{
	BOOST_CHECK_EQUAL(NChooseK<double>(1,0),1);
	BOOST_CHECK_EQUAL(NChooseK<double>(2,0),1);
	BOOST_CHECK_EQUAL(NChooseK<double>(3,0),1);
	BOOST_CHECK_EQUAL(NChooseK<double>(4,0),1);
	BOOST_CHECK_EQUAL(NChooseK<double>(30,0),1);

	BOOST_CHECK_EQUAL(NChooseK<double>(1,1),1);
	BOOST_CHECK_EQUAL(NChooseK<double>(2,1),2);
	BOOST_CHECK_EQUAL(NChooseK<double>(3,1),3);
	BOOST_CHECK_EQUAL(NChooseK<double>(4,1),4);
	BOOST_CHECK_EQUAL(NChooseK<double>(30,1),30);

	BOOST_CHECK_EQUAL(NChooseK<double>(1,1),1);
	BOOST_CHECK_EQUAL(NChooseK<double>(2,2),1);
	BOOST_CHECK_EQUAL(NChooseK<double>(3,3),1);
	BOOST_CHECK_EQUAL(NChooseK<double>(4,4),1);
	BOOST_CHECK_EQUAL(NChooseK<double>(30,30),1);


	BOOST_CHECK_EQUAL(NChooseK<double>(4,2),6);
	BOOST_CHECK_EQUAL(NChooseK<double>(5,2),10);

}

BOOST_AUTO_TEST_SUITE_END()
