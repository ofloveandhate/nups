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




// \todo Make the DYN_LINK change depending on the targeted architecture.  some need it, others don't.
//if used, this BOOST_TEST_DYN_LINK appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_DYN_LINK 1

//this #define MUST appear before #include <boost/test/unit_test.hpp>
#define BOOST_TEST_MODULE "NUPS Testing"
#include <boost/test/unit_test.hpp>

#include <ctime>

#include <iostream>
#include "nups/nups.hpp"

using std::abs;
// #define NO_RANDOM 1

typedef std::complex<double> dbl;
using namespace nups;



template<typename T>
void RandomTestOctic(std::vector<T> & coefficients, std::vector<T> & solutions)
{
	solutions.resize(8);
	for (int ii=0; ii<8; ii++)
		solutions[ii] = Random<T>::Generate();


	const T& s0 = solutions[0];
	const T& s1 = solutions[1];
	const T& s2 = solutions[2];
	const T& s3 = solutions[3];
	const T& s4 = solutions[4];
	const T& s5 = solutions[5];
	const T& s6 = solutions[6];
	const T& s7 = solutions[7];

	coefficients.resize(8); // omitting the 1, so it's monic

	coefficients[0] = s0*s1*s2*s3*s4*s5*s6*s7;
	coefficients[1] = - s7*(s6*(s5*(s4*(s3*(s2*(s0 + s1) + s0*s1) + s0*s1*s2) + s0*s1*s2*s3) + s0*s1*s2*s3*s4) + s0*s1*s2*s3*s4*s5) - s0*s1*s2*s3*s4*s5*s6;
	coefficients[2] = s7*(s6*(s4*(s3*(s2*(s0 + s1) + s0*s1) + s0*s1*s2) + s5*(s4*(s2*(s0 + s1) + s0*s1 + s3*(s0 + s1 + s2)) + s3*(s2*(s0 + s1) + s0*s1) + s0*s1*s2) + s0*s1*s2*s3) + s5*(s4*(s3*(s2*(s0 + s1) + s0*s1) + s0*s1*s2) + s0*s1*s2*s3) + s0*s1*s2*s3*s4) + s6*(s5*(s4*(s3*(s2*(s0 + s1) + s0*s1) + s0*s1*s2) + s0*s1*s2*s3) + s0*s1*s2*s3*s4) + s0*s1*s2*s3*s4*s5;
	coefficients[3] = - s7*(s4*(s3*(s2*(s0 + s1) + s0*s1) + s0*s1*s2) + s5*(s4*(s2*(s0 + s1) + s0*s1 + s3*(s0 + s1 + s2)) + s3*(s2*(s0 + s1) + s0*s1) + s0*s1*s2) + s6*(s5*(s2*(s0 + s1) + s0*s1 + s3*(s0 + s1 + s2) + s4*(s0 + s1 + s2 + s3)) + s4*(s2*(s0 + s1) + s0*s1 + s3*(s0 + s1 + s2)) + s3*(s2*(s0 + s1) + s0*s1) + s0*s1*s2) + s0*s1*s2*s3) - s6*(s4*(s3*(s2*(s0 + s1) + s0*s1) + s0*s1*s2) + s5*(s4*(s2*(s0 + s1) + s0*s1 + s3*(s0 + s1 + s2)) + s3*(s2*(s0 + s1) + s0*s1) + s0*s1*s2) + s0*s1*s2*s3) - s5*(s4*(s3*(s2*(s0 + s1) + s0*s1) + s0*s1*s2) + s0*s1*s2*s3) - s0*s1*s2*s3*s4;
	coefficients[4] = s7*(s5*(s2*(s0 + s1) + s0*s1 + s3*(s0 + s1 + s2) + s4*(s0 + s1 + s2 + s3)) + s4*(s2*(s0 + s1) + s0*s1 + s3*(s0 + s1 + s2)) + s6*(s2*(s0 + s1) + s0*s1 + s5*(s0 + s1 + s2 + s3 + s4) + s3*(s0 + s1 + s2) + s4*(s0 + s1 + s2 + s3)) + s3*(s2*(s0 + s1) + s0*s1) + s0*s1*s2) + s4*(s3*(s2*(s0 + s1) + s0*s1) + s0*s1*s2) + s5*(s4*(s2*(s0 + s1) + s0*s1 + s3*(s0 + s1 + s2)) + s3*(s2*(s0 + s1) + s0*s1) + s0*s1*s2) + s6*(s5*(s2*(s0 + s1) + s0*s1 + s3*(s0 + s1 + s2) + s4*(s0 + s1 + s2 + s3)) + s4*(s2*(s0 + s1) + s0*s1 + s3*(s0 + s1 + s2)) + s3*(s2*(s0 + s1) + s0*s1) + s0*s1*s2) + s0*s1*s2*s3;
	coefficients[5] = - s5*(s2*(s0 + s1) + s0*s1 + s3*(s0 + s1 + s2) + s4*(s0 + s1 + s2 + s3)) - s4*(s2*(s0 + s1) + s0*s1 + s3*(s0 + s1 + s2)) - s7*(s2*(s0 + s1) + s0*s1 + s5*(s0 + s1 + s2 + s3 + s4) + s3*(s0 + s1 + s2) + s6*(s0 + s1 + s2 + s3 + s4 + s5) + s4*(s0 + s1 + s2 + s3)) - s6*(s2*(s0 + s1) + s0*s1 + s5*(s0 + s1 + s2 + s3 + s4) + s3*(s0 + s1 + s2) + s4*(s0 + s1 + s2 + s3)) - s3*(s2*(s0 + s1) + s0*s1) - s0*s1*s2;
	coefficients[6] = s2*(s0 + s1) + s0*s1 + s7*(s0 + s1 + s2 + s3 + s4 + s5 + s6) + s5*(s0 + s1 + s2 + s3 + s4) + s3*(s0 + s1 + s2) + s6*(s0 + s1 + s2 + s3 + s4 + s5) + s4*(s0 + s1 + s2 + s3);
	coefficients[7] = - s0 - s1 - s2 - s3 - s4 - s5 - s6 - s7;
}

BOOST_AUTO_TEST_SUITE(NUPS_timing)

BOOST_AUTO_TEST_CASE(solve_random_octic_xxx_times)
{
	#ifndef NO_RANDOM
	srand((unsigned) time(NULL));
	#endif
	unsigned num_solves = 100;

	std::clock_t start = std::clock();
   

	std::vector<std::complex<double> > solutions;
	std::vector<double> expected_solutions, coefficients;

	double accuracy(1e-5);
	unsigned num_misses = 0;

	nups::solver::Octic<nups::predict::RKKC45> octic_solver(accuracy);


	for (unsigned solve_counter = 0; solve_counter < num_solves; solve_counter++)
	{
		RandomTestOctic(coefficients, expected_solutions);

		octic_solver.Solve(solutions, coefficients);

		for (int jj=0; jj<8; jj++)
		{
			unsigned soln_counter(0);
			for (int kk=0; kk<8; kk++)
			{
				
				if (abs(solutions[kk]/expected_solutions[jj]-1.) < accuracy
				    ||
				    abs(solutions[kk]-expected_solutions[jj]) < accuracy
				    )
					soln_counter++;
			}
			if (soln_counter!=1)
				num_misses++;
		}
	}

	std::cout << "solving " << num_solves << " random " << (NumTraits<double>::IsComplex==1? "complex" : "real") << " octics in double precision took " << (std::clock()-start ) / (double) CLOCKS_PER_SEC << " seconds \n";
	std::cout << "missed " << num_misses << " solutions, to an accuracy of " << accuracy << '\n';
}


BOOST_AUTO_TEST_CASE(solve_random_octic_xxx_times_single_precision)
{
	#ifndef NO_RANDOM
	srand((unsigned) time(NULL));
	#endif
	unsigned num_solves = 100;

	std::clock_t start = std::clock();

	std::vector<std::complex<float> > solutions;
	std::vector<float> expected_solutions, coefficients;

	float accuracy(1e-3f);
	unsigned num_misses = 0;

	for (unsigned solve_counter = 0; solve_counter < num_solves; solve_counter++)
	{
		RandomTestOctic(coefficients, expected_solutions);

		nups::solver::Octic<> octic_solver(accuracy);
		octic_solver.Solve(solutions, coefficients);

		for (int jj=0; jj<8; jj++)
		{
			unsigned soln_counter(0);
			for (int kk=0; kk<8; kk++)
			{
				
				if (float(abs(solutions[kk]/expected_solutions[jj])-1.) < accuracy)
					soln_counter++;
			}
			if (soln_counter!=1)
				num_misses++;
		}
	}

	std::cout << "solving " << num_solves << " random " << (NumTraits<float>::IsComplex==1? "complex" : "real") << " octics in single precision took " << (std::clock()-start ) / (double) CLOCKS_PER_SEC << " seconds \n";
	std::cout << "missed " << num_misses << " solutions, to an accuracy of " << accuracy << '\n';
}


BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(NUPS_linear_solvers)

 // -1.631156917435833 - 0.165573023171453i
 // -2.773590646243862 - 1.286827745115013i
 //  0.723863891338603 + 0.586140323200161i
 // -0.474105293913687 + 1.325907158894314i

 // -0.889662440444532 - 0.997784222450224i
 //  1.498528263775562 + 0.138648791678056i
 // -0.808347937232812 - 1.110773307588731i
 //  2.260707593911436 - 3.399377350698251i

 //  0.287478899831107 - 1.075179992397824i
 //  0.173647195021945 + 0.994050782471957i
 //  1.263365201250299 + 0.927310007556103i
 // -0.898300927134314 + 0.502004640423452i

 // -0.097002427585921 + 1.344007829606197i
 //  0.857260131416047 + 2.808087909665921i
 //  1.945996333272493 - 0.218837249349544i
 // -1.001455205130960 - 0.345932869010513i


// b =
//   0.865925787347672 - 1.240316296365440i
//   0.533826978448580 + 0.889176782061189i
//  -0.200582102021233 - 1.409805928876607i
//  -0.744875095706259 + 0.583146128751567i

// x = M\b;

// x =
//  -1.279349503109161 - 0.652455114334128i
//   0.062551338952494 + 0.872986846927085i
//  -1.241450919227002 - 0.354763623795993i
//  -0.107946252286943 + 0.771816786734144i
BOOST_AUTO_TEST_CASE(linear_solve_4x4)
{
	std::vector<std::complex<double> > M(16);
	M[0] = std::complex<double>(-1.631156917435833, + 0.165573023171453);
	M[1] = std::complex<double>(-2.773590646243862, + 1.286827745115013);
	M[2] = std::complex<double>( 0.723863891338603, - 0.586140323200161);
	M[3] = std::complex<double>(-0.474105293913687, - 1.325907158894314);

	M[4] = std::complex<double>(-0.889662440444532, + 0.997784222450224);
	M[5] = std::complex<double>( 1.498528263775562, - 0.138648791678056);
	M[6] = std::complex<double>(-0.808347937232812, + 1.110773307588731);
	M[7] = std::complex<double>( 2.260707593911436, + 3.399377350698251);

	M[8] = std::complex<double>( 0.287478899831107, + 1.075179992397824);
	M[9] = std::complex<double>( 0.173647195021945, - 0.994050782471957);
	M[10] = std::complex<double>( 1.263365201250299, - 0.927310007556103);
	M[11] = std::complex<double>(-0.898300927134314, - 0.502004640423452);

	M[12] = std::complex<double>(-0.097002427585921, - 1.344007829606197);
	M[13] = std::complex<double>( 0.857260131416047, - 2.808087909665921);
	M[14] = std::complex<double>( 1.945996333272493, + 0.218837249349544);
	M[15] = std::complex<double>(-1.001455205130960, + 0.345932869010513);

	std::vector<std::complex<double> > b(4);
	b[0] = dbl( 0.865925787347672, - 1.240316296365440);
	b[1] = dbl( 0.533826978448580, + 0.889176782061189);
	b[2] = dbl(-0.200582102021233, - 1.409805928876607);
	b[3] = dbl(-0.744875095706259, + 0.583146128751567);

	std::vector<std::complex<double> > solution;
	solver::QuarticDenseLinear::Solve(solution, M, b);

	BOOST_CHECK_EQUAL(solution.size(),4);

	std::vector<std::complex<double> > expected_soln(4);
	expected_soln[0] = dbl(-1.279349503109161, - 0.652455114334128);
	expected_soln[1] = dbl( 0.062551338952494, + 0.872986846927085);
	expected_soln[2] = dbl(-1.241450919227002, - 0.354763623795993);
	expected_soln[3] = dbl(-0.107946252286943, + 0.771816786734144);

	for (unsigned ii=0; ii<4; ii++)
		BOOST_CHECK(abs(expected_soln[ii]-solution[ii])<1e-10);

}	

BOOST_AUTO_TEST_CASE(linear_solve_8x8)
{

	std::vector<std::complex<double> > coeffs(8);

coeffs[0] = std::complex<double>(-0.3151061157083938,-1.3414324935486135);
coeffs[1] = std::complex<double>(0.5130251277984067,-0.9492387398685719);
coeffs[2] = std::complex<double>(-0.7175701576689141,-0.4503636396955308);
coeffs[3] = std::complex<double>(-1.5706863737242078,-1.2271908671607061);
coeffs[4] = std::complex<double>(-0.5266673941452300,-1.5078131850342797);
coeffs[5] = std::complex<double>(0.1343329839984859,1.3085501610686945);
coeffs[6] = std::complex<double>(-0.1757466132541140,-1.0279008815303350);
coeffs[7] = std::complex<double>(1.4615184079698247,1.3451143295000334);



std::vector<std::complex<double> > b(8);

b[0] = std::complex<double>(-0.2784213415124605,1.2672413358057288);
b[1] = std::complex<double>(0.1843749327615239,1.2396624590242373);
b[2] = std::complex<double>(-0.5963537326230182,-0.6881297149121884);
b[3] = std::complex<double>(-0.4971761886226711,-0.7590994249738381);
b[4] = std::complex<double>(0.7538997879176960,-0.8276617521751226);
b[5] = std::complex<double>(0.6463047638648889,0.5219204693674823);
b[6] = std::complex<double>(-0.4756291783769999,-1.0239159797433899);
b[7] = std::complex<double>(2.2253382909024038,-0.9128527205472096);



std::vector<std::complex<double> > expected_solution(8);

expected_solution[0] = std::complex<double>(2.6969283330365772,0.4184387508306890);
expected_solution[1] = std::complex<double>(-1.0855360585087270,0.0210950704334326);
expected_solution[2] = std::complex<double>(0.1993241582055479,-4.9826677670055286);
expected_solution[3] = std::complex<double>(3.7270121880678162,-3.4359598399925901);
expected_solution[4] = std::complex<double>(-2.9753496745490384,0.8488025849750411);
expected_solution[5] = std::complex<double>(-0.0167945039196800,1.7816416871906171);
expected_solution[6] = std::complex<double>(-2.8884713502448576,-3.6373103065930361);
expected_solution[7] = std::complex<double>(3.3801630124113076,-2.0651567154747932);



	std::vector<std::complex<double> > solution;
	solver::OcticLinear::Solve(solution, coeffs, b);

	BOOST_CHECK_EQUAL(solution.size(),8);

	for (unsigned ii=0; ii<8; ii++)
		BOOST_CHECK(abs(expected_solution[ii]/solution[ii]-1.)<1e-14);
}


BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(NUPS_factorizations)


BOOST_AUTO_TEST_CASE(factor_8_into_4by4)
{
	typedef nups::factor::Octic<nups::predict::RKKC45, nups::factor::UnitCoefficient > FactorT;

	std::vector<double > solution_r(4);
	std::vector<double > solution_s(4);

	solution_r[0] = double(0.4693906410582058);
	solution_r[1] = double(0.0119020695012414);
	solution_r[2] = double(0.3371226443988815);
	solution_r[3] = double(0.1621823081932428);

	solution_s[0] = double(0.7942845406839070);
	solution_s[1] = double(0.3112150420448049);
	solution_s[2] = double(0.5285331355062127);
	solution_s[3] = double(0.1656487294997809);


	std::vector<double> coefficients(8); // omitting the 1, so it's monic

	FactorT::EvaluateF(coefficients,solution_r,solution_s);


	std::vector<dbl> r, s;
	
	FactorT factorizer(20,5,5);
	factorizer.Factor(r,s,coefficients);

	BOOST_CHECK_EQUAL(r.size(),4);
	BOOST_CHECK_EQUAL(s.size(),4);

	std::vector<dbl> computed_prod_coeffs;
	FactorT::EvaluateF(computed_prod_coeffs,r,s);

	for (unsigned ii=0; ii<8; ++ii)
	{
		BOOST_CHECK(abs(computed_prod_coeffs[ii]-coefficients[ii])< 1e-14);
	}

}

BOOST_AUTO_TEST_SUITE_END()









BOOST_AUTO_TEST_SUITE(NUPS_solvers)


// a = rand_complex; s = zeros(8,1);
// for ii = 1:8
// s(ii) = rand_complex;
// a = a*(x-s(ii));
// end
// format long
// expand(vpa(a))
// s
// ans =
// x^8*(0.5102771270608572606874986377079 - 1.6528224144013987650936314821593i) + 
// x^7*(0.90080116997308669452903423530711 + 7.1654836636778492636366262896879i) - 
// x^6*(3.2421015927003609632114702833982 + 16.842917356017265417578826356585i) + 
// x^5*(19.115437042678237805862192791065 + 35.697059470028409128493009212141i) - 
// x^4*(61.198835175702023123328733919335 + 25.625494282276837277630792392328i) + 
// x^3*(72.192744000978711865061943895201 - 28.463232401133662670713566770604i) - 
// x^2*(35.643299557675072119635294047226 - 60.438826166075916406691691175645i) - 
//   x*(0.6695843516918776007318869197661 + 44.16834581802318752729978192059i) + 
//     (6.9362163213480266725291948886633 + 11.36510350183261168011047654863i)
// s =
//   1.020468846389024 - 0.994630198215936i
//   0.075088202598739 + 2.509347597404298i
//   0.537318199875311 + 0.631321535977640i
//   0.177740405586594 - 1.360379480013326i
//  -0.964554166387889 - 0.666930213646549i
//   0.965393415487517 - 0.696220582676325i
//   0.954957261427937 - 0.779301681708935i
//   1.038011326693907 - 0.362757386567121i
// BOOST_AUTO_TEST_CASE(solve_octic_0_real_roots)
// {
// 	#ifndef NO_RANDOM
// 	srand(time(NULL));
// 	#endif

// 	std::vector<std::complex<double> > solutions;
// 	std::vector<std::complex<double> > coefficients(9);

// 	coefficients[0] = std::complex<double>(0.5102771270608572606874986377079, -1.6528224144013987650936314821593);
// 	coefficients[1] = std::complex<double>(0.90080116997308669452903423530711, 7.1654836636778492636366262896879);
// 	coefficients[2] = std::complex<double>(-3.2421015927003609632114702833982, -16.842917356017265417578826356585);
// 	coefficients[3] = std::complex<double>(19.115437042678237805862192791065, 35.697059470028409128493009212141);
// 	coefficients[4] = std::complex<double>(-61.198835175702023123328733919335, -25.625494282276837277630792392328);
// 	coefficients[5] = std::complex<double>(72.192744000978711865061943895201, -28.463232401133662670713566770604);
// 	coefficients[6] = std::complex<double>(-35.643299557675072119635294047226, 60.438826166075916406691691175645);
// 	coefficients[7] = std::complex<double>(-0.6695843516918776007318869197661, -44.16834581802318752729978192059);
// 	coefficients[8] = std::complex<double>(6.9362163213480266725291948886633, 11.36510350183261168011047654863);



// 	nups::solver::Octic<
// 						// nups::predict::Euler
// 											>::Solve(solutions, coefficients);
// 	//
// 	BOOST_CHECK(solutions.size()==8);

// 	std::vector<std::complex<double> > expected_solutions(8);

// 	expected_solutions[0] = std::complex<double>(1.020468846389024, -0.994630198215936);
// 	expected_solutions[1] = std::complex<double>(0.075088202598739, 2.509347597404298);
// 	expected_solutions[2] = std::complex<double>(0.537318199875311, 0.631321535977640);
// 	expected_solutions[3] = std::complex<double>(0.177740405586594, -1.360379480013326);
// 	expected_solutions[4] = std::complex<double>(-0.964554166387889, -0.666930213646549);
// 	expected_solutions[5] = std::complex<double>(0.965393415487517, -0.696220582676325);
// 	expected_solutions[6] = std::complex<double>(0.954957261427937, -0.779301681708935);
// 	expected_solutions[7] = std::complex<double>(1.038011326693907, -0.362757386567121);

// 	for (int jj=0; jj<8; jj++)
// 		std::cout << solutions[jj] << "\t" << expected_solutions[jj] << std::endl;

// 	for (int jj=0; jj<8; jj++)
// 	{
// 		unsigned soln_counter(0);
// 		for (int ii=0; ii<8; ii++)
// 		{
// 			// std::cout << abs(solutions[ii] - expected_solutions[jj]) << std::endl;
// 			if (abs(solutions[ii]-expected_solutions[jj])<1e-5)
// 				soln_counter++;
// 		}
// 		BOOST_CHECK_EQUAL(soln_counter,1);
// 	}
// }





// ans =
// 0.95788953015050193329216199344955*x^8 - 4.616048869845871114437097477657*x^7 + 9.4118690334462455604524464083315*x^6 - 10.520090852574764399782523930068*x^5 + 6.9607710697687273625964490097213*x^4 - 2.7299923602832188539429250393216*x^3 + 0.59209079987668368905301009559429*x^2 - 0.057341579637820076072414341720623*x + 0.00086658848141621393294743877255483
// s =
//    0.533165284973017
//    0.691877113950473
//    0.315515631006063
//    0.686500927681584
//    0.834625671897373
//    0.018288277344192
//    0.750144314944967
//    0.988861088906495
BOOST_AUTO_TEST_CASE(solve_octic_8_real_roots)
{
	#ifndef NO_RANDOM
	srand((unsigned) time(NULL));
	#endif
	std::vector<std::complex<double> > solutions;
	std::vector<std::complex<double> > coefficients(9);

	double accuracy(1e-5);

	coefficients[8] = std::complex<double>(0.95788953015050193329216199344955);
	coefficients[7] = std::complex<double>(-4.616048869845871114437097477657);
	coefficients[6] = std::complex<double>(9.4118690334462455604524464083315);
	coefficients[5] = std::complex<double>(- 10.520090852574764399782523930068);
	coefficients[4] = std::complex<double>(6.9607710697687273625964490097213);
	coefficients[3] = std::complex<double>(- 2.7299923602832188539429250393216);
	coefficients[2] = std::complex<double>(0.59209079987668368905301009559429);
	coefficients[1] = std::complex<double>(- 0.057341579637820076072414341720623);
	coefficients[0] = std::complex<double>(0.00086658848141621393294743877255483);



	nups::solver::Octic<> octic_solver(accuracy);
	octic_solver.Solve(solutions, coefficients);
	//
	BOOST_CHECK(solutions.size()==8);

	std::vector<std::complex<double> > expected_solutions(8);

	expected_solutions[0] = std::complex<double>(0.533165284973017);
	expected_solutions[1] = std::complex<double>(0.691877113950473);
	expected_solutions[2] = std::complex<double>(0.315515631006063);
	expected_solutions[3] = std::complex<double>(0.686500927681584);
	expected_solutions[4] = std::complex<double>(0.834625671897373);
	expected_solutions[5] = std::complex<double>(0.018288277344192);
	expected_solutions[6] = std::complex<double>(0.750144314944967);
	expected_solutions[7] = std::complex<double>(0.988861088906495);

	// for (int jj=0; jj<8; jj++)
	// 	std::cout << solutions[jj] << "\t" << expected_solutions[jj] << std::endl;

	for (int jj=0; jj<8; jj++)
	{
		unsigned soln_counter(0);
		for (int ii=0; ii<8; ii++)
		{
			
			if (abs(solutions[ii]-expected_solutions[jj])<accuracy)
				soln_counter++;
		}
		BOOST_CHECK_EQUAL(soln_counter,1);
	}
}


// > (x-1)*(x-2)*(x-3)*(x-4)
// ans =
// (x - 1)*(x - 2)*(x - 3)*(x - 4)
// >> expand((x-1)*(x-2)*(x-3)*(x-4))
// ans =
// x^4 - 10*x^3 + 35*x^2 - 50*x + 24
BOOST_AUTO_TEST_CASE(solve_quartic_4_real_roots)
{
	std::vector<std::complex<double> > solutions;
	std::vector<std::complex<double> > coefficients(5);

	coefficients[4] = 1;
	coefficients[3] = -10;
	coefficients[2] = 35;
	coefficients[1] = -50;
	coefficients[0] = 24;
	nups::solver::Quartic quartic_solver;
	quartic_solver.Solve(solutions, coefficients);

	BOOST_CHECK(solutions.size()==4);

	std::vector<std::complex<double> > expected_solutions(4);

	expected_solutions[0] = std::complex<double>(1);
	expected_solutions[1] = std::complex<double>(2);
	expected_solutions[2] = std::complex<double>(3);
	expected_solutions[3] = std::complex<double>(4);

	for (int jj=0; jj<4; jj++)
	{
		unsigned soln_counter(0);
		for (int ii=0; ii<4; ii++)
		{
			if (abs(solutions[ii]-expected_solutions[jj])<1e-10)
				soln_counter++;
		}
		BOOST_CHECK_EQUAL(soln_counter,1);
	}
}







// BOOST_AUTO_TEST_CASE(solve_quartic_2_real_roots)
// {
// 	BOOST_CHECK(false && "implemented");
// }








// some matlab code for generating tests:
// a = rand_complex; s = zeros(4,1);
// for ii = 1:4
// s(ii) = rand_complex;
// a = a*(x-s(ii));
// end
// format long
// expand(vpa(a))
// s
//
//
//
// ans =
// x^4*(0.96304910599625548339730585212237 + 1.1715233183482309797796006023418i) + x^3*(0.53733279234624897505035665648996 + 1.7912943630496996702086268119871i) - x^2*(4.0115274206549127964284879992643 - 2.3549546567909716367121327400579i) - x*(2.9430655987982489554880706116498 + 0.83577829284967472428411276018548i) + (1.4515802519268858209061932953758 - 4.5689863685090177137745906415721i)
// s =
//  -1.358118174149949 + 1.089664271820386i
//   0.562412945988231 - 1.223429851379966i
//  -1.312303777711468 - 0.005634770866368i
//   0.970574191947381 - 0.336964720241924i
BOOST_AUTO_TEST_CASE(solve_quartic_0_real_roots)
{
	std::vector<std::complex<double> > solutions;
	std::vector<std::complex<double> > coefficients(5);

	coefficients[4] = std::complex<double>(0.96304910599625548339730585212237, 1.1715233183482309797796006023418);
	coefficients[3] = std::complex<double>(0.53733279234624897505035665648996,  1.7912943630496996702086268119871);
	coefficients[2] = std::complex<double>(-4.0115274206549127964284879992643, 2.3549546567909716367121327400579);
	coefficients[1] = std::complex<double>(-2.9430655987982489554880706116498,  -0.83577829284967472428411276018548);
	coefficients[0] = std::complex<double>(1.4515802519268858209061932953758, - 4.5689863685090177137745906415721);
	
	nups::solver::Quartic quartic_solver;
	quartic_solver.Solve(solutions, coefficients);


	BOOST_CHECK(solutions.size()==4);

	std::vector<std::complex<double> > expected_solutions(4);

	expected_solutions[0] = std::complex<double>(-1.358118174149949, 1.089664271820386);
	expected_solutions[1] = std::complex<double>(0.562412945988231, - 1.223429851379966);
	expected_solutions[2] = std::complex<double>(-1.312303777711468, - 0.005634770866368);
	expected_solutions[3] = std::complex<double>(0.970574191947381, - 0.336964720241924);

	for (int jj=0; jj<4; jj++)
	{
		unsigned soln_counter(0);
		for (int ii=0; ii<4; ii++)
		{
			if (abs(solutions[ii]-expected_solutions[jj])<1e-10)
				soln_counter++;
		}
		BOOST_CHECK_EQUAL(soln_counter,1);
	}

	std::complex<double> x(0.56243623425,-1.746127614234);


	std::complex<double> v4(- 4.4126034727066186051181830480481,9.1840269188548498886544028400651);
	BOOST_CHECK((abs(nups::solver::Quartic::EvaluatePolyNonMonic(x,coefficients)/v4-1.))<1e-14);

	std::complex<double> d4(- 29.396529641921596442672263599669, - 11.17737384922451863185169287606);
	BOOST_CHECK((abs(nups::solver::Quartic::EvaluateDerivNonMonic(x,coefficients)/d4-1.))<1e-14);

	// this actually tests degree 5 functionality...
	std::complex<double> v5(16.361532830119800589854309702301, 8.9193373061803859503872508096568);
	BOOST_CHECK((abs(nups::solver::Quartic::EvaluatePolyMonic(x,coefficients)/v5-1.))<1e-14);

	// this actually tests degree 5 functionality...
	std::complex<double> d5(- 11.350110773711822693095622755286, 42.496092442820758968461413713892);
	BOOST_CHECK((abs(nups::solver::Quartic::EvaluateDerivMonic(x,coefficients)/d5-1.))<1e-14);
}





BOOST_AUTO_TEST_CASE(solve_quadratic_2_real_roots)
{	
	std::vector<std::complex<double> > solutions;
	std::vector<std::complex<double> > coefficients(3);

	// this is the polynomial 2x^2+3x+1
	coefficients[2] = 2;
	coefficients[1] = 3;
	coefficients[0] = 1;
	
	nups::solver::Quadratic quadratic_solver;
	quadratic_solver.Solve(solutions, coefficients);


	BOOST_CHECK(solutions.size()==2);

	std::complex<double> a(-1), b(-0.5);
	unsigned soln_counter(0);

	for (int ii=0; ii<2; ii++)
	{
		if (abs(solutions[ii]-a)<1e-10)
			soln_counter++;
	}

	BOOST_CHECK_EQUAL(soln_counter,1);
	
	soln_counter = 0;

	for (int ii=0; ii<2; ii++)
	{
		if (abs(solutions[ii]-b)<1e-10)
			soln_counter++;
	}

	BOOST_CHECK_EQUAL(soln_counter,1);
}

BOOST_AUTO_TEST_CASE(solve_quadratic_0_real_roots)
{
	std::vector<std::complex<double> > solutions;
	std::vector<std::complex<double> > coefficients(3);

	// this is the polynomial 2x^2 + 3x + 3 
	coefficients[2] = 2;
	coefficients[1] = 3;
	coefficients[0] = 3;
	
	nups::solver::Quadratic quadratic_solver;
	quadratic_solver.Solve(solutions, coefficients);

	BOOST_CHECK(solutions.size()==2);

	std::complex<double> a(-0.75, -0.9682458365518542212948163499456), b(-0.75, +0.9682458365518542212948163499456);
	unsigned soln_counter(0);

	for (int ii=0; ii<2; ii++)
	{
		if (abs(solutions[ii]-a)<1e-10)
			soln_counter++;
	}

	BOOST_CHECK_EQUAL(soln_counter,1);
	
	soln_counter = 0;

	for (int ii=0; ii<2; ii++)
	{
		if (abs(solutions[ii]-b)<1e-10)
			soln_counter++;
	}

}


BOOST_AUTO_TEST_SUITE_END()



BOOST_AUTO_TEST_SUITE(NUPS_solvers)



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


