/*
 * Copyright (C) 2012 uebb.tu-berlin.de.
 *
 * This file is part of tnp
 *
 * tnp is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * tnp is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with tnp. If not, see <http://www.gnu.org/licenses/>.
 */

#include <tnp/polynomial.hpp>

#include <boost/test/included/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "polynomialIdentities.hpp"

/**
 * test cases involving polynomial evaluation
 */
namespace tnp {
  namespace test {

    using namespace std;
    using namespace boost::unit_test;

    class PolynomialFixture {
    public:
      StdPolynomial stdPoly;
    
      virtual ~PolynomialFixture() {}

      PolynomialFixture(StdPolynomial p) : stdPoly(p) {}

      virtual double eval(const vector<double>& args) const = 0;
    };

    class Polynomial1 : public PolynomialFixture {
      /* x_1²*x_2³ + x_2 */
    public:
      Polynomial1() : PolynomialFixture(StdPolynomial(
						      (var(0)^2)*(var(1)^3) + var(1))) {}
    
      virtual double eval(const vector<double>& args) const {
	return args[0] * args[0] * args[1]*args[1]*args[1] + args[1];
      }
    };

    class PolynomialTest {

      static vector<StdPolynomial> makeStdPolys() {
	vector <StdPolynomial> res;
	/* x_1²*x_2³ + x_2 */
	res.push_back(StdPolynomial(var(0)^2)*(var(1)^3) + var(1));
	/* x_0² + 2 */
	res.push_back(StdPolynomial(var(0)^2) + 2);		      
	return res;
      }

    public:
      static const vector<StdPolynomial>& stdPolys() {
	static vector <StdPolynomial> res(makeStdPolys());
	return res;
      }

      PolynomialTest(vector<double> a, const PolynomialFixture* f) : args(a), fixture(f) {}

      ~PolynomialTest() {
	delete fixture;
      }

      const vector<double> args;
      const PolynomialFixture* const fixture;      

    };

    class PolynomialTestSuite : public test_suite {
      static vector<PolynomialFixture*> fixtures() {
	PolynomialFixture* p1 = new Polynomial1();
      
	vector<PolynomialFixture*> fs({p1});
	return fs;
      }

      static vector<vector<double>> inputs() {
	vector<vector<double>> res;
	res.push_back(vector<double>{1,1,1,1,1,1});
	res.push_back(vector<double>{0,1,2,3,4,5});
	res.push_back(vector<double>{0.1,1.1,.2,.3,.4,.5});
	return res;
      }

      static vector<PolynomialTest*> mkTestCases() {
	vector<PolynomialTest*> res;
	for (vector<double> i : inputs()) {
	  for (PolynomialFixture* f : fixtures()) {
	    res.push_back(new PolynomialTest(i, f));
	  }
	}
	return res;
      }

      static void testPolynomial(const PolynomialTest* const test) {
	BOOST_CHECK_EQUAL(test->fixture->stdPoly.eval(test->args), test->fixture->eval(test->args));
      }

      static void testHornerPolynomial(const PolynomialTest* const test) {
	HornerPolynomial h(test->fixture->stdPoly);
	const double r1 = h.eval(test->args);
	const double r2 = test->fixture->eval(test->args);
	BOOST_CHECK_MESSAGE(r1 == r2, h << " != " << test->fixture->stdPoly << " because " << r1 << " != " << r2 << " at (" << test->args[0] << ", " << test->args[1] << ")");
      }

    public:
      vector<PolynomialTest*> testCases;

      PolynomialTestSuite() : test_suite("Polynomials"), testCases(mkTestCases()) {
	add( BOOST_PARAM_TEST_CASE( &testPolynomial, testCases.begin(), testCases.end() ) );

	add( BOOST_PARAM_TEST_CASE( &testHornerPolynomial, testCases.begin(), testCases.end() ) );

	const std::vector<tnp::StdPolynomial> testPolys = PolynomialTest::stdPolys();

	add( BOOST_PARAM_TEST_CASE( &testPolyEquality, testPolys.begin(), testPolys.end() ) );
  
	add( BOOST_PARAM_TEST_CASE( &testPolyAdditionWithZero, testPolys.begin(), testPolys.end() ) );
  
	add( BOOST_PARAM_TEST_CASE( &testPolyMultiplicationWithConstantOne, testPolys.begin(), testPolys.end() ) );

	add( BOOST_PARAM_TEST_CASE( &testPolyMultiplicationWithConstantTwo, testPolys.begin(), testPolys.end() ) );

	add( BOOST_PARAM_TEST_CASE( &testPolyMultiplicationWithOne, testPolys.begin(), testPolys.end() ) );

	add( BOOST_PARAM_TEST_CASE( &testPolyMultiplicationWithTwo, testPolys.begin(), testPolys.end() ) );
      }

      ~PolynomialTestSuite() {
	for (PolynomialTest* t : testCases)
	  delete t;
      }
    };

  }
}
