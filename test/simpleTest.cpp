#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include "polynomialTest.hpp"
#include "polynomialIdentities.hpp"
#include "unaryAlgebraic.hpp"
#include "analyticFunctions.hpp"
#include "numberGenerator.hpp"


using namespace boost::unit_test;

using namespace tnp::test;

test_suite*
init_unit_test_suite( int argc, char* argv[] ) {
 
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testEquality, testNumbers().begin(), testNumbers().end() ) );

  
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testSelfSubtraction, testNumbers().begin(), testNumbers().end() ) );

  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testAdditionWithZero, testNumbers().begin(), testNumbers().end() ) );

  
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testMultiplicationWithConstantOne, testNumbers().begin(), testNumbers().end() ) );

  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testMultiplicationWithConstantTwo, testNumbers().begin(), testNumbers().end() ) );

  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testMultiplicationWithOne, testNumbers().begin(), testNumbers().end() ) );

  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testMultiplicationWithTwo, testNumbers().begin(), testNumbers().end() ) );
  
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testManyMultiplicationsWithOne, testDimensions.begin(), testDimensions.end() ) );

  const std::vector<UnaryAnalyticFunctionTest>& analyticTestCases = UnaryAnalyticFunctionTest::testCases();

  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testUnaryAnalyticFunction, analyticTestCases.begin(), analyticTestCases.end() ) );

  const std::vector<PolynomialTest*>& polynomialTestCases = PolynomialTest::testCases();

  /*framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testPolynomial, polynomialTestCases.begin(), polynomialTestCases.end() ) );*/

  const std::vector<tnp::StdPolynomial> testPolys = PolynomialTest::stdPolys();

  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testPolyEquality, testPolys.begin(), testPolys.end() ) );
  
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testPolyAdditionWithZero, testPolys.begin(), testPolys.end() ) );
  
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testPolyMultiplicationWithConstantOne, testPolys.begin(), testPolys.end() ) );

  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testPolyMultiplicationWithConstantTwo, testPolys.begin(), testPolys.end() ) );

  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testPolyMultiplicationWithOne, testPolys.begin(), testPolys.end() ) );

  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testPolyMultiplicationWithTwo, testPolys.begin(), testPolys.end() ) );
  
  return 0;
}

