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
  
  /*  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testManyMultiplicationsWithOne, testDimensions.begin(), testDimensions.end() ) );
  */
  const std::vector<UnaryAnalyticFunctionTest>& analyticTestCases = UnaryAnalyticFunctionTest::testCases();

  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testUnaryAnalyticFunction, analyticTestCases.begin(), analyticTestCases.end() ) );

  PolynomialTestSuite* polynomialSuite = new PolynomialTestSuite();
  framework::master_test_suite().add( polynomialSuite );

  PolynomialFunctionTestSuite* polynomialFunctionSuite = new PolynomialFunctionTestSuite();
  framework::master_test_suite().add( polynomialFunctionSuite );
  
  return 0;
}

