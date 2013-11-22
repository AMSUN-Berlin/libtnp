#include <boost/test/included/unit_test.hpp>
#include <boost/test/parameterized_test.hpp>

#include "unaryAlgebraic.hpp"
#include "numberGenerator.hpp"

using namespace boost::unit_test;

using namespace tnp::test;

test_suite*
init_unit_test_suite( int argc, char* argv[] ) {

  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testEquality, testNumbers.begin(), testNumbers.end() ) );

  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testSelfSubtraction, testNumbers.begin(), testNumbers.end() ) );

 
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testAdditionWithZero, testNumbers.begin(), testNumbers.end() ) );

  /*
  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testMultiplicationWithOne, testNumbers.begin(), testNumbers.end() ) );

  framework::master_test_suite().
    add( BOOST_PARAM_TEST_CASE( &testSelfDivision, testNumbers.begin(), testNumbers.end() ) );
  */

  return 0;
}

