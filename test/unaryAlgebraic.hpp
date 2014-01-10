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
#ifndef TNP_TEST_UNARY_ALGEBRAIC_ID_HPP
#define TNP_TEST_UNARY_ALGEBRAIC_ID_HPP 1

#include <tnp/npnumber.hpp>

#include <utility>

#include <boost/test/floating_point_comparison.hpp>

/**
 * These test cases take single NPNumbers and assert a 
 * algebraic identity upon a unary operation (e.g. adding a constant) 
 * on them.
 */
namespace tnp {
  namespace test {

    const unsigned int MANY_ITERATIONS = 2000000;

    NPNumber multiplyLoop(NPNumber result, unsigned int n);

    NPNumber bellLoop(NPNumber result, unsigned int n);

    NPNumber assignLoop(NPNumber result, unsigned int n);

    NPNumber squareLoop(NPNumber result, unsigned int n);

    void checkClose(const NPNumber& expected, const NPNumber& actual) {
      BOOST_CHECK_EQUAL(expected.params(), actual.params());
      BOOST_CHECK_EQUAL(expected.order(), actual.order());
      BOOST_CHECK_EQUAL(expected.data().size(), actual.data().size());

      for (int i = 0; i < expected.data().size(); ++i)
	if (expected.data()[i] == 0.0) {
	  BOOST_CHECK_SMALL(actual.data()[i], 1e-300);
	} else {
	  BOOST_CHECK_CLOSE(expected.data()[i], actual.data()[i], 1e-3);
	}
    }

    void testEquality(const NPNumber& in) {
      BOOST_CHECK_EQUAL(in, in);
    }

    void testAdditionWithZero(const NPNumber& in) {
      const NPNumber zero(in.params(), in.order(), 0.0);      
      BOOST_CHECK_EQUAL(in, in + zero);
    }

    void testMultiplicationWithOne(const NPNumber& in) {
      const NPNumber one(in.params(), in.order(), 1.0);      
      BOOST_CHECK_EQUAL(in * one, in);
    }

    void testManyMultiplicationsWithOne(const std::pair<unsigned int, unsigned int> sizes) {
      const NPNumber one(sizes.first, sizes.second, 1.0);
      NPNumber result = multiplyLoop(one, MANY_ITERATIONS);

      BOOST_CHECK_EQUAL(result, one);
    }

    void testManyAssignments(const std::pair<unsigned int, unsigned int> sizes) {
      const NPNumber zero(sizes.first, sizes.second, 0.0);
      cout << "Testing assignments for " << sizes << endl;
      NPNumber result = assignLoop(zero, MANY_ITERATIONS);

      BOOST_CHECK_EQUAL(result, zero);
    }

    void testManyBellsOfZero(const std::pair<unsigned int, unsigned int> sizes) {
      const NPNumber zero(sizes.first, sizes.second, 0.0);
      NPNumber result = bellLoop(zero, MANY_ITERATIONS);

      BOOST_CHECK_EQUAL(result, zero);
    }

    void testManySquaresOfOne(const std::pair<unsigned int, unsigned int> sizes) {
      const NPNumber one(sizes.first, sizes.second, 1.0);
      NPNumber result = squareLoop(one, MANY_ITERATIONS);

      BOOST_CHECK_EQUAL(result, one);
    }

    void testMultiplicationWithTwo(const NPNumber& in) {
      const NPNumber two(in.params(), in.order(), 2.0);      
      BOOST_CHECK_EQUAL(in * two, in + in);
    }

    void testMultiplicationWithConstantOne(const NPNumber& in) {
      BOOST_CHECK_EQUAL(in * 1, in);
    }

    void testMultiplicationWithConstantTwo(const NPNumber& in) {
      BOOST_CHECK_EQUAL(in * 2, in + in);
    }

    void testSelfSubtraction(const NPNumber& in) {
      const NPNumber zero(in.params(), in.order(), 0.0);      
      BOOST_CHECK_EQUAL(zero, in - in);
    }

    void testSelfDivision(const NPNumber& in) {
      const NPNumber one(in.params(), in.order(), 1.0);      
      BOOST_CHECK_EQUAL(one, in / in);
    }

    void testPowOne(const NPNumber& in) {
      BOOST_CHECK_EQUAL(in.pow(1), in);
    }

    void testPowZero(const NPNumber& in) {
      const NPNumber one(in.params(), in.order(), 1.0);      
      BOOST_CHECK_EQUAL(in.pow(0), one);
    }

    void testSquare(const NPNumber& in) {
      BOOST_CHECK_MESSAGE(in.pow(2) == in*in,
			  "\nExpected pow(2) (" << in << ") to be\n" <<
			  in * in << "\n" <<
			  "but got: " << in.pow(2) << "\n"
			  );
    }

    void testCube(const NPNumber& in) {
      BOOST_CHECK_EQUAL(in.pow(3), in*in*in);
    }
  
  }
}
#endif
