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

#include <vector>

#include <boost/test/floating_point_comparison.hpp>

/**
 * These test cases take single NPNumbers and assert a 
 * algebraic identity upon a unary operation (e.g. adding a constant) 
 * on them.
 */
namespace tnp {
  namespace test {

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
      BOOST_CHECK_EQUAL(in, in * one);
    }

    void testSelfSubtraction(const NPNumber& in) {
      const NPNumber zero(in.params(), in.order(), 0.0);      
      BOOST_CHECK_EQUAL(zero, in - in);
    }

    void testSelfDivision(const NPNumber& in) {
      const NPNumber one(in.params(), in.order(), 1.0);      
      BOOST_CHECK_EQUAL(one, in / in);
    }
  
  }
}
#endif
