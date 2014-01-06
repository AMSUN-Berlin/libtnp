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
#ifndef TNP_TEST_UNARY_POLYNOMIAL_ID_HPP
#define TNP_TEST_UNARY_POLYNOMIAL_ID_HPP 1

#include <tnp/npnumber.hpp>

#include <utility>

#include <boost/test/floating_point_comparison.hpp>

namespace tnp {
  namespace test {

    using namespace boost;

    void testPolyPacking(const StdPolynomial& in) {
      HornerPolynomial h(in);
      PackedPolynomial packed(0);
      packInto(packed, &h);

      vector<double> args({13, 7});
      BOOST_CHECK_EQUAL(in.eval(args), eval(packed, args));
      BOOST_CHECK_MESSAGE(in.eval(args) == eval(packed, args), 
			  in << " != " << packed);
    }

    void testPolyFactorization(const StdPolynomial& in) {
      HornerPolynomial h(in);
      vector<double> args({42, 21});
      BOOST_CHECK_EQUAL(in.eval(args), h.eval(args));
    }

    void testPolyEquality(const StdPolynomial& in) {
      BOOST_CHECK_EQUAL(in, in);
      testPolyFactorization(in);
    }

    void testPolyAdditionWithZero(const StdPolynomial& in) {
      const StdPolynomial zero;      
      BOOST_CHECK_EQUAL(in, in + zero);
      testPolyFactorization(in + zero);
    }

    void testPolyMultiplicationWithOne(const StdPolynomial& in) {
      const StdPolynomial one(1);      
      BOOST_CHECK_EQUAL(in * one, in);
      testPolyFactorization(in * one);
    }

    void testPolyMultiplicationWithTwo(const StdPolynomial& in) {
      const StdPolynomial two(2);      
      BOOST_CHECK_EQUAL(in * two, in + in);
      testPolyFactorization(in * two);
      testPolyFactorization(in + in);
    }

    void testPolyMultiplicationWithConstantOne(const StdPolynomial& in) {
      BOOST_CHECK_EQUAL(in * Term(1), in);
    }

    void testPolyMultiplicationWithConstantTwo(const StdPolynomial& in) {
      BOOST_CHECK_EQUAL(in * Term(2), in + in);
    }
 
  }
}

#endif
