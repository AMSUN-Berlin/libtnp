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

#ifndef NP_POLYNOMIAL_HPP
#define NP_POLYNOMIAL_HPP 1

#include <tnp/polynomial.hpp>
#include <tnp/npnumber.hpp>

#include <vector>
#include <utility>

namespace tnp {

  using namespace std;

  class NPPolynomial {

    NPNumber eval(const Term& t, const vector<NPNumber>& arg) const {
      NPNumber res(arg[0].params()+1, constant(t.factor, arg[0].data().size()));
      for (auto e : t.monomial) {
	res *= arg[get<0>(e)].pow(get<1>(e));
      }
      return res;
    }

  public:
    StdPolynomial poly;

    NPPolynomial(const StdPolynomial& p) : poly(p) {}

    NPNumber eval(const vector<NPNumber>& arg) const {
      NPNumber res(arg[0].params()+1, constant(0.0, arg[0].data().size()));
      for (Term t : poly.terms)
	res += eval(t, arg);
      return res;
    }

  };  

}

#endif
