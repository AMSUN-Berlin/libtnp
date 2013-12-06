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

#include <math.h>

namespace tnp {

  using namespace boost::heap;

  bool Term::operator <(const Term& o) const {
    return o.monomial < monomial;
  }

  double Term::eval(const std::vector<double>& arg) const {
    double res = factor;
    for(auto e : monomial) {
      res *= pow(arg[std::get<0>(e)], std::get<1>(e));
    }
    return res;
  }

  Term Term::operator^(unsigned int p) const {
    if (p == 0)
      return Term(1);

    Monomial m2(monomial);
    for (auto e : m2)
      m2[std::get<0>(e)] += (p-1);
    return Term(pow(factor, p), m2);
  }

  Term Term::operator*(const int f) const {
    return Term(factor*f, monomial);
  }

  Term Term::operator*(const Term t) const {
    Monomial m(monomial);
    for (auto e : t.monomial) {
      unsigned int key = std::get<0>(e);
      if (m.find(key) != m.end()) {
	m[key] += t.monomial.at(key);
      } else
	m[key] = t.monomial.at(key);
    }
    return Term(factor * t.factor, m);
  }

  StdPolynomial Term::operator+(const Term& t) const {
    if (t.monomial == monomial)
      return StdPolynomial(Term(factor + t.factor, monomial));
    else {
      priority_queue<Term> ts;
      ts.push(t);
      ts.push(*this);
      return StdPolynomial(ts);
    }
  }    

  Term var(unsigned int idx) {
    return Term(1, Monomial({{idx, 1}}));
  }
      
  double StdPolynomial::eval(const std::vector<double>& arg) const {
    double res = 0.0;
    for (Term t : terms)
      res += t.eval(arg);
    return res;
  }
}
