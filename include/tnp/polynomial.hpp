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

#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP 1

#include <boost/heap/priority_queue.hpp>
#include <utility>
#include <map>

namespace tnp {

  using namespace boost::heap;

  class StdPolynomial;

  /**
   * standard representation of polynomials
   */
  typedef std::map<unsigned int, unsigned int> Monomial;

  class Term {
  public:
    int factor;
    Monomial monomial;

    Term(int f) : factor(f) {}

    Term(int f, Monomial m) : factor(f), monomial(m) {}

    inline double eval(const std::vector<double>& arg) const;

    Term operator ^(unsigned int p) const;

    bool operator <(const Term& o) const;
    
    Term operator*(const int f) const;

    Term operator*(const Term f) const;

    StdPolynomial operator+(const Term& t) const; 
  };

  Term var(unsigned int idx);

  class StdPolynomial {
  public:
    priority_queue<Term> terms;

    StdPolynomial(Term t) { terms.push(t); }

    StdPolynomial(priority_queue<Term> t) : terms(t) {}

    double eval(const std::vector<double>& arg) const;
  };
}

#endif
