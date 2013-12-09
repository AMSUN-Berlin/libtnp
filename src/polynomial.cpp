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

  using namespace std;
  using namespace boost;

  void addTerm(set<Term>& set, const Term& t) {
    if (t.factor == 0) 
      return;

    auto i = set.find(t);
    if (i != set.end()) {
      set.erase(t);
      Term t2(t);
      t2.factor++;
      set.insert(t2);
    } else
      set.insert(t);
  }

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

  Term Term::operator*(const Term& t) const {
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
      set<Term> ts;
      addTerm(ts, t);
      addTerm(ts, *this);

      return StdPolynomial(ts);
    }
  }    

  Term Term::operator/(const Term& t) const {
    Monomial div(monomial);

    for (auto e : t.monomial) {
      const unsigned int var = get<0>(e);
      const unsigned int pow = get<1>(e);
      
      auto contains = div.find(var);

      if (contains != div.end()) {
	if (div[var] > pow) {
	  div[var] -= pow;
	} else if (div[var] == pow) {
	  div.erase(var);
	} else {
	  return 0;
	}
      } else {
	return 0;
      }
    }

    return Term(factor / t.factor, div);
  }

  Term Term::operator%(const Term& t) const {
    Monomial div(monomial);

    for (auto e : t.monomial) {
      const unsigned int var = get<0>(e);
      const unsigned int pow = get<1>(e);
      
      auto contains = div.find(var);

      if (contains != div.end()) {
	if (div[var] < pow) {
	  return *this;
	}
      } else {
	return *this;
      }
    }
    
    return 0;
  }

  Term var(unsigned int idx) {
    return Term(1, Monomial({{idx, 1}}));
  }  

  StdPolynomial StdPolynomial::operator*(const StdPolynomial& f) const {
    set<Term> pTerms;
    
    for (const Term t : terms) 
      for (const Term t2 : f.terms) {
	addTerm(pTerms, t * t2);
      }
    
    return StdPolynomial(pTerms);
  }

  StdPolynomial StdPolynomial::operator+(const StdPolynomial& t) const {
    set<Term> sTerms(terms);
    for (Term t : t.terms)
      addTerm(sTerms, t);
    return StdPolynomial(sTerms);
  }

  StdPolynomial StdPolynomial::operator/(const Term& t) const {
    set<Term> dTerms;
    for (Term t2 : terms)
      addTerm(dTerms, t2 / t);
    return StdPolynomial(dTerms);
  }

  StdPolynomial StdPolynomial::operator%(const Term& t) const {
    set<Term> rTerms;
    for (Term t2 : terms)
      addTerm(rTerms, t2 % t);
    return StdPolynomial(rTerms);
  }
      
  double StdPolynomial::eval(const std::vector<double>& arg) const {
    double res = 0.0;
    for (Term t : terms)
      res += t.eval(arg);
    return res;
  }

  /* return the first variable, if any exists */
  optional<unsigned int> nextVarIn(const set<Term>& terms) {
    if (terms.size() > 0) {
      auto begin = terms.begin();
      if (begin->monomial.size() > 0)
	return optional<unsigned int>(get<0>(*(begin->monomial.begin())));
    }
    return none;
  }

  /**
   * Factorize an ordered list of terms
   * Returns a newly allocated Horner polynomial, if successful.
   */
  optional<HornerPolynomial*> HornerPolynomial::factorize(const set<Term>& terms) {
    std::cout << "Factorizing: " << StdPolynomial(terms) << std::endl;
    if (terms.size() == 0)
      return none;

    set<Term> ps;
    set<Term> qs;

    const optional<unsigned int> v = nextVarIn(terms);

    if (v) {
      std::cout << "Var: " << *v << std::endl;
      for (Term t : terms) {
	addTerm(ps, t / var(*v));
	addTerm(qs, t % var(*v));
      }
      HornerPolynomial* h = new HornerPolynomial(1, *v);

      std::cout << "== " << *v << " * " << StdPolynomial(ps) << " + " << StdPolynomial(qs) << std::endl;

      if (ps.size() > 0)
	h->hp = factorize(ps);

      if (qs.size() > 0)
	h->hq = factorize(qs);

      return optional<HornerPolynomial*>(h);
    }
    
    return none;
  }

  double HornerPolynomial::eval(const std::vector<double>& arg) const {
    double res = 0.0;

    if (hp)
      res += hp.get()->eval(arg);
    
    res *= factor * arg[variable];

    if (hq)
      res += hq.get()->eval(arg);

    return res;
  }

}
