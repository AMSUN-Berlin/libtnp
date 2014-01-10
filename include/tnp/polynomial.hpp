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

#include <iostream>
#include <boost/optional.hpp>
#include <utility>
#include <map>
#include <set>

#include "prettyprint.hpp"

namespace tnp {

  using namespace boost;
  using namespace std;

  inline double powi (double base, unsigned int exp) {
    double res = 1;
    while (exp) {
        if (exp & 1)
            res *= base;
        exp >>= 1;
        base *= base;
    }
    return res;
  }


  class StdPolynomial;

  /**
   * standard representation of polynomials
   */
  typedef map<unsigned int, unsigned int> Monomial;

  class Term {
    Term deriveTotal(unsigned int var, unsigned int width) const;
  public:
    Monomial monomial;
    int factor;

    Term() : monomial(), factor(0) {}

    Term(const Term& o) : monomial(o.monomial), factor(o.factor) {}

    Term(const int f) : monomial(), factor(f) {}

    Term(int f, Monomial m) : monomial(m), factor(f) {}

    double eval(const vector<double>& arg) const;

    double eval(const vector<double>& arg, const unsigned int width) const;

    Term operator ^(unsigned int p) const;

    bool operator <(const Term& o) const;
    
    Term operator*(const int f) const;

    Term operator*(const Term& f) const;

    Term operator/(const Term& term) const;

    Term operator%(const Term& term) const;

    StdPolynomial operator+(const Term& t) const; 

    bool operator==(const Term& o) const { return factor == o.factor && monomial == o.monomial; }

    Term partialDerivative(unsigned int var) const;

    set<Term> totalDerivative(unsigned int width) const;

    unsigned int variables() const { return monomial.size() > 0 ? monomial.rbegin()->first + 1 : 0; }

    friend std::ostream& operator<<(std::ostream& out, const Term& t) {
      out << t.factor << " * (";
      for (auto e : t.monomial)
	out << "x_" << get<0>(e) << "^" << get<1>(e) << " ";
      out << ")";
      return out;
    }
  };

  Term var(unsigned int idx);

  void addTerm(set<Term>& set, const Term& t);

  struct AddDelims { static const pretty_print::delimiters_values<char> values; }; 

  class StdPolynomial {
  public:
    set<Term> terms;

    StdPolynomial() : terms() {}
    
    StdPolynomial(const StdPolynomial& o) : terms(o.terms) {}

    StdPolynomial(const Term& t) : terms() { addTerm(terms, t); }

    StdPolynomial(const set<Term>& t) : terms(t) {}

    StdPolynomial(const unsigned int i) : terms() { addTerm(terms, Term(i)); }

    unsigned int variables() const { return terms.empty() ? 0 : terms.begin()->variables(); }

    double eval(const vector<double>& arg) const;

    double eval(const vector<double>& arg, const unsigned int width) const;

    StdPolynomial operator*(const StdPolynomial& f) const;

    StdPolynomial operator/(const Term& term) const;

    StdPolynomial operator%(const Term& term) const;

    StdPolynomial operator+(const StdPolynomial& t) const; 

    bool operator==(const StdPolynomial& o) const { return terms == o.terms; }

    StdPolynomial& operator+=(const StdPolynomial& o);

    StdPolynomial partialDerivative(const unsigned int var) const;

    StdPolynomial totalDerivative(unsigned int width) const;
    
    friend std::ostream& operator<<(std::ostream& out, const StdPolynomial& p) {
      out << p.terms; //pretty_print::custom_delims<AddDelims>(p.terms);
      return out;
    }
  };

  struct Product {
    int factor;
    vector<int> fields;
  };

  class SumOfProducts {
  public:
    static long lookups;
    static long evals;
    vector<Product> sum;

    SumOfProducts(const StdPolynomial& poly) {      
      
      for (const Term& t : poly.terms) {
	Product p;
	p.factor = t.factor;
	if (p.factor != 0)
	  for (auto pair : t.monomial) {
	    const int field = get<0>(pair);
	    const int  power = get<1>(pair);

	    for (int i = 0; i < power; ++i)
	      p.fields.push_back(field);
	  }
	sum.push_back(p);
      }
    }

    inline double eval(const double* arg, const int width) const {
      double res = 0.0;
      for (const Product& p : sum) {
	evals++;
	double prod = p.factor;
	for (int f : p.fields) {
	  prod *= arg[f*width];
	}
	res += prod;
      }
      return res;
    };
  };

  struct DerProductField {
    bool is_der;
    int key;
  };

  struct DerProduct {
    int factor;
    vector<DerProductField> fields;
  };

  class DerSumOfProducts {
  public:
    vector<DerProduct> sum;

    DerSumOfProducts(const StdPolynomial& poly) {
      const int vars = poly.variables();
      /* derive poly, mark maximum var as derivative */
      const StdPolynomial& derP = poly.totalDerivative(vars);

      for (const Term& t : derP.terms) {
	DerProduct p;
	p.factor = t.factor;
	if (p.factor != 0)
	  for (auto pair : t.monomial) {
	    const int var = get<0>(pair);
	    const int power = get<1>(pair);
	    
	    DerProductField field;
	    
	    if (var > vars) {	      
	      field.key = var - vars;
	      field.is_der = true;
	    } else {
	      field.key = var;
	      field.is_der = false;	      
	    }

	    for(int i = 0; i < power; ++i)
	      p.fields.push_back(field);
	  }
	sum.push_back(p);
      }
    }

    inline double eval(const double* arg, const int width, const int der) const {
      double res = 0.0;
      for (const DerProduct& p : sum) {
	SumOfProducts::evals++;
	double prod = p.factor;
	for (DerProductField f : p.fields) {
	  prod *= arg[f.key*width + (f.is_der ? der : 0)];
	}
	res += prod;
      }
      return res;
    };
  };
 
  /**
   * f * var^power * (p_1 + .. + p_n) + q_1 + .. + q_m
   */
  class HornerPolynomial {
    /**
     * Factorize an ordered list of terms
     * Returns a newly allocated Horner polynomial, if successful.
     */
    optional<HornerPolynomial*> factorize(const set<Term>& terms);

  public:
    unsigned int variable;
    unsigned int power;
    int factor;
    optional<HornerPolynomial*> hp;
    optional<HornerPolynomial*> hq;

    HornerPolynomial(const StdPolynomial& p) {
      optional<HornerPolynomial*> ho = factorize(p.terms);

      if (ho) {
	HornerPolynomial* h = *ho;
	variable = h->variable;
	power = h->power;
	factor = h->factor;
	hp = h->hp;
	hq = h->hq;
	
	h->hp = none;
	h->hq = none;
	delete h;	
      }
    }

    HornerPolynomial(const HornerPolynomial& o) : variable(o.variable), power(o.power), factor(o.factor), 
						  hp(o.hp ? optional<HornerPolynomial*>(new HornerPolynomial(**o.hp)) : none), 
						  hq(o.hq ? optional<HornerPolynomial*>(new HornerPolynomial(**o.hq)) : none) {}

    HornerPolynomial() : variable(0), power(0), factor(0) {}
    
    HornerPolynomial(int f) : variable(0), power(0), factor(f) {}

    HornerPolynomial(int f, unsigned int v) : variable(v), power(1), factor(f) {}

    HornerPolynomial(int f, unsigned int v, unsigned int pwr) : variable(v), power(pwr), factor(f) {}

    HornerPolynomial(int f, unsigned int v, unsigned int pwr, HornerPolynomial p,  
		     HornerPolynomial q) :
      variable(v), power(pwr), factor(f), hp(new HornerPolynomial(p)), hq(new HornerPolynomial(q)) {}    

    ~HornerPolynomial() {
      if (hp)
	delete *hp;
      if (hq)
	delete *hq;
    }

    double eval(const vector<double>& arg) const;    

    double eval(const vector<double>& arg, const unsigned int width) const;    

    double evalDer(const vector<double>& arg, const unsigned int der, const unsigned int width) const;

    HornerPolynomial* der(const unsigned int der) const;

    friend std::ostream& operator<<(std::ostream& out, const HornerPolynomial& p) {
      out << p.factor << " x_" << p.variable << "^" << p.power;
      if (p.hp)
	out << " * (" << (**(p.hp)) << ")";
      if (p.hq)
	out << " + " << (**(p.hq));
      return out;
    }
  };

  enum PCode { CONST, LOAD, DER, MULT, ADD};

  struct PInst {
    PCode code;
    int carry_1;
    int carry_2;
  };  

  ostream& operator<<(ostream& out, const PInst& e);

  typedef vector<PInst> PackedPolynomial;

  void packInto(PackedPolynomial& packed, const HornerPolynomial* p);

  void packDerInto(PackedPolynomial& der, const HornerPolynomial* p);

  double eval(const PackedPolynomial& p, const vector<double>& arg);

  double eval(const PackedPolynomial& p, const vector<double>& arg, const unsigned int width);

  double eval(const PackedPolynomial& p, const vector<double>& arg, const unsigned int der, const unsigned int width);

}

#endif
