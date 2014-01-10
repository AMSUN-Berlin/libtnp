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

  long SumOfProducts::lookups = 0; 
  long SumOfProducts::evals = 0; 
  const pretty_print::delimiters_values<char> AddDelims::values = { "", " + ", "" };

  void addTerm(set<Term>& set, const Term& t) {
    if (t.factor == 0) 
      return;

    auto i = set.find(t);
    if (i != set.end()) {
      const unsigned int factor = i->factor;
      Term t2(t);
      set.erase(t);
      t2.factor += (factor);
      set.insert(t2);
    } else
      set.insert(t);
  }

  bool Term::operator <(const Term& o) const {
    return o.monomial < monomial;
  }

  double Term::eval(const std::vector<double>& arg) const {
    return eval(arg, 1);
  }

  double Term::eval(const std::vector<double>& arg, const unsigned int width) const {
    double res = factor;
    for(auto const& e : monomial) {
      res *= powi(arg[std::get<0>(e) * width], std::get<1>(e));
    }
    return res;
  }

  Term Term::operator^(unsigned int p) const {
    if (p == 0)
      return Term(1);
    
    Monomial m2(monomial);
    for (auto const& e : m2)
      m2[std::get<0>(e)] += (p-1);
    return Term(pow(factor, p), m2);
  }

  Term Term::operator*(const int f) const {
    return Term(factor*f, monomial);
  }

  Term Term::operator*(const Term& t) const {
    Monomial m(monomial);
    for (auto const& e : t.monomial) {
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

    for (auto const& e : t.monomial) {
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

    for (auto const& e : t.monomial) {
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

  Term Term::partialDerivative(unsigned int var) const {
    auto x = monomial.find(var);
    if (x == monomial.end()) {
      return Term(); // == zero
    } else {
      Monomial der(monomial);
      const unsigned int deg = get<1>(*x);
      if (deg == 1) {
	der.erase(var);
	return Term(factor, der);
      } else {
	der[var] = deg - 1;
	return Term(factor * deg, der);
      }      
    }
  }

  Term Term::deriveTotal(unsigned int var, unsigned int width) const {
    Monomial der(monomial);
    unsigned int f = 1;
    if (der[var] == 1) {
      der.erase(var);
    } else {
      f = der[var];
      der[var] -= 1;
    }

    if (monomial.find(var + width) == monomial.end()) 
      der[var + width] = 1;
    else
      der[var + width] += 1;
    
    return Term(factor*f, der);
  }

  set<Term> Term::totalDerivative(unsigned int width) const {
    set<Term> dTerms;

    /* derive for every key */
    for (auto const& e : monomial) {
      addTerm(dTerms, deriveTotal(e.first, width));
    }

    return dTerms;
  }

  Term var(unsigned int idx) {
    return Term(1, Monomial({{idx, 1}}));
  }  

  StdPolynomial StdPolynomial::operator*(const StdPolynomial& f) const {
    set<Term> pTerms;
    
    for (const Term& t : terms) 
      for (const Term& t2 : f.terms) {
	addTerm(pTerms, t * t2);
      }
    
    return StdPolynomial(pTerms);
  }

  StdPolynomial& StdPolynomial::operator+=(const StdPolynomial& o) {
    for (Term t : o.terms)
      addTerm(terms, t);    
    return *this;
  }

  StdPolynomial StdPolynomial::operator+(const StdPolynomial& t) const {
    StdPolynomial p(terms);
    p += t;
    return p;
  }

  StdPolynomial StdPolynomial::operator/(const Term& t) const {
    set<Term> dTerms;
    for (const Term& t2 : terms)
      addTerm(dTerms, t2 / t);
    return StdPolynomial(dTerms);
  }

  StdPolynomial StdPolynomial::operator%(const Term& t) const {
    set<Term> rTerms;
    for (const Term& t2 : terms)
      addTerm(rTerms, t2 % t);
    return StdPolynomial(rTerms);
  }
      
  double StdPolynomial::eval(const std::vector<double>& arg) const {
    return eval(arg, 1);
  }

  double StdPolynomial::eval(const std::vector<double>& arg, const unsigned int width) const {
    double res = 0.0;
    for (const Term& t : terms)
      res += t.eval(arg, width);
    return res;
  }

  StdPolynomial StdPolynomial::partialDerivative(const unsigned int var) const {
    set<Term> dterms;
    for (const Term& t : terms)
      addTerm(dterms, t.partialDerivative(var));
    return dterms;
  }

  StdPolynomial StdPolynomial::totalDerivative(const unsigned int width) const {
    set<Term> dterms;
    for (const Term& t : terms)
      for (const Term& dt : t.totalDerivative(width))
	addTerm(dterms, dt);
    return dterms;
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

  unsigned int maxCommonPwr(const unsigned int var, const set<Term>& sum) {
    unsigned int pwr = 0;
    for (const Term& t : sum)
      if (t.monomial.find(var) != t.monomial.end())
	if (t.monomial.at(var) < pwr || pwr == 0)
	  pwr = t.monomial.at(var);
    return pwr;
  }

  /**
   * Factorize an ordered list of terms
   * Returns a newly allocated Horner polynomial, if successful.
   */
  optional<HornerPolynomial*> HornerPolynomial::factorize(const set<Term>& terms) {
    if (terms.size() == 0)
      return none;

    set<Term> ps;
    set<Term> qs;

    const optional<unsigned int> v = nextVarIn(terms);

    if (v) {
      const unsigned int pwr = maxCommonPwr(*v, terms);
      for (const Term& t : terms) {
	addTerm(ps, t / (var(*v)^pwr));
	addTerm(qs, t % (var(*v)^pwr));
      }
      HornerPolynomial* h = new HornerPolynomial(1, *v, pwr);

      h->hp = factorize(ps);

      h->hq = factorize(qs);
      
      return optional<HornerPolynomial*>(h);
    }
    
    if (terms.size() > 0) {
      int sum = 0;
      for (const Term& t : terms)
	sum += t.factor;
        
      return optional<HornerPolynomial*>(new HornerPolynomial(sum));
    }
    return none;
  }

  double HornerPolynomial::eval(const std::vector<double>& arg) const {
    return eval(arg, 1);
  }

  double HornerPolynomial::eval(const std::vector<double>& arg, const unsigned int width) const {
    double res = 0.0;

    if (hp) {
      res += hp.get()->eval(arg, width);
      res *= factor * powi(arg[variable*width], power);
    } else {
      res = factor * powi(arg[variable*width], power);
    }
    
    if (hq)
      res += hq.get()->eval(arg, width);

    return res;
  }

  double HornerPolynomial::evalDer(const vector<double>& arg, const unsigned int der, const unsigned int width) const {
    double res = 0.0;
    if (hp) {      
      res += hp.get()->eval(arg, width); // p
      res *= factor * arg[variable*width+der]; // d g_i / dx
      if (power > 1) 
      	res *= power * powi(arg[variable*width], power - 1); // d ( f*g_i^p ) / dx 
      
    
      double sres = hp.get()->evalDer(arg, der, width);
      sres *= factor * powi(arg[variable*width], power);
      res += sres;
    }

    if (hq)
      res += hq.get()->evalDer(arg, der, width);

    return res;
  };

  ostream& operator<<(ostream& out, const PInst& inst) {
    switch (inst.code) {
    case CONST : out << "CONST " << inst.carry_1; break;
    case LOAD : out << "LOAD {var=" << inst.carry_1 << ", power=" << inst.carry_2 << "}"; break;
    case DER : out << "DER {var=" << inst.carry_1 << "}"; break;
    case MULT : out << "MULT"; break;
    case ADD : out << "ADD"; break;
    }
    return out;
  }

  void packInto(PackedPolynomial& packed, const HornerPolynomial* p) {    
    packed.push_back(PInst());
    if (p->power > 0) {
      packed.back().carry_1 = p->variable;
      packed.back().carry_2 = p->power;
      packed.back().code = LOAD;

      if (p->factor != 1) {
	packed.push_back(PInst());
	packed.back().code = MULT;
      }      
    } else {
      packed.back().carry_1 = p->factor;
      packed.back().code = CONST;
    }

    if (p->hp) {
      packInto(packed, *(p->hp));      
      packed.push_back(PInst());
      packed.back().code = MULT;
    }

    if (p->hq) {
      packInto(packed, *(p->hq)); 
      packed.push_back(PInst());
      packed.back().code = ADD;
    }
  }

  void packDerInto(PackedPolynomial& packed, const HornerPolynomial* p) {    
    packed.push_back(PInst());
    if (p->power > 0) {
      /* f*k*x^(k-1)*x'*(p) + f*x^k*p' + q' */

      packed.back().carry_1 = p->variable; // x'
      packed.back().code = DER;

      // f * k * x^(k-1)
      if (p->power > 1) {
	packed.push_back(PInst());
	packed.back().carry_1 = p->variable;
	packed.back().carry_2 = p->power - 1;
	packed.back().code = LOAD;

	packed.push_back(PInst());
	packed.back().code = MULT;

	packed.push_back(PInst());
	packed.back().code = CONST;
	packed.back().carry_1 = p->power * p->factor;

	packed.push_back(PInst());
	packed.back().code = MULT;
      } else if (p->factor != 1) {
	packed.push_back(PInst());
	packed.back().code = CONST;
	packed.back().carry_1 = p->power * p->factor;

	packed.push_back(PInst());
	packed.back().code = MULT;
      }

      // p
      if (p->hp) {
	packInto(packed, *(p->hp));      
	packed.push_back(PInst());
	packed.back().code = MULT;    

	packed.push_back(PInst());

	// f*x^k
	if (p->power > 0) {
	  packed.back().carry_1 = p->variable;
	  packed.back().carry_2 = p->power;
	  packed.back().code = LOAD;

	  if (p->factor != 1) {
	    packed.push_back(PInst());
	    packed.back().code = MULT;
	  }      
	} else {
	  packed.back().carry_1 = p->factor;
	  packed.back().code = CONST;
	}

	//p'
	packDerInto(packed, *(p->hp));
	packed.push_back(PInst());
	packed.back().code = MULT;  

	packed.push_back(PInst());
	packed.back().code = ADD;
      }

      //q'
      if (p->hq) {
	packDerInto(packed, *(p->hq));
	packed.push_back(PInst());
	packed.back().code = ADD;	
      }
    } else {
      packed.back().carry_1 = 0;
      packed.back().code = CONST;
    }
  }

  inline double eval(const PackedPolynomial& packed, const vector<double>& arg, const unsigned int der, const unsigned int width) {
    vector<double> data;
    data.reserve(packed.size());

    for (const PInst& inst : packed) {
      switch (inst.code) {
      case CONST : data.push_back(inst.carry_1); break;
      case LOAD : {
	data.push_back(powi(arg.at(inst.carry_1 * width), inst.carry_2)); break;
      }
      case DER : data.push_back(arg.at(inst.carry_1 * width + der)); break;
      case MULT : {
	double last = data.back();
	data.pop_back();
	data.back() *= last;
	break;
      }
      case ADD : {
	double last = data.back(); 
	data.pop_back();
	data.back() += last;
	break;
      }
      }
    }
    
    return data.back();
  }

  double eval(const PackedPolynomial& p, const vector<double>& arg) {
    return eval(p, arg, 0, 1);
  }

  double eval(const PackedPolynomial& p, const vector<double>& arg, const unsigned int width) {
    return eval(p, arg, 0, width);
  }
}
