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
#ifndef TNP_OPS_COMP_HPP
#define TNP_OPS_COMP_HPP 1

#include <vector>
#include <tuple>

#include <boost/ptr_container/ptr_vector.hpp>

#include <tnp/ops/multiplication.hpp>
#include <tnp/polynomial.hpp>
#include <boost/math/special_functions/factorials.hpp>

namespace tnp {
  namespace ops {
    using namespace std;
    using namespace boost::ptr_container;

    /**
     * Efficient implementation of Composition based on pre-computed polynomials
     */
    class Composition { 
    public:
      const unsigned int order;

    private:
      /* contains (x <> .. <> x), the i-th entry is the order-th field of the i-th convolution */
      vector<StdPolynomial> convolutes;

      const vector<double> binomial;

      const vector<SumOfProducts> bell_polynomials;
      const vector<DerSumOfProducts> der_bell_polynomials;
      const Composition* last;

      const StdPolynomial& convolute(unsigned int n, unsigned int k);
      const StdPolynomial& getConvolute(unsigned int k);
      const StdPolynomial makeConvolute(unsigned int k);

      vector<SumOfProducts> compilePolynomials(const unsigned int order);

      vector<DerSumOfProducts> compileDerPolynomials(const unsigned int order);

    public:
      const vector<SumOfProducts>& bell() const { return bell_polynomials; }

      Composition(Composition* smaller) : order(smaller->order+1),
					  binomial(Multiplication::compileBinomial(smaller->order+1)), 
					  bell_polynomials(compilePolynomials(smaller->order+1)), 
					  der_bell_polynomials(compileDerPolynomials(smaller->order+1)), last(smaller)  {}

      Composition() : order(0), binomial(Multiplication::compileBinomial(0)), last(NULL) {}
      
      void apply(const vector<double>& a, const vector<double>& b,
		 vector<double>& target, unsigned int width) const;

      void apply(const double* f, const double* b,
		 double* target, unsigned int width) const;

    };

    class CompositionCache {
      
      static CompositionCache instance;

      vector<Composition*> cache;
    public:
      CompositionCache() {
	cache.push_back(new Composition());       
      }
      
      ~CompositionCache() {
	for (Composition* c : cache)
	  delete c;
      }

      inline static Composition* staticGetInstance(int order) {
	return instance.getInstance(order);
      };

      Composition* getInstance(int order) {
	if (cache.size() > order)
	  return cache[order];
	else {
	  while(cache.size() <= order)
	    cache.push_back(new Composition(cache.back()));
	  
	  return cache.back();
	}	  
      };
    };
  }
}
#endif
