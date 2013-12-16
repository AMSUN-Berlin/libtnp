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

      const unsigned int order;

      /* contains (x <> .. <> x), the i-th entry is the order-th field of the i-th convolution */
      vector<StdPolynomial> convolutes;

      const vector<double> binomial;

      const ptr_vector<HornerPolynomial> bell_polynomials;    

      const StdPolynomial& convolute(unsigned int k, unsigned int n) {
	if (n < order) {
	  Composition& comp = ensureExistance(n);
	  return comp.getConvolute(k);
	} else {
	  return getConvolute(k);
	}
      }

      const StdPolynomial& getConvolute(unsigned int k) {
	while (convolutes.size() < k) {
	  convolutes.push_back(makeConvolute(k));
	}
	return convolutes[k-1];	
      }

      const StdPolynomial makeConvolute(unsigned int k) {
	StdPolynomial p;
	if (k == 1)
	  return var(order);
	else {
	  for (int j = 1; j < order; j++)
	    p += convolute(order - j, k - 1) * var(j) * ((unsigned int)binomial[j]);
	  return p;
	}
      }

      ptr_vector<HornerPolynomial> compilePolynomials(const unsigned int order) {
	ptr_vector<HornerPolynomial> b;
	for (unsigned int k = 1; k <= order; k++)
	  b.push_back(new HornerPolynomial(convolute(k, order) / ((unsigned int) boost::math::factorial<double>(k))));
	return b;
      }

      void evalPartialDerivative(const vector<double>& f, const vector<double>& a,
				 vector<double>& target,
				 const unsigned int width, const unsigned int j) const;

      void evalValue(const vector<double>& f, const vector<double>& a,
		     vector<double>& target, const unsigned int width) const;

    public:
      /**
       * return the cacheVector without initialisation
       */
      inline static vector<Composition>& cacheVector() {
	static vector<Composition> cache({Composition(0)});
	return cache;
      }

      /**
       * Gets the cache vector and ensures that it is filled up to the given order
       */
      static vector<Composition>& cacheVectorInitialized(const unsigned int upTo);

      static Composition& ensureExistance(const unsigned int order) {
	if (cacheVector().size() <= order) {
	  cacheVectorInitialized(order);
	}
	return cacheVector()[order];
      }
        
      Composition(unsigned int o) : order(o), binomial(Multiplication::compileBinomial(o)), 
				    bell_polynomials(compilePolynomials(o)) {}
      
      void apply(const vector<double>& a, const vector<double>& b,
		 vector<double>& target, unsigned int width) const;
    };
  }
}
#endif
