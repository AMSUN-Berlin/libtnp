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
#ifndef TNP_OPS_MULT_HPP
#define TNP_OPS_MULT_HPP 1

#include <vector>
#include <tuple>

namespace tnp {

  using namespace std;
  
  /**
   * Efficient implementation of Multiplication based on stored triples 
   */
  class Multiplication { 

    const unsigned int order;
    const vector<double> binomial;

    void evalPartialDerivative(const double* a, const double* b,
			       double* target,
			       const unsigned int width, const unsigned int j) const;

    void evalValue(const double* a, const double* b,
		   double* target, const unsigned int width) const;

  public:
    
    /*
     * Compiles a vector containing all binomials \frac{k}{n} \forall k \in 1 \ldots n
     */
    static vector<double> compileBinomial(const unsigned int order);

    /**
     * return the cacheVector without initialisation
     */
    inline static vector<Multiplication>& cacheVector() {
      static vector<Multiplication> cache({Multiplication(0)});

      return cache;
    }

    /**
     * Gets the cache vector and ensures that it is filled up to the given order
     */
    static vector<Multiplication>& cacheVectorInitialized(const unsigned int upTo);

    static Multiplication& ensureExistance(const unsigned int order) {
      if (cacheVector().size() <= order) {
	cacheVectorInitialized(order+1);
      }
      return cacheVector()[order];
    }
        
    Multiplication(unsigned int o) : order(o), //valueSum(compileValueSum(o)), partialDerSum(compileDerSum(o)), 
				     binomial(compileBinomial(o)) {}

    void apply(const vector<double>& a, const vector<double>& b,
	       vector<double>& target, unsigned int width) const;

    void apply(const double* a, const double* b,
	       double* target, unsigned int width) const;

  };
}
#endif
