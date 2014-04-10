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

#include <tnp/ops/composition.hpp>

namespace tnp {
  namespace ops {

    using namespace std;
    /*
    void Composition::evalValue(const vector<double>& f, const vector<double>& bell,
				vector<double>& target, const unsigned int width) const {
      target[order*width] = 0.0;
      
      for (int k = 0; k < order; k++)
	target[order*width] += f[k+1] * bell[k];
    }

    inline void Composition::evalPartialDerivative(const vector<double>& f, const vector<double>& a,
						   const vector<double>& bell,
						   vector<double>& target,
						   const unsigned int width, const unsigned int j) const {
      target[order * width + j] = 0.0;
      for (int k = 0; k < order; k++) 
	target[order*width + j] += f[k+2] * a[j] * bell[k] + 
	  f[k+1] * evalDer(bell_polynomials[k], a, j, width);
      
    }
    */

    CompositionCache CompositionCache::instance;

    const StdPolynomial& Composition::convolute(unsigned int n, unsigned int k) {
      if (n < order) {
	Composition* comp = CompositionCache::staticGetInstance(n);
	return comp->getConvolute(k);
      } else {
	return getConvolute(k);
      }
    }

    const StdPolynomial& Composition::getConvolute(unsigned int k) {
      while (convolutes.size() < k) {
	convolutes.push_back(makeConvolute(k));
      }
      return convolutes[k-1];	
    }

    const StdPolynomial Composition::makeConvolute(unsigned int k) {
      StdPolynomial p;
      if (k == 1)
	return var(order);
      else {
	for (int j = 1; j < order; j++) {
	  p += convolute(order - j, k - 1) * var(j) * ((unsigned int)binomial[j]);
	}
	return p;
      }
    }

    vector<SumOfProducts> Composition::compilePolynomials(const unsigned int order) {
      vector<SumOfProducts> b;
      for (unsigned int k = 1; k <= order; k++) {
	const StdPolynomial p = convolute(order, k) / ((unsigned int) boost::math::factorial<double>(k));
	const SumOfProducts s(p);	  
	b.push_back(s);
      }
      return b;
    }

    vector<DerSumOfProducts> Composition::compileDerPolynomials(const unsigned int order) {
      vector<DerSumOfProducts> b;
      for (unsigned int k = 1; k <= order; k++) {
	const StdPolynomial p = convolute(order, k) / ((unsigned int) boost::math::factorial<double>(k));
	b.push_back(DerSumOfProducts(p));
      }
      return b;
    }
    
    void Composition::apply(const vector<double>& f, const vector<double>& a,
			    vector<double>& target, const unsigned int width) const {
      apply(f.data(), a.data(), target.data(), width);
    }
    
    void Composition::apply(const double* f, const double* a,
		 double* target, unsigned int width) const {

      const unsigned int params = width - 1;
      if (order > 0) {
	last->apply(f, a, target, width);

	for (int j = 0 ; j <= params; ++j)
	  target[order*width + j] = 0;

	for (int k = 0; k < order; k++) {
	  const SumOfProducts& bellK = bell_polynomials[k];
	  const DerSumOfProducts& dBellK = der_bell_polynomials[k];
	  const double bell = bellK.eval(a, width);
	
	  target[order*width] += f[k+1] * bell;
	  
	  for (int j = 1; j <= params; ++j) {
	    target[order*width + j] += f[k+2] * a[j] * bell + f[k+1] * dBellK.eval(a, width, j);
	  }
	}
      } else {
	target[0] = f[0];
	for (int j = 1; j <= params; ++j) {
	  target[order*width + j] = f[1] * a[order*width+j];
	}
      }
    }
  }
}
