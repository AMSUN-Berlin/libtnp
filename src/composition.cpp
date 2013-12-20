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
    
    void Composition::evalValue(const vector<double>& f, const vector<double>& a,
				vector<double>& target, const unsigned int width) const {
      target[order*width] = 0.0;
      
      for (int k = 0; k < order; k++)
	target[order*width] += f[k+1] * bell_polynomials[k].eval(a, width);
    }

    inline void Composition::evalPartialDerivative(const vector<double>& f, const vector<double>& a,
						   vector<double>& target,
						   const unsigned int width, const unsigned int j) const {
      target[order * width + j] = 0.0;
      for (int k = 0; k < order; k++) 
	target[order*width + j] += f[k+2] * a[j] * bell_polynomials[k].eval(a, width) + 
	  f[k+1] * bell_polynomials[k].evalDer(a, j, width);
      
    }

    void Composition::apply(const vector<double>& f, const vector<double>& a,
			    vector<double>& target, unsigned int width) const {
      const unsigned int params = width - 1;
      if (order > 0) {
	cacheVector()[order-1].apply(f, a, target, width);
	evalValue(f, a, target, width);

	for (int j = 1; j <= params; ++j) {
	  evalPartialDerivative(f, a, target, width, j);
	}
      } else {
	target[0] = f[0];
	for (int j = 1; j <= params; ++j) {
	  target[order*width + j] = f[1] * a[order*width+j];
	}
      }
    }

    vector<Composition>& Composition::cacheVectorInitialized(const unsigned int upTo) {
      while (upTo >= cacheVector().size()) {
	cacheVector().push_back(Composition(cacheVector().size()));      
      }
      return cacheVector();
    }
  }
}
