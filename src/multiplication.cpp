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

#include <tnp/ops/multiplication.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include <vector>
#include <map>

namespace tnp {
  
  using namespace std;

  vector<double> Multiplication::compileBinomial(const unsigned int order) {
    vector<double> b;
    for (int k = 0; k <= order; ++k)
      b.push_back(boost::math::binomial_coefficient<double>(order, k));
    return b;
  }

  void Multiplication::evalPartialDerivative(const vector<double>& a, const vector<double>& b,
					     vector<double>& target,
					     const unsigned int width, const unsigned int j) const {
    /*
    double d = 0;
    for (Product p : partialDerSum) {                          
      d += p.factor * a[get<0>(p.keys) * width + get<2>(p.keys) * j]
	* b[get<1>(p.keys) * width + (1 - get<2>(p.keys)) * j];
    }
    target[order*width + j] = d;
    */

    double d= 0;
    //EXPERIMENTAL, TEST (better prove)!
    for (int k = 0; k <= order; ++k) {
      const double c = binomial[k];
      d += c * a[(order - k)*width +j] * b[k * width];
      d += c * a[(order - k)*width] * b[k * width +j];
    }
    target[order*width+j] = d;
  }
  
  void Multiplication::evalValue(const vector<double>& a, const vector<double>& b,
				 vector<double>& target, const unsigned int width) const {
    /*
    double d = 0;
    for (Product p : valueSum) {
      const unsigned int x = get<0>(p.keys);
      const unsigned int y = get<1>(p.keys);      
      d += p.factor * a[x * width] * b[y * width];
    }
    target[order*width] = d;
    */

    double d = 0;
    for (int k = 0; k <= order; ++k)
      d += binomial[k] * a[(order - k)*width] * b[k * width];
    target[order*width] = d;

  }

  void Multiplication::apply(const vector<double>& a, const vector<double>& b,
			     vector<double>& target, unsigned int width) const {

    unsigned int params = width - 1;
    
    if (order > 0)
      cacheVector()[order-1].apply(a, b, target, width);
      
    evalValue(a, b, target, width);
      
    for (int j = 1; j <= params; ++j) {
      evalPartialDerivative(a, b, target, width, j);
    }
  }

  vector<Multiplication>& Multiplication::cacheVectorInitialized(const unsigned int upTo) {
    while (upTo >= cacheVector().size()) {
      cacheVector().push_back(Multiplication(cacheVector().size()));      
    }
    return cacheVector();
  }
  
  /*
  vector<Product> Multiplication::compileDerSum(const unsigned int order) {
    vector<Product> sum;

    if (order > 0) {
      const Multiplication smaller = cacheVector()[order - 1];
      
      for (Product p : smaller.partialDerSum) {
	// D(a) * I(b)
	sum.push_back(Product(p.factor, make_tuple(get<0>(p.keys) + 1, get<1>(p.keys), get<2>(p.keys))));

        // I(a) * D(b)
	sum.push_back(Product(p.factor, make_tuple(get<0>(p.keys), get<1>(p.keys) + 1, get<2>(p.keys))));
      }      
    } else {
      // (a*b)[i] == a[0] * der(b)[i] + der(a)[i] * b[0]
      sum.push_back(Product(1, make_tuple(0,0,1)));
      sum.push_back(Product(1, make_tuple(0,0,0)));
    }

    return merge(sum);
  }

  vector<Product> Multiplication::compileValueSum(const unsigned int order) {
    vector<Product> sum;
    
    if (order > 0) {
      const Multiplication smaller = cacheVector()[order - 1];

      // D|[a*b]| = D(a) * I(b) + I(a) * D(b)
      
      for (Product p : smaller.valueSum) {
	// D(a) * I(b)
	sum.push_back(Product(p.factor, make_tuple(get<0>(p.keys) + 1, get<1>(p.keys), 0)));
	
	// I(a) * D(b)
	sum.push_back(Product(p.factor, make_tuple(get<0>(p.keys), get<1>(p.keys) + 1, 0)));
      }            
    } else {
      // |[a*b]|(0,0) == a[0] * b[0]
      sum.push_back(Product(1, make_tuple(0, 0, 0) ) );      
    }

    return merge(sum);
  }
  
  vector<Product> Multiplication::merge(const vector<Product>& products) {    
    map<IntTriple, Product> m;
    
    for (Product p : products) {
      if (m.count(p.keys) != 0)
	m[p.keys] = m.find(p.keys)->second.plusOne();
      else
	m[p.keys] = p;
    }
    
    vector<Product> array;
    array.reserve(m.size());

    for (auto& element : m) {
      const Product p = element.second;
      array.push_back(p);
    }
    return array;
  }
*/

}
