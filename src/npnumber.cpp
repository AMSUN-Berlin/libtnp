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

#include <tnp/npnumber.hpp>
#include <functional>
#include <algorithm>
#include <cmath>

#include "prettyprint.hpp"

using namespace std;

namespace tnp {

  const vector<double> constant(double val, unsigned int size) {
    vector<double> v(size);
    v[0] = val;
    return v;
  }

  const vector<double> variable(double val, unsigned int n, unsigned int size) {
    vector<double> v(size);
    v[0] = val;
    v[n] = 1.0;
    return v;
  }


  NPNumber NPNumber::plus(const NPNumber& o) const {
    vector<double> c(values.size());
    transform(values.begin(),values.end(),o.values.begin(),c.begin(), std::plus<double>());

    return NPNumber(width, c);
  }

  NPNumber NPNumber::plus(const double o) const {
    vector<double> c(values.size());
    for (unsigned int i = 0; i < c.size(); ++i)
      c[i] = values[i] + o;
    return NPNumber(width, c);
  }

  NPNumber NPNumber::minus(const double o) const {
    vector<double> c(values.size());
    for (unsigned int i = 0; i < c.size(); ++i)
      c[i] = values[i] - o;
    return NPNumber(width, c);
  }

  NPNumber NPNumber::minus(const NPNumber& o) const {
    vector<double> c(values.size());   
    transform(values.begin(),values.end(),o.values.begin(),c.begin(), std::minus<double>());
    return NPNumber(width, c);
  }

  NPNumber NPNumber::times(const NPNumber& o) const {
    vector<double> c(values.size());   
    mult().apply(values, o.values, c, width);
    return NPNumber(width, c);
  }

  NPNumber NPNumber::times(const double f) const {
    vector<double> c(values.size());   
    for (unsigned int i = 0; i < c.size(); ++i)
      c[i] = values[i] * f;
    return NPNumber(width, c);
  }

  NPNumber NPNumber::pow(int n) const {
    // create the power function value and derivatives
    // [x^n, nx^(n-1), n(n-1)x^(n-2), ... ]
    vector<double> f(_order + 2);

    if (n >= 0) {
      // strictly positive power
      const int maxOrder = min(_order + 1, (unsigned int)n);
      double xk = std::pow(values[0], n - maxOrder);
      for (int i = maxOrder; i > 0; --i) {
	f[i] = xk;
	xk *= values[0];
      }
      f[0] = xk;
    } else {
      // strictly negative power
      const double inv = 1.0 / values[0];
      double xk = std::pow(inv, -n);
      for (int i = 0; i <= _order + 1; ++i) {
	f[i] = xk;
	xk *= inv;
      }
    }

    double coefficient = n;
    for (int i = 1; i <= _order + 1; ++i) {
      f[i] *= coefficient;
      coefficient *= n - i;
    }

    NPNumber newNum(params(), order());
    comp()->apply(f, values, newNum.values, width);
    
    return newNum;
  }

  NPNumber NPNumber::operator*=(const NPNumber& o) {
    vector<double> c(values);
    mult().apply(c, o.values, values, width);
    return *this;
  }

  NPNumber NPNumber::operator*=(const double o) {
    for (int i = 0; i < values.size(); i++)
      values[i] *= o;
    return *this;
  }

  NPNumber NPNumber::operator+=(const double o) {
    values[0] += o;
    return *this;
  }

  NPNumber NPNumber::operator+=(const NPNumber& o) {
    for (int i = 0; i < values.size(); i++)
      values[i] += o.values[i];
    return *this;
  }  

  std::ostream& operator<<(std::ostream& out, const NPNumber& n) {
    out << "npnumber{width=" << n.width << ", order=" << n.order() << ", values=" << n.values << "}";
    return out;
  }   
}
