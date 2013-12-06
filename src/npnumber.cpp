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

using namespace std;

namespace tnp {

  const std::vector<double> constant(double val, unsigned int size) {
    std::vector<double> v(size);
    v[0] = val;
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


}
