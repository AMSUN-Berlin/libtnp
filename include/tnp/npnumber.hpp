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
#ifndef TNP_NPNUMBER_HPP
#define TNP_NPNUMBER_HPP 1

#include <iostream>
#include <vector>
#include <iterator>

namespace tnp {


  const std::vector<double> constant(double val, unsigned int size) {
    std::vector<double> v(size);
    v[0] = val;
    return v;
  }

  /**
   * An AD-number of arbitrary width and depth
   */
  class NPNumber {
    const unsigned int width;
    const std::vector<double> values;
    
  public:
    NPNumber(unsigned int params, unsigned int order) : width(params+1), values(width * order+1, 0.0) {}

    NPNumber(unsigned int width, const std::vector<double>& values) : width(width), values(values) {}

    NPNumber(unsigned int params, unsigned int order, double value) : width(params+1), values(constant(value, width * (order+1))) {}

    ~NPNumber() {}

    NPNumber plus(const NPNumber& o) const;

    NPNumber minus(const NPNumber& o) const;
    
    unsigned int order() const { return values.size() / width - 1; }
    
    unsigned int params() const { return width - 1; }

    const std::vector<double>& data() const { return values; }

    bool operator==(const NPNumber& o) const {
      return values == o.values;
    }

    NPNumber operator+(const NPNumber& o) const {
      return plus(o);
    }

    NPNumber operator-(const NPNumber& o) const {
      return minus(o);
    }

    NPNumber operator*(const NPNumber& o) const {
      //TODO
      return *this;
    }

    NPNumber operator/(const NPNumber& o) const {
      //TODO
      return *this;
    }

    friend std::ostream& operator<<(std::ostream& out, const NPNumber& n) {
      out << "npnumber{width=" << n.width << ", order=" << n.order() << ", values=[";
      std::copy(n.values.begin(), n.values.end(), std::ostream_iterator<double>(out, ", "));
      out << "]}";
      return out;
    }   
  };

}

#endif
