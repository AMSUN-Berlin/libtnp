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

#include <tnp/ops/multiplication.hpp>

namespace tnp {

  const std::vector<double> constant(double val, unsigned int size);

  /**
   * An AD-number of arbitrary width and depth
   */
  class NPNumber {

    inline const Multiplication& mult() const {
      return Multiplication::cacheVector()[_order];
    }

    unsigned int _order;
    unsigned int width;
    std::vector<double> values;    
  public:

    NPNumber(unsigned int params, unsigned int order) : _order(order), width(params+1), 
							values(width * order+1, 0.0) {
      Multiplication::ensureExistance(_order);
    }
							
    NPNumber(unsigned int width, const std::vector<double>& values) : _order(values.size() / width - 1), 
								      width(width), 
								      values(values) {
      Multiplication::ensureExistance(_order);
    }

    NPNumber(unsigned int params, unsigned int order, double value) : _order(order), width(params+1),
								      values(constant(value, width * (order+1))) {
      Multiplication::ensureExistance(_order);
    }
    
    ~NPNumber() {}

    /* addition */
    NPNumber plus(const NPNumber& o) const;
    NPNumber plus(const unsigned int s) const;
    NPNumber plus(const double s) const;

    NPNumber minus(const NPNumber& o) const;
    NPNumber minus(const unsigned int s) const;
    NPNumber minus(const double s) const;

    /* multiplication */
    NPNumber times(const NPNumber& o) const;
    NPNumber times(const unsigned int f) const;
    NPNumber times(const double f) const;
    
    unsigned int order() const { return _order; }
    unsigned int params() const { return width - 1; }

    const std::vector<double>& data() const { return values; }

    inline double der(const unsigned int param, const unsigned int order) const 
    { return values[width*order + param]; }

    /* operators */
    bool operator==(const NPNumber& o) const {
      return values == o.values;
    }

    NPNumber operator+(const NPNumber& o) const {
      return plus(o);
    }

    NPNumber operator+(const double o) const {
      return plus(o);
    }

    NPNumber operator-(const NPNumber& o) const {
      return minus(o);
    }

    NPNumber operator-(const double o) const {
      return minus(o);
    }

    NPNumber operator*(const NPNumber& o) const {
      return times(o);
    }

    NPNumber operator*(const double f) const {
      return times(f);
    }

    NPNumber operator/(const NPNumber& o) const {
      //TODO
      return *this;
    }
    
    NPNumber asParameter(unsigned int param) const {
      vector<double> args(values);
      /* dx/dx = 1.0, ddx/ddx = 1.0, ... */
      while (param < args.size()) {
	args[param] = 1.0;
	param += width;
      }
      return NPNumber(width, args);
    }

    /**
     * returns the np-number that represents the free variable t of a n-ary composed 
     * function f(x_i(t))
     */
    static NPNumber freeVar(unsigned int params, unsigned int order, double value) {
      NPNumber num(params, order, value);
      if (order > 0)
	num.values[num.width] = 1.0; //dt/dt = 1
      return num;
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
