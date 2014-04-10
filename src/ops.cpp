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

#include <tnp/ops.h>
#include <tnp/ops/multiplication.hpp>
#include <tnp/ops/composition.hpp>

#include <algorithm>
#include <cmath>

extern "C" {
  
  void op_prepare(int order) {
    tnp::Multiplication::ensureExistance(order);
  }
  
  size_t tnp_number_payload_size(int params, int order) {
    return sizeof(double) * (params+1) * (order+1);
  }
  
  void op_tnp_number_to_zero(int params, int order, double* a) {
    memset(a, 0, tnp_number_payload_size(params, order));
  }

  void op_tnp_number_write_variable(int params, int order, double* a, double val, int n) {
    op_tnp_number_to_zero(params, order, a);
    a[0] = val;
    a[n] = 1.0;
  }

  void op_tnp_number_write_constant(int params, int order, double* a, double val) {
    op_tnp_number_to_zero(params, order, a);
    a[0] = val;
  }

  double op_tnp_number_total_derivative(int params, double* nr, int order) {
    return op_tnp_number_mixed_derivative(params, nr, order, 0);
  }

  double op_tnp_number_partial_derivative(int params, double* nr, int param) {
    return op_tnp_number_mixed_derivative(params, nr, 0, param);
  }

  double op_tnp_number_mixed_derivative(int params, double* nr, int order, int param) {
    const int width = params + 1;
    return nr[order*width + param];        
  }

  void op_tnp_number_add(int params, int order, double* target, double* a, double* b) {
    const size_t size = (params + 1) * (order + 1);
    for (size_t i = 0; i < size; i++)
      target[i] = a[i] + b[i];
  }

  void op_tnp_number_dadd(int params, int order, double* target, double* a, double b) {
    const size_t size = (params + 1) * (order + 1);
    std::copy(a, a + size, target);
    a[0] += b;
  }

  void op_tnp_number_mult(int params, int order, double* target, double* a, double* b) {
    tnp::Multiplication::cacheVector()[order].apply(a, b, target, params + 1);
  }

  void op_tnp_number_dmult(int params, int order, double* target, double* a, double b) {
    const size_t size = (params + 1) * (order + 1);
    for (size_t i = 0; i < size; i++)
      target[i] = a[i] + b;
  }

  void op_tnp_number_sub(int params, int order, double* target, double* a, double* b) {
    const size_t size = (params + 1) * (order + 1);
    for (size_t i = 0; i < size; i++)
      target[i] = a[i] - b[i];
  }

  void op_tnp_number_dsub(int params, int order, double* target, double* a, double b) {
    const size_t size = (params + 1) * (order + 1);
    std::copy(a, a + size, target);
    a[0] -= b;
  }

  /*
    double* double_div(int params, int order, double* a, double* b);

    double* double_ddiv(int params, int order, double* a, double b);
  */
  void op_tnp_number_pow(int params, int order, double* target, double* a, int n) {
    // create the power function value and derivatives
    // [x^n, nx^(n-1), n(n-1)x^(n-2), ... ]
    double f[order + 2];

    if (n >= 0) {
      // strictly positive power
      const int maxOrder = std::min(order + 1, n);
      double xk = std::pow(a[0], n - maxOrder);
      for (int i = maxOrder; i > 0; --i) {
	f[i] = xk;
	xk *= a[0];
      }
      f[0] = xk;
    } else {
      // strictly negative power
      const double inv = 1.0 / a[0];
      double xk = std::pow(inv, -n);
      for (int i = 0; i <= order + 1; ++i) {
	f[i] = xk;
	xk *= inv;
      }
    }

    double coefficient = n;
    for (int i = 1; i <= order + 1; ++i) {
      f[i] *= coefficient;
      coefficient *= n - i;
    }

    tnp::ops::CompositionCache::staticGetInstance(order)->apply(f, a, target, params+1);
  }

}

