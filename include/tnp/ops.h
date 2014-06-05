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

#ifndef TNP_OPS_HPP
#define TNP_OPS_HPP 1

#ifdef __cplusplus
extern "C" {
#endif

  #include <stddef.h>

  void op_prepare(int order);
  
  size_t tnp_number_payload_size(int params, int order);
  
  void op_tnp_number_to_zero(int params, int order, double* a);

  void op_tnp_number_write_variable(int params, int order, double* a, double val, int n);

  void op_tnp_number_write_constant(int params, int order, double* a, double val);

  double op_tnp_number_total_derivative(int params, double* nr, int order);

  double op_tnp_number_partial_derivative(int params, double* nr, int param);

  double op_tnp_number_mixed_derivative(int params, double* nr, int order, int param);

  void op_tnp_number_add(int params, int order, double* target, double* a, double* b);

  void op_tnp_number_dadd(int params, int order, double* target, double* a, double b);

  void op_tnp_number_mult(int params, int order, double* target, double* a, double* b);

  void op_tnp_number_dmult(int params, int order, double* target, double* a, double b);

  void op_tnp_number_sub(int params, int order, double* target, double* a, double* b);

  void op_tnp_number_dsub(int params, int order, double* target, double* a, double b);
  /*
    double* double_div(int params, int order, double* a, double* b);

    double* double_ddiv(int params, int order, double* a, double b);
  */
  void op_tnp_number_pow(int params, int order, double* target, double* a, int power);

#ifdef __cplusplus
}
#endif

#endif


