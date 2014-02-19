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

#include <tnp.hpp>
#include <tnp.h>

using namespace tnp;

extern "C" {

  struct tnp_number : public NPNumber {};
  
  struct tnp_number* tnp_number_create(int params, int order) {
    return static_cast<tnp_number*>(new NPNumber(params, order));
  }
  
  void tnp_number_delete(struct tnp_number* nr) {
    delete nr;
  }
  
  int tnp_number_params(struct tnp_number* nr) {
    return nr->params();
  }

  int tnp_number_order(struct tnp_number* nr) {
    return nr->order();
  }

  double tnp_number_value(struct tnp_number* nr) {
    return nr->data()[0];
  }
  
  double tnp_number_total_derivative(struct tnp_number* nr, int order) {
    return nr->der(0, order);
  }
  
  double tnp_number_partial_derivative(struct tnp_number* nr, int param) {
    return nr->der(param, 0);
  }

  double tnp_number_mixed_derivative(struct tnp_number* nr, int order, int param) {
    return nr->der(param, order);
  }

}
