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
#ifndef TNP_TEST_NUMBERS_HPP
#define TNP_TEST_NUMBERS_HPP 1

#include <vector>
#include <tnp/npnumber.hpp>

/**
 * These test cases take single NPNumbers and assert a 
 * algebraic identity upon a unary operation (e.g. adding a constant) 
 * on them.
 */
namespace tnp {
  namespace test {

    const std::vector<unsigned int> testParams({1,2,3,4,5});

    const std::vector<unsigned int> testOrders({1,2,3,4,5,6,7,8,9,10});

    const std::vector<double> constantsField({
      0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    });

    const std::vector<double> testField({
      0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
    });

    std::vector<NPNumber> makeNumbers() {
      std::vector<NPNumber> v;
      
      for (unsigned int p : testParams) {
	for (int o : testOrders) {
	  
	  /* always test on constants */
	  for (double d : constantsField) {
	    v.push_back(NPNumber(p, o, d));
	  }

	  /* Try to make some numbers */
	  for (int n = 0; n < testField.size(); n+=(p+1)*(o+1)) {
	    const std::vector<double> vals(testField.begin() + n, testField.begin() + n + (p+1)*(o+1));
	    v.push_back(NPNumber(p+1, vals));
	  }
	}
      }
     
      return v;
    }

    const std::vector<NPNumber> testNumbers = makeNumbers();
  }
}

#endif
