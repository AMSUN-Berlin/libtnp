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
 * Test - instances of NPNumber
 */
namespace tnp {
  namespace test {

    std::vector<std::pair<unsigned int, unsigned int>> makeSizes() {
      std::vector<std::pair<unsigned int, unsigned int>> sizes({
	  std::make_pair(0,0),
	    std::make_pair(0,2),
	    std::make_pair(2,2),
	    std::make_pair(10,5)});
      return sizes;
    }

    const std::vector<std::pair<unsigned int, unsigned int>> testDimensions(makeSizes());

    std::vector<NPNumber> makeNumbers() {
      std::vector<NPNumber> v;
      
      const std::vector<unsigned int> testParams({1,2,3,4,5,6,7,8,9,10});
      
      const std::vector<unsigned int> testOrders({1,2,3,4,5});

      const std::vector<double> constantsField({
	  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 });

      const std::vector<double> testField({
	  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 });

      for (unsigned int p : testParams) {
	for (int o : testOrders) {	  
	  /* always test on constants */
	  for (double d : constantsField) {
	    v.push_back(NPNumber(p, o, d));
	  }
	  
	  const unsigned int s = (p+1)*(o+1);
	  /* Try to make some numbers */
	  for (int n = 0; (n+s) <= testField.size(); n+=s) {
	    const std::vector<double> vals(testField.begin() + n, testField.begin() + n + s);
	    v.push_back(NPNumber(p+1, vals));
	  }
	}
      }
      
      std::cout << "returning " << v.size() << " test numbers of size: " << sizeof(NPNumber) << std::endl;
      return v;
    }

    const std::vector<NPNumber>& testNumbers() {
      static const std::vector<NPNumber> nums(makeNumbers());
      return nums;
    }
  }
}

#endif
