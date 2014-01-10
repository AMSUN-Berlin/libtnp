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
#include <boost/timer/timer.hpp>

#include <tnp/polynomial.hpp>

namespace tnp {
  namespace test {

    NPNumber bellLoopTimed(const SumOfProducts& p, NPNumber result, unsigned int n) {
      boost::timer::auto_cpu_timer t;
      double res = 0.0;

      SumOfProducts::evals = 0;
      for (unsigned int i = 0; i < n; ++i) 
	res += p.eval(result.data().data(), result.params() + 1);      
      cout << SumOfProducts::evals << " evaluations" << endl;
      return result + res;
    }

    NPNumber bellLoop(NPNumber result, unsigned int n) {
      cout << "Running bell performance evaluation " << result.order() << endl;
	for (auto p : tnp::CompositionCache::staticGetInstance(result.order())->bell()) {
	  cout << p.sum.size() << " products " << endl;
	  bellLoopTimed(p, result, n);
	}

      return result;
    }

    NPNumber assignLoop(NPNumber result, unsigned int n) {
      boost::timer::auto_cpu_timer t;
      for(unsigned int i = 0; i < n; ++i)
	result = result;
      return result;
    }

    NPNumber multiplyLoop(NPNumber result, unsigned int n) {
      boost::timer::auto_cpu_timer t;

      for (unsigned int i = 0; i < n; ++i)
	result = result * result;
      return result;
    }

    NPNumber squareLoop(NPNumber result, unsigned int n) {
      cout << "Running pow() performance evaluation, order=" << result.order() << ", params=" << result.params() << endl;
      cout << "Composition order: " << result.comp()->order << endl;
      boost::timer::auto_cpu_timer t;
      SumOfProducts::evals = 0;
      for (unsigned int i = 0; i < n; ++i)
	result = result.pow(2);
      return result;
      cout << SumOfProducts::evals << " evaluations" << endl; 
    }

  }
}
