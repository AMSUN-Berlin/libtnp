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
#ifndef TNP_TEST_ANALYTIC_FUNCTIONS_HPP
#define TNP_TEST_ANALYTIC_FUNCTIONS_HPP 1

#include <boost/test/floating_point_comparison.hpp>

#include <tnp/npnumber.hpp>

#include <iostream>
#include <math.h> 

#include "prettyprint.hpp"
#include "numberGenerator.hpp"

namespace tnp {
  namespace test {
    using namespace std;

    class UnaryAnalyticFunction {
      UnaryAnalyticFunction* der = NULL;

      virtual UnaryAnalyticFunction* _derivative() const = 0;

    public:

      virtual ~UnaryAnalyticFunction() {
	if (der != NULL)
	  delete der;
      }

      UnaryAnalyticFunction* derivative() {
	if (der == NULL)
	  der = _derivative();

	return der;
      }
            
      virtual NPNumber eval(const unsigned int order, const double argValue) const = 0;
      
      virtual double eval(const double arg) const = 0;      

      virtual std::ostream& to(std::ostream&) const = 0;
    };
    
    std::ostream& operator<<(std::ostream& os, const UnaryAnalyticFunction& obj) {
      return obj.to(os);
    }

    class Constant : public UnaryAnalyticFunction {
      double value;
      
      UnaryAnalyticFunction* _derivative() const {
	return new Constant(0.0);
      }

    public:
      Constant(double x) : value(x) {}
      
      NPNumber eval(const unsigned int order, const double argValue) const {
	return NPNumber(0, order, value);
      }
      
      double eval(const double arg) const {
	return value;
      }

      std::ostream& to(std::ostream& o) const { 
	return o << "{f(x) = " << value << "}";
      }

    };

    class Linear : public UnaryAnalyticFunction {
      double factor;

      UnaryAnalyticFunction* _derivative() const {
	return new Constant(factor);
      }

    public:
      Linear() : factor(1.0) {}
      
      Linear(double d) : factor(d) {}

      NPNumber eval(const unsigned int order, const double arg) const {
	return NPNumber::freeVar(0, order, arg) * factor;
      }
      
      double eval(const double arg) const {
	return arg * factor;
      }
           
      std::ostream& to(std::ostream& o) const { 
	return o << "{f(x) = " << factor << " * x}";
      }
    };

    class Polynom : public UnaryAnalyticFunction {
      int power;
      double factor;

      UnaryAnalyticFunction* _derivative() const {
	if (power == 2)
	  return new Linear(power * factor);
	else
	  return new Polynom(power - 1, power*factor);
      }
    public:
      Polynom() : power(1), factor(1.0) {}
      Polynom(int p, double f) : power(p), factor(f) {}
      

      NPNumber eval(const unsigned int order, const double arg) const {
	const NPNumber b = NPNumber::freeVar(0, order, arg);
	NPNumber n(0, order, 1.0);
	
	if (power > 0) {
	  //TODO: use ^ operator
	  for (int i = 0; i < power; ++i)
	    n = n * b;
	} else {
	  //TODO: use ^ operator
	  for (int i = 0; i < power; ++i)
	    n = n / b;
	}
	return n * factor;
      }
      
      double eval(const double arg) const {
	return pow(arg, power) * factor;
      }
     
      std::ostream& to(std::ostream& o) const { 
	return o << "{f(x) = " << factor << " * x^" << power << "}";
      }
    };
    
    const vector<UnaryAnalyticFunction*> testFunctions() {
      static Polynom p33(3,3);
      static Polynom p21(3,3);
      static Linear l4(4);
      
      static vector<UnaryAnalyticFunction*> fs({&p33, &p21, &l4});
      return fs;
    }

    class UnaryAnalyticFunctionTest {
      static vector<double> values(const vector<NPNumber>& args) {
	vector<double> dArgs(args.size());
	for (int i = 0; i < args.size(); ++i)
	  dArgs[i] = args[i].data()[0];
	return dArgs;
      }


      static vector<UnaryAnalyticFunctionTest> makeTestCases() {
	vector<UnaryAnalyticFunctionTest> tests;
	for (UnaryAnalyticFunction* pf : testFunctions()) {	  
	  for (int i = 0; i < testNumbers().size(); ++i) {
	    tests.push_back(UnaryAnalyticFunctionTest(pf, testNumbers()[i].order(), testNumbers()[i].data()[0]));
	  }	    	  
	}
	return tests;
      }

    public:
      static const vector<UnaryAnalyticFunctionTest>& testCases() {
	static vector<UnaryAnalyticFunctionTest> tests(makeTestCases());
	return tests;
      }

      const int order;
      const double arg;
      UnaryAnalyticFunction* fun;

      UnaryAnalyticFunctionTest(UnaryAnalyticFunction* f, int o, double a) : order(o), arg(a), fun(f) {}
    };

    void testUnaryAnalyticFunction(const UnaryAnalyticFunctionTest& test) {
      const NPNumber npRes = test.fun->eval(test.order, test.arg);
      
      UnaryAnalyticFunction* f = test.fun;

      for (int o = 0; o <= test.order; o++) {
	const double res = f->eval(test.arg);
	BOOST_CHECK_MESSAGE(boost::test_tools::check_is_close(res, npRes.der(0, o), 
								boost::test_tools::percent_tolerance(1e-10)), 
			    "Testing function: " << (*test.fun) << "\n" <<
			    "Derivation: " << o << " = ideal function: " << (*f) << "\n" <<
			    "argument: " << test.arg << " expected: " << res << "\n" <<
			    "result: " << npRes << "\n") ;
	f = f->derivative();
      }
    }
  }
}

#endif
