/*
/ Copyright (c) 2015 Chris Rook
/
/ Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of
/ the License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by applicable law or agreed to in writing, software distributed under the License
/ is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language
/ governing permissions and limitations under the License.
/
/ Filename: GetCDF.cpp
/
/ Function: GetCDF()
/
/ Summary: Single function that returns any CDF value required by this application. All CDF calls derived in this application can be expressed as a linear
/ combination of CDF calls from normal and gamma random variables. The constants and random variables are passed as parameters and this function
/ combines them and computes then returns the CDF value for the value x (set in last argument). CDF calls are only needed when estimation is via
/ dynamic programming. When simulation is the estimation method the PDFs are used. Having one function return any CDF value that is required
/ allows for only one function to compute needed probabilities via dynamic programming.
/
/ Parameters:
/
/ 1.) Array of precomputed moments & constants for the given time-point and asset allocation (to save processing time): m, mp, v, vp, vpp, mva, kh1.
/ 2.) Vector of constants. All CDFs in this application can be expressed as a linear combination of known CDF calls. This is the constant vector for
/ that linear combination.
/ 3.) A normal random variable (class type=boost distribution) reflecting the inflation/expense-adjusted return for this time-point and asset allocation.
/ 4.) A vector of gamma distributed random variables (class type=boost distribution) needed to build the CDF for a given time-point and asset allocation.
/ The CDF G() of g(), and the CDF H1() of h1() both use 2 gamma RVs. The CDF H2() of h2() uses 4 gamma RVs.
/ 5.) A value on the axis of the random variable of interest. The value returned is the probability of being to the left of this value.
/
/ Return Value:
/
/ This function returns the CDF of x, where the CDF is a linear combination (specific form) of the constants and random variables passed as arguments.
/ This function can return CDF values for F() of f(), G() of g(), H1() of h1(), and H2() of h2().
/------------------------------------------------------------------------------------------------------------------------------------------------------------*/
using System;
using System.Collections.Generic;
using Accord.Statistics.Distributions.Univariate;

namespace OptimalEquity.Model
{
    public static class GetCDF
    {
        public static double Run(double[] mts, List<double> cLocal, NormalDistribution n, List<GammaDistribution> g,
            double x)
        {
            double retval = 0;
            double sgn = x <= mts[0] ? 1 : -1;
            int hival = 1;

            // Array to hold Gamma CDF values.
            var f = new double[g.Count + 1];
            f[0] = 0.00;

            // Derive & return the CDF value for the given parameters.
            for (int i = 1; i <= g.Count; ++i)
            {
                f[i] = 1.00 - Math.Pow(sgn, i) *
                       g[i - 1].DistributionFunction(Math.Pow((x - mts[0]) / Math.Sqrt(2.00 * mts[2]), 2));
                hival = Math.Abs(f[i] - (1.0 + Math.Pow(-1.0, i + 1))) < 1e-15 ? hival : 0;
                retval = retval + cLocal[i] * f[i];
            }

            // Return the CDF value, making sure the highest possible value is 1.00.
            if (hival == 1 && Math.Abs(n.DistributionFunction(x) - 1.00) < 1e-15)
            {
                return 1.00;
            }
            return cLocal[0] * (retval + cLocal[cLocal.Count - 1] * n.DistributionFunction(x));
        }
    }
}
