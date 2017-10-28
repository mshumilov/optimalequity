/*
/ Copyright (c) 2015 Chris Rook
/
/ Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of
/ the License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by applicable law or agreed to in writing, software distributed under the License
/ is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language
/ governing permissions and limitations under the License.
/
/ Filename: GetGamma.cpp
/
/ Function: GetGamma()
/
/ Summary:
/
/ Each CDF call needed by this application can be expressed as a linear combination of known CDF calls to normal and gamma distributed random variables.
/ The constants and random variables needed to build the CDF call depend on the current time-point and the current special-density indicator array. Once
/ these are known the vector of gamma random variables can be built. This is combined with the appropriate normal random variable and required constants
/ in the function GetCDF().
/
/ Parameters:
/
/ 1.) A 4-element indicator array that specifies which densities should represent inflation/expense-adjusted returns at the given time-point. The
/ settings for this array are:
/ a.) prtls[] = {-1,-1,-1,-1} for all standard real/expense adjusted return densities.
/ b.) prtls[] = {i,-1,-1,-1} for gradient element "i" which uses gradient density g() for the real/expense-adjusted return at the i-th time-point.
/ c.) prtls[] = {i,j,-1,-1} for the off-diagonal Hessian elements which uses the gradient density g() at both time-points "i" and "j".
/ d.) prtls[] = {-1,-1,i,-1} for the diagonal Hessian elements which uses h1() as the real/expense-adjusted return at time-point "i".
/ e.) prtls[] = {-1,-1,-1,i} for the diagonal Hessian elements which uses h2() as the real/expense-adjusted return at time-point "i".
/ 2.) The time-point currently being processed. This value is compared with the indicator array passed in argument #1 to determine which gamma random
/ variables are needed for the CDF call.
/
/ Return Value:
/
/ This function returns a vector of gamma random variables (with appropriate shape/scale settings) to use when building the CDF required at a specific
/ time-point and special-density setting.
/------------------------------------------------------------------------------------------------------------------------------------------------------------*/
using System.Collections.Generic;
using Accord.Statistics.Distributions.Univariate;

namespace OptimalEquity.Model
{
    public static class GetGamma
    {
        public static List<GammaDistribution> Run(int[] partls, int y)
        {
            List<GammaDistribution> g = new List<GammaDistribution>();

            // Build vector of random variables.
            if (partls[3] == y)
            {
                g.Add(new GammaDistribution(2.50, 1.00));
                g.Add(new GammaDistribution(2.00, 1.00));
            }
            if (partls[0] == y || partls[1] == y || partls[2] == y || partls[3] == y)
            {
                g.Add(new GammaDistribution(1.50, 1.00));
                g.Add(new GammaDistribution(1.00, 1.00));
            }
            return g;
        }
    }
}
