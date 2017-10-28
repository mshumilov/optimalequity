/*
/ Copyright (c) 2015 Chris Rook
/
/ Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of
/ the License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by applicable law or agreed to in writing, software distributed under the License
/ is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language
/ governing permissions and limitations under the License.
/
/ Filename: GetConst.cpp
/
/ Function: GetConst()
/
/ Summary:
/
/ Each CDF call used in this application can be expressed as a linear combination of CDF calls to known random variables. The random variables are
/ normal and gamma distributed. The value of the constants in the linear combination depend on whether or not special densities are used to represent
/ inflation/expense-adjusted returns at each time-point. This function takes the information for the given time-point, including the moments, the time-
/ point, and the special density indicator array and builds the appropriately sized vector of constants. A similar function will build the appropriately
/ sized vector of random variables, and the GetCDF() function combines the constants and random variables to return the needed CDF value.
/
/ Parameters:
/
/ 1.) Array of precomputed moments & constants for the given time-point and asset allocation (to save processing time): m, mp, v, vp, vpp, mva, kh1.
/ 2.) A 4-element indicator array that specifies which densities should represent inflation/expense-adjusted returns at the given time-point.
/ The settings for this array are:
/ a.) prtls[] = {-1,-1,-1,-1} for all standard real/expense adjusted return densities.
/ b.) prtls[] = {i,-1,-1,-1} for gradient element "i" which uses gradient density g() for the real/expense-adjusted return at the i-th
/ time-point.
/ c.) prtls[] = {i,j,-1,-1} for the off-diagonal Hessian elements which uses the gradient density g() at both time-points "i" and "j".
/ d.) prtls[] = {-1,-1,i,-1} for the diagonal Hessian elements which uses h1() as the real/expense-adjusted return at time-point "i".
/ e.) prtls[] = {-1,-1,-1,i} for the diagonal Hessian elements which uses h2() as the real/expense-adjusted return at time-point "i".
/ 3.) The time-point currently being processed. This value is compared with the indicator array passed in argument #2 to determine which
/ constants are needed for the CDF call.
/
/ Return Value:
/
/ This function returns the vector of constants needed to build the CDF for a given asset allocation and time-point. The size of the vector and
/ constant values returned depend on the density being used to represent inflation/expense-adjusted returns for the given time-point.
/-----------------------------------------------------------------------------------------------------------------------------------------------------------*/
using System;
using System.Collections.Generic;

namespace OptimalEquity.Model
{
    public static class GetConst
    {
        public static List<double> Run(double[] mts, int[] partls, int y)
        {
            List<double> c = new List<double>();

            // Build vector of constants.
            if (partls[0] != y && partls[1] != y && partls[2] != y && partls[3] != y)
            {
                c.Add(1.00);
            }
            else if (partls[0] == y || partls[1] == y)
            {
                c.Add(1.00 / (Math.Pow(mts[1], 2) * Math.Pow(mts[2], 2) + mts[2] * Math.Pow(mts[3], 2)));
                c.Add(Math.Pow(mts[3], 2) * mts[2] / 2.00);
                c.Add(-mts[3] * mts[1] * mts[2] * Math.Sqrt(2.00 * mts[2] / Math.PI));
                c.Add(Math.Pow(mts[1], 2) * Math.Pow(mts[2], 2));
            }
            else if (partls[2] == y)
            {
                c.Add(1.00 / (mts[2] + Math.Pow(mts[6], 2)));
                c.Add(mts[2] / 2.00);
                c.Add(-Math.Sqrt((2.00 * mts[2]) / Math.PI) * mts[6]);
                c.Add(Math.Pow(mts[6], 2));
            }
            else if (partls[3] == y)
            {
                c.Add(2.00 / (Math.Pow(mts[3], 2) + 2.00 * mts[2] * Math.Pow(mts[1], 2)));
                c.Add(3.00 * Math.Pow(mts[3], 2) / 8.00);
                c.Add(-Math.Sqrt(2.00 * mts[2] / Math.PI) * mts[3] * mts[1]);
                c.Add((2.00 * Math.Pow(mts[1], 2) * mts[2] - Math.Pow(mts[3], 2)) / 4.00);
                c.Add(Math.Sqrt(mts[2] / (2.00 * Math.PI)) * mts[3] * mts[1]);
                c.Add(Math.Pow(mts[3], 2) / 4.00);
            }
            return c;
        }
    }
}
