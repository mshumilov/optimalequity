/*
/ Copyright (c) 2015 Chris Rook
/
/ Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of
/ the License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by applicable law or agreed to in writing, software distributed under the License
/ is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language
/ governing permissions and limitations under the License.
/
/ Filename: GetPNR.cpp
/
/ Function: GetPNR()
/
/ Summary:
/
/ This function computes a single probability of avoiding ruin by passing needed parameters to either the simulation specific or dynamic programming
/ specific functions that perform and return the actual calculation. The estimation type (simulation or dynamic programming) is specified by the user
/ in the control file.
/
/ Parameters:
/
/ 1.) Estimation type: Either dynamic programming (dp) or simulation (sim).
/ 2.) Six member long double array of parameter settings for stock mean, stock variance, bond mean, bond variance, stock-bond covariance, expense ratio.
/ 3.) Array of long doubles holding the glide-path that is the basis for the probability being requested.
/ 4.) Number of time-points for this application. This represents a fixed number of years to consider for the retirement horizon (i.e., 30).
/ 5.) The fixed inflation-adjusted withdrawal rate that the retiree is using during decumulation (i.e., 0.04 for a 4% withdrawal rate).
/ 6.) The sample size to use when simulating the ruin/success probability. (Note that the simulated ruin probability is the number of times ruin occurs
/ divided by the number of trial runs. The success probability is 1 minus the ruin probability.) (This parameter is only used when type="sim".)
/ 7.) The number of buckets to use when discretizing the ruin factor dimension. This value equals RFMax*(Discretization Precision), where these 2 values
/ are set by the user in the control file. (This parameter is only used when type="dp".)
/ 8.) The number of buckets to use when discretizing each ruin factor unit value. For example, if this value=1000 then each ruin factor unit is
/ represented by 1000 buckets. This value is specified by the user in the control file when type="dp".
/ 9.) A 4-element indicator array that specifies which densities should represent inflation/expense-adjusted returns at the given time-point. The
/ settings for this array are:
/ a.) prtls[] = {-1,-1,-1,-1} for all standard real/expense adjusted return densities.
/ b.) prtls[] = {i,-1,-1,-1} for gradient element "i" which uses gradient density g() for the real/expense-adjusted return at the i-th time-point.
/ c.) prtls[] = {i,j,-1,-1} for the off-diagonal Hessian elements which uses the gradient density g() at both time-points "i" and "j".
/ d.) prtls[] = {-1,-1,i,-1} for the diagonal Hessian elements which uses h1() as the real/expense-adjusted return at time-point "i".
/ e.) prtls[] = {-1,-1,-1,i} for the diagonal Hessian elements which uses h2() as the real/expense-adjusted return at time-point "i".
/ 10.) The number of independent processing units on the computer running the application or alternatively the number of parallel processes to use during
/ execution as specified by the user. (If user specified it is set as a parameter when invoking the application.)
/
/ Return Value:
/
/ This function returns the computed probability.
/------------------------------------------------------------------------------------------------------------------------------------------------------------*/

namespace OptimalEquity.Model
{
    public static class GetPNR
    {
        public static double Run(string type, double[] prms, double[] a, int fxTD, double rf0, long n, int nbuckets,
            int prec, int[] partls, int plproc)
        {
            double retvar = 0;
            if (type == "sim")
                retvar = ThrdPNRsim.Run(prms, a, fxTD, rf0, n, partls, plproc);
            else if (type == "dp")
                retvar = ThrdPNRdyn.Run(prms, a, fxTD, rf0, nbuckets, partls, plproc, prec);
            return retvar;
        }
    }
}
