/*
/ Copyright (c) 2015 Chris Rook
/
/ Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of
/ the License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by applicable law or agreed to in writing, software distributed under the License
/ is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language
/ governing permissions and limitations under the License.
/
/ Filename: ThrdPNRsim.cpp
/
/ Function: ThrdPNRsim()
/
/ Summary:
/
/ This function serves as a wrapper for the PNRsim() function and it splits the job of simulating a ruin probability across the various independent
/ processing units on the computer running the application. This function is passed the number of independent processing units on the computer and it
/ divides the simulation sample size by this number and launches simultaneous calls to PNRsim() then averages the resulting probability across calls for
/ a final value that is returned to the calling program. When invoking this function expect a single success probability to be returned.
/
/ Parameters:
/
/ 1.) Six member long double array of parameter settings for stock mean, stock variance, bond mean, bond variance, stock-bond covariance, expense ratio.
/ 2.) Array of long doubles holding the glide-path to use for this ruin calculation (as pointer). The corresponding array must have exactly TD elements.
/ (Note that we only consider alphas between the low-volatility portfolio and 1.00. If any alpha provided in this array is not between the low-
/ volatility portfolio and 1.00 it is forced to the nearest acceptable value.)
/ 3.) Number of time-points for this implementation. This represents a fixed number of years to consider for the retirement horizon (i.e., 30).
/ 4.) The fixed inflation-adjusted withdrawal rate that the retiree is using during decumulation (i.e., 0.04 for a 4% withdrawal rate).
/ 5.) The sample size to use when simulating the ruin/success probability. (Note that the ruin probability is the number of times ruin occurs divided
/ by the number of trial runs. The success probability is 1 minus the ruin probability.)
/ 6.) A 4-element indicator array that specifies exactly which densities to use for the requested ruin calculation. The settings for this array are:
/ a.) prtls[] = {-1,-1,-1,-1} for all standard real/expense adjusted return densities.
/ b.) prtls[] = {i,-1,-1,-1} for gradient element "i" which uses the gradient density g() for the real/expense-adjusted rtrn at the ith time point.
/ c.) prtls[] = {i,j,-1,-1} for the off-diagonal Hessian elements which uses the gradient density g() at both time-points "i" and "j".
/ d.) prtls[] = {-1,-1,i,-1} for the diagonal Hessian elements which uses h1() as the real/expense-adjusted return at time-point "i".
/ e.) prtls[] = {-1,-1,-1,i} for the diagonal Hessian elements which uses h2() as the real/expense-adjusted return at time-point "i".
/ 7.) The number of independent processing units on the computer running the application or alternatively the number of parallel processes to use during
/ execution as specified by the user. (If user specified it is set as a parameter when invoking the application.)
/
/ Return Value:
/
/ This function returns the simulated success probability.
/------------------------------------------------------------------------------------------------------------------------------------------------------------*/
using System;
using System.Diagnostics;
using System.Threading.Tasks;
using OptimalEquity.Helpers;

namespace OptimalEquity.Model
{
    public static class ThrdPNRsim
    {
        public static double Run(double[] prms, double[] a, int fxTD, double rf0, long n, int[] partls, int plproc)
        {
            double probnr = 0;
            // Check that all equity allocations are between MVA and 1.00.
            for (int y = 0; y < fxTD; ++y)
            {
                if (a[y] < Funcs.mva(prms) + 0.0001 || a[y] > 1.00)
                {
                    if (a[y] > 1.00)
                    {
                        a[y] = 1.00;
                    }
                    else if (a[y] < Funcs.mva(prms) + 0.0001)
                    {
                        a[y] = Funcs.mva(prms) + 0.0001;
                    }
                }
            }
            var t = new Task<double>[plproc];
            // Separate processing when using simulation vs DP to determine PNR.
            for (int i = 0; i < plproc; ++i)
            {
                var idx = i + 1;
                t[i] = Task.Run(() =>
                {
                    try
                    {
                        return PNRsim.Run(prms, a, fxTD, rf0, n / plproc, partls);
                    }
                    catch (Exception ex)
                    {
                        Trace.WriteLine(
                            $"Calculation of #{idx} task failed with error: {ex.Message}. Stacktrace: {ex.StackTrace}");
                        return double.NaN;
                    }
                });
            }

            for (int i = 0; i < plproc; ++i)
                if (!double.IsNaN(t[i].Result))
                    // average probabilities across threads
                    probnr += t[i].Result / plproc;

            // the single simulated PNR is returned.
            return probnr;
        }
    }
}
