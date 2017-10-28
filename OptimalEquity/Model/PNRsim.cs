/*
/ Copyright (c) 2015 Chris Rook
/
/ Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of
/ the License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by applicable law or agreed to in writing, software distributed under the License
/ is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language
/ governing permissions and limitations under the License.
/
/ Filename: PNRsim.cpp
/
/ Function: PNRsim()
/
/ Summary:
/
/ This function derives the probability of no ruin (PNR), i.e. the success probability for a given glide-path using simulation. We assume that
/ stock/bond real returns are iid normal RVs in this implementation. This is not a requirement for the model but an assumption that we make. To
/ derive an element of the gradient we need a success probability where the normal RV is replaced by an RV that follows a different distribution and
/ this is specified by the parameter prtls[4]. When deriving the elements of the Hessian matrix other success probabilities are needed and those also
/ are specified using the array prtls[4]. The array prtls[4] uses indicators to request a specific success probability (more below). This function
/ accepts parameters that set the stock/bond means, variances, covariance along with the withdrawal rate, sample size, and glide-path. At each time-
/ point we update the ruin factor RF(t). Ruin occurs if-and-only-if RF(t) becomes negative. Note that the ruin probability is not returned by this
/ function but placed into a long double array that is supplied to the function. This enables multi-threading of the ruin calculation to save processing
/ time, which is done in the wrapper function ThrdPNRsim().
/
/ Parameters:
/
/ 1.) Six member long double array of parameter settings for stock mean, stock variance, bond mean, bond variance, stock-bond covariance, expense ratio.
/ 2.) Array of long doubles holding the glide-path to use for this ruin calculation (as pointer). The corresponding array must have exactly TD elements.
/ 3.) Number of time-points for this implementation. This represents a fixed number of years to consider for the retirement horizon (i.e., 30).
/ 4.) The fixed inflation-adjusted withdrawal rate that the retiree is using during decumulation (i.e., 0.04 for a 4% withdrawal rate).
/ 5.) The sample size to use when simulating the ruin/success probability. (Note that the ruin probability is the number of times ruin occurs divided
/ by the number of trial runs. The success probability is 1 minus the ruin probability.)
/ 6.) The thread ID. PNRsim() is invoked from a wrapper function that determines the number of independent processing units on the PC running the
/ application and splits the ruin derivation into smaller jobs then runs them concurrently in separate threads and combines the results when
/ finished.
/ 7.) The array to return the calculated success probability in using the thread ID as the array element index.
/ 8.) A 4-element indicator array that specifies exactly which densities to use for the requested ruin calculation. The settings for this array are:
/ a.) prtls[] = {-1,-1,-1,-1} for all standard real/expense adjusted return densities.
/ b.) prtls[] = {i,-1,-1,-1} for gradient element "i" which uses the gradient density g() for the real/expense-adjusted return at the i-th time-
/ point.
/ c.) prtls[] = {i,j,-1,-1} for the off-diagonal Hessian elements which uses the gradient density g() at both time-points "i" and "j".
/ d.) prtls[] = {-1,-1,i,-1} for the diagonal Hessian elements which uses h1() as the real/expense-adjusted return at time-point "i".
/ e.) prtls[] = {-1,-1,-1,i} for the diagonal Hessian elements which uses h2() as the real/expense-adjusted return at time-point "i".
/
/ Return Value:
/
/ This function returns the success probability in the array supplied as probs[tnum-1].
/------------------------------------------------------------------------------------------------------------------------------------------------------------*/
using System;
using System.Diagnostics;
using Accord.Statistics.Distributions.Univariate;
using OptimalEquity.Helpers;

namespace OptimalEquity.Model
{
    public static class PNRsim
    {
        private static readonly Random Gen = new Random();
        public static double Run(double[] prms, double[] a, int fxTD, double rf0, long n, int[] partls)
        {
            // Declare local variables.
            //==========================
            double rf;
            double rtrn;
            double xlow = 0;
            double xhigh = 0;
            double yhigh = 0;
            bool ruin;
            long cntr = 0;

            // Get means/standard deviations for each timepoint and build the normal RV generator.
            var mn = new double[fxTD];
            var std = new double[fxTD];
            var rrtrn = new NormalDistribution[fxTD];

            // Check that no value of the partls[] array is equal to or greater than fxTD.
            for (int i = 0; i < 4 && partls[i] >= fxTD; ++i)
            {
                Trace.Write("ERROR: Invalid (>=TD) partial derivative specification partls[");
                Trace.Write(i);
                Trace.Write("]=");
                Trace.Write(partls[i]);
                Trace.WriteLine(".");
                Trace.WriteLine("EXITING...PNRsim()...");
                Console.Read();
                Environment.Exit(1);
            }

            // Iterate over valid time points and set parameters needed during simulation.
            for (int y = 0; y < fxTD; ++y)
            {
                // Check that only one Hessian special density is specified.
                // Set the lower/upper range and density bounds for all special densities.
                if (partls[2] != -1 && partls[3] != -1)
                {
                    Trace.WriteLine("ERROR: Only 1 Hessian special density should be used per simulation.");
                    Trace.WriteLine("EXITING...PNRsim()...");
                    Console.Read();
                    Environment.Exit(1);
                }
                else if (partls[0] != -1 || partls[1] != -1) // Gradient density g() ranges. (Confirmed for our distributional assumptions.)
                {
                    xlow = -0.15;
                    xhigh = 2.20;
                    yhigh = 6.20;
                }
                else if (partls[2] != -1) // Hessian density h1() ranges. (Confirmed for our distributional assumptions.)
                {
                    xlow = -0.10;
                    xhigh = 2.30;
                    yhigh = 5.70;
                }
                else if (partls[3] != -1) // Hessian density h2() ranges. (Confirmed for our distributional assumptions.)
                {
                    xlow = -0.15;
                    xhigh = 2.40;
                    yhigh = 5.85;
                }

                // Build arrays for means/variances based on the incoming glidepath.
                mn[y] = Funcs.m(prms, a[y]);
                std[y] = Math.Sqrt(Funcs.v(prms, a[y]));
                rrtrn[y] = new NormalDistribution(mn[y], std[y]);
            }

            // Estimate the probability.
            if (partls[0] == -1 && partls[1] == -1 && partls[2] == -1 && partls[3] == -1) // Processing for no special densities.
            {
                for (long i = 1; i <= n; ++i)
                {
                    rf = rf0;
                    ruin = false;
                    for (int y = 0; y < fxTD && ruin == false; ++y)
                    {
                        rtrn = rrtrn[y].Generate(Gen); // All returns from real/expense adjusted PDF.
                        if (rf > 0 && rtrn > rf)
                        {
                            rf = rf / (rtrn - rf);
                        }
                        else
                        {
                            ruin = true;
                        }
                    }
                    cntr += ruin ? 1 : 0;
                }
            }
            else // Processing for special densities.
            {
                for (long i = 1; i <= n; ++i)
                {
                    rf = rf0;
                    ruin = false;
                    for (int y = 0; y < fxTD && ruin == false; ++y)
                    {
                        if (partls[0] != y && partls[1] != y && partls[2] != y && partls[3] != y)
                        {
                            rtrn = rrtrn[y].Generate(Gen); // Regular density at this time point.
                        }
                        else
                        {
                            rtrn = double.NaN;
                            do
                            {
                                double unx = (xhigh - xlow) * Gen.NextDouble() + xlow;
                                double uny = yhigh * Gen.NextDouble();
                                if ((partls[0] == y || partls[1] == y) && uny <= Funcs.g(prms, a[y], unx) ||
                                    partls[2] == y && uny <= Funcs.h1(prms, a[y], unx) ||
                                    partls[3] == y && uny <= Funcs.h2(prms, a[y], unx)
                                ) // Generate RV from special density.
                                {
                                    rtrn = unx;
                                }
                            } while (double.IsNaN(rtrn));
                        }
                        if (rf > 0 && rtrn > rf)
                            rf = rf / (rtrn - rf);
                        else
                            ruin = true;
                    }
                    cntr += ruin ? 1 : 0;
                }
            }

            // Get the probability of no ruin and add it to the incoming array.
            return 1.00 - (double)cntr / n;
        }
    }
}
