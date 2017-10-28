/*
/ Copyright (c) 2015 Chris Rook
/
/ Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of
/ the License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by applicable law or agreed to in writing, software distributed under the License
/ is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language
/ governing permissions and limitations under the License.
/
/ Filename: BldGrad.cpp
/
/ Function: BldGrad()
/
/ Summary:
/
/ This function computes each gradient element and when using simulation tests the null hypothesis (Ho) that it equals zero using the type I error alpha
/ set in the control file (2nd alpha). If simulation is used and this hypothesis test is not rejected then the element is set to zero in a 2nd vector
/ termed the adjusted gradient. When simulation is used, each gradient element is subject to sampling error and if the noise is retained as the
/ direction of steepest ascent then we will move in this direction when climbing which will need to be reversed during future iterations. By setting the
/ gradient element to zero when we are at the optimal allocation for a given time point during a given iteration we do not modify that equity ratio when
/ climbing. The incoming appropriately sized array is populated with the gradient (when estimating via dp) or adjusted gradient (when estimating via
/ simulation) entries and the maximum absolute effective value of this array is the return value for the function. The effective gradient entry is the
/ value for each time point that will not drive the equity ratio outside of the feasible region. This maximum effective gradient value can be used to
/ determine convergence. Note that when we are operating along the border region, the gradient will continue to point in the direction of steepest
/ ascent even if we cannot climb any further in that direction. Therefore, without using the maximum effective gradient value, the procedure would fail
/ to converge. The goal is to drive each element of the gradient to a value of zero since this implies that there is no direction we can move in to
/ improve the success probability, thus we are at a local optimum. Both the adjusted (when using simulation) and unadjusted gradient vectors are printed
/ to the screen in block form. (Note that the goal is to drive the unadjusted gradient vector to zero in all dimensions and any technique that achieves
/ this goal faster is valid as long as it works, which is the reason for testing each element against zero and eliminating values that are not
/ significantly different from zero.) When this function is invoked, an estimation method is provided and determines how each gradient element is
/ estimated. The two estimation choices are dynamic program (dp) or simulation (sim). When simulation is used each quantity estimated is subject to
/ sampling error, and when dynamic programming is used each quantity is subject to approximation error.
/
/ Parameters:
/
/ 1.) Estimation type: Either dynamic programming (dp) or simulation (sim).
/ 2.) Six member long double array of parameter settings for stock mean, stock variance, bond mean, bond variance, stock-bond covariance, expense ratio.
/ 3.) Array of long doubles holding the glide-path that is the basis for computing the gradient (as pointer). The corresponding array must have exactly
/ TD elements.
/ 4.) Number of time-points for this implementation. This represents a fixed number of years to consider for the retirement horizon (i.e., 30).
/ 5.) The fixed inflation-adjusted withdrawal rate that the retiree is using during decumulation (i.e., 0.04 for a 4% withdrawal rate).
/ 6.) The sample size to use when simulating the ruin/success probability. (Note that the ruin probability is the number of times ruin occurs divided
/ by the number of trial runs. The success probability is 1 minus the ruin probability.) Here the probability will use the special density g() at
/ the given time point since we are building gradient elements. (This parameter is only used when type="sim".)
/ 7.) The number of buckets to use when discretizing the ruin factor dimension. This value equals RFMax*(Discretization Precision), where these 2 values
/ are set by the user in the control file. (This parameter is only used when type="dp".)
/ 8.) The number of buckets to use when discretizing each ruin factor unit value. For example, if this value=1000 then each ruin factor unit is
/ represented by 1000 buckets. This value is specified by the user in the control file when type="dp".
/ 9.) The success probability (i.e., probability of no ruin) for the incoming glide-path. (Use a large sample size when estimating this if type="sim"
/ since it is subject to sampling error and used when deriving each gradient value. A poor representation will negatively impact each gradient
/ element.)
/ 10.) The sample size used to simulate the success probability passed in with the above parameter. (We need to know the sample size when testing the
/ hypothesis that each gradient element equals zero. (This parameter is only used when type="sim".)
/ 11.) The alpha value used to test the hypothesis that the gradient element equals zero. If this hypothesis holds then the value can be considered
/ sampling noise and we set it to zero when climbing. (The reason for doing this is to achieve convergence faster. This value is set by the user in
/ the control file, and it is only used when type="sim".)
/ 12.) The number of independent processing units on the computer running the application or alternatively the number of parallel processes to use during
/ execution as specified by the user. (If user specified it is set as a parameter when invoking the application.)
/ 13.) A pointer to an empty array of appropriate size (i.e., TD elements) to populate with the gradient values.
/
/ Return Value:
/
/ This function populates an empty array passed to it with the gradient elements and it returns the maximum absolute effective gradient value. (The
/ gradient element with largest effective magnitude can be used to determine convergence (i.e., stopping criteria).) Effective here means it takes into
/ account the boundary values since we have a constrained optimization problem.
/------------------------------------------------------------------------------------------------------------------------------------------------------------*/
using System;
using System.Diagnostics;
using Accord.Statistics.Distributions.Univariate;
using OptimalEquity.Helpers;

namespace OptimalEquity.Model
{
    public static class BldGrad
    {
        public static double Run(string type, double[] prms, double[] a, int fxTD, double rf0, long n, int nbuckets,
            int prec, double stdpnr, long npnr, double alpha, int plproc, double[] gradvctr)
        {
            int[] prtls = {-1, -1, -1, -1};
            double[] k = new double[fxTD];
            double maxval = 0;
            NormalDistribution normdist = new NormalDistribution(0.00, 1.00);

            // Iterate over each time point deriving each gradient entry.
            Trace.WriteLine("");
            Trace.Write("Building gradient ");
            for (int g = 0; g < fxTD; ++g)
            {
                // Construct the constant needed for gradient entries.
                k[g] = Funcs.vp(prms, a[g]) / (2.00 * Funcs.v(prms, a[g])) +
                       Math.Pow(Funcs.mp(prms), 2) / (2.00 * Funcs.vp(prms, a[g]));

                // Populate the gradient vector for this time point.
                prtls[0] = g;
                double grdpnr = GetPNR.Run(type, prms, a, fxTD, rf0, n, nbuckets, prec, prtls, plproc);
                gradvctr[g] = k[g] * (grdpnr - stdpnr);
                Trace.Write(".");

                // Maximum effective absolute value of this gradient vector.
                if (a[g] + gradvctr[g] > 1.00)
                {
                    if (maxval < 1.00 - a[g])
                        maxval = 1.00 - a[g];
                }
                else if (a[g] + gradvctr[g] < Funcs.mva(prms) + 0.0001)
                {
                    if (a[g] - (Funcs.mva(prms) + 0.0001) > maxval)
                    {
                        maxval = a[g] - (Funcs.mva(prms) + 0.0001);
                    }
                }
                else if (Math.Abs(gradvctr[g]) > maxval)
                {
                    maxval = Math.Abs(gradvctr[g]);
                }
            }
            Trace.WriteLine(" (Done)");

            // Print the unadjusted gradient vector entries.
            Trace.WriteLine("");
            Trace.Write("Gradient (no adjustment):");
            WrtAry.Run(-1, gradvctr, "Grd", fxTD);

            // If using simulation, test each element for equality with zero. If test holds set the element to zero.
            if (type == "sim" && alpha < 1.00)
            {
                maxval = 0;
                for (int g = 0; g < fxTD; ++g)
                {
                    double cmbpnr = (n * (stdpnr + gradvctr[g] / k[g]) + npnr * stdpnr) / (n + npnr);
                    double ts = (stdpnr + gradvctr[g] / k[g] - stdpnr) /
                                Math.Sqrt(cmbpnr * (1.00 - cmbpnr) * (1.00 / n + 1.00 / npnr));
                    double pval = 2.00 * Math.Min(normdist.DistributionFunction(ts),
                                      1.00 - normdist.DistributionFunction(ts));
                    if (pval > alpha)
                        gradvctr[g] = 0.00; // The element is not different from zero at significance level alpha.

                    // Maximum effective absolute value of this gradient vector.
                    if (a[g] + gradvctr[g] > 1.00)
                    {
                        if (maxval < 1.00 - a[g])
                            maxval = 1.00 - a[g];
                    }
                    else if (a[g] + gradvctr[g] < Funcs.mva(prms) + 0.0001)
                    {
                        if (a[g] - (Funcs.mva(prms) + 0.0001) > maxval)
                            maxval = a[g] - (Funcs.mva(prms) + 0.0001);
                    }
                    else if (Math.Abs(gradvctr[g]) > maxval)
                    {
                        maxval = Math.Abs(gradvctr[g]);
                    }
                }

                // Print the adjusted gradient vector entries.
                Trace.WriteLine("");
                Trace.Write("ADJUSTED gradient:");
                WrtAry.Run(-1, gradvctr, "Adj-Grd", fxTD);
            }

            // Display the maximum effective absolute value of the gradient vector (used to determine convergence).
            Trace.WriteLine("");
            Trace.WriteLine($"Maximum effective absolute value of this gradient: {maxval:F10}");

            // This function returns the maximum absolute value of the gradient elements.
            // (Which is used to define the stopping/convergence criteria.)
            return maxval;
        }
    }
}
