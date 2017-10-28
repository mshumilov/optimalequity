/*
/ Copyright (c) 2015 Chris Rook
/
/ Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of
/ the License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by applicable law or agreed to in writing, software distributed under the License
/ is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language
/ governing permissions and limitations under the License.
/
/ Filename: Climb.cpp
/
/ Function: Climb()
/
/ Summary:
/
/ This function steps in the direction of the gradient until no further progress can be made. When using simulation as the estimation type, the ruin
/ probabilities are subject to sampling error and we test the hypothesis that the updated GP is not inferior to the best performing GP encountered thus
/ far in the iteration. When using simulation, climbing ends when this test is rejected at an alpha level specified by the user (the first alpha in the
/ control file). Since the test is for non-inferiority and alpha is the probability of making a Type I Error (reject Ho when it is true), a small alpha
/ is more likely to result in climbing past the true stopping point at each iteration, whereas a large alpha will do the exact opposite, namely result in
/ stopping too soon. Since the gradient is an expensive computation (in terms of runtime) the alpha decision is an important one. The hypothesis test
/ is:
/ Ho: Probability of ruin using new GP >= Probability of ruin using base GP (base GP is the best performing thus far in the iteration)
/ Ha: Probability of ruin using new GP < Probability of ruin using base GP (base GP is the best performing thus far in the iteration)
/
/ and only applies when simulation is the estimation technique. Simulated probabilities are subject to sampling error which is a function of the sample
/ size, which is taken into account by the test above. Quantities that are estimated using dynamic programming (dp) are subject to approximation error
/ and climbing proceeds until the probability declines from the last value. No hypothesis test is used when the estimation type is dynamic programming.
/ The step size has been set as a function of the largest gradient (absolute) value. The user may archive better performance by changing the step size
/ rules. A small step size leads to longer run times, whereas a larger step size leads to imprecise stopping criteria which must be corrected at a
/ future iteration.
/
/ Parameters:
/
/ 1.) Estimation type: Either dynamic programming (dp) or simulation (sim).
/ 2.) Six member long double array of parameter settings for stock mean, stock variance, bond mean, bond variance, stock-bond covariance, expense ratio.
/ 3.) Array of long doubles holding the glide-path to start climbing at (as pointer). The corresponding array must have exactly TD elements and is
/ updated with new values when the hypothesis test above is not rejected or climbing continues when using DP approximation.
/ 4.) Number of time-points for this implementation. This represents a fixed number of years to consider for the retirement horizon (i.e., 30).
/ 5.) The fixed inflation-adjusted withdrawal rate that the retiree is using during decumulation (i.e., 0.04 for a 4% withdrawal rate).
/ 6.) The sample size to use when simulating the ruin/success probability. (Note that the ruin probability is the number of times ruin occurs divided
/ by the number of trial runs. The success probability is 1 minus the ruin probability.) This parameter is only used when type="sim".
/ 7.) The number of buckets to use when discretizing the ruin factor dimension. This value equals RFMax*(Discretization Precision), where these 2 values
/ are set by the user in the control file. (This parameter is only used when type="dp".)
/ 8.) The number of buckets to use when discretizing each ruin factor unit value. For example, if this value=1000 then each ruin factor unit is
/ represented by 1000 buckets. This value is specified by the user in the control file when type="dp".
/ 9.) The number of independent processing units on the computer running the application or alternatively the number of parallel processes to use during
/ execution as specified by the user. (If user specified it is set as a parameter when invoking the application.)
/ 10.) The largest absolute gradient value used to determine the step size. Larger gradient entries result in smaller step sizes and vice-versa.
/ 11.) The success probability for the glide-path provided which is the starting point for climbing. To proceed in the direction of the gradient the new
/ success probability must be non-inferior to the current value (using the hypothesis test above) when using type="sim" or just greater than when
/ using type="dp".
/ 12.) The sample size used to simulate the success probability passed via the previous parameter. (This parameter is only used when type="sim".)
/ 13.) The gradient vector as a pointer to an array of exactly TD elements. This is the direction of steepest ascent.
/ 14.) The alpha value used to test the hypothesis that the new GP is non-inferior to the best prior GP. We continue climbing until this Ho is rejected.
/ (This parameter is only used when type="sim". This is the first alpha level set in the control file.)
/ 15.) The maximum number of climbing iterations before returning to the invoking program. (For example, use with initial stepping to avoid a lengthy
/ climb that can occur with a poorly selected starting point.)
/
/ Return Value:
/
/ This function returns the optimal success probability achieved while climbing. The updated glide-path is placed directly into the glide-path array
/ that is passed to this function as a pointer in argument #3.
/------------------------------------------------------------------------------------------------------------------------------------------------------------*/
using System;
using System.Diagnostics;
using Accord.Statistics.Distributions.Univariate;
using OptimalEquity.Helpers;

namespace OptimalEquity.Model
{
    public static class Climb
    {
        public static double Run(string type, double[] prms, double[] gpath, int fxTD, double rf0, long n, int nbuckets, int prec, int plproc, double mxgrd, double stdpnr, long npnr, double[] grdnt, double alpha, int nitrs = 0)
        {
            double maxpnr = stdpnr;
            double iter = 1.00;
            int[] prtls = { -1, -1, -1, -1 };
            long maxn = npnr;
            int cont = 1;
            int fstimpr = 0;
            int iindx = 0;
            int tryup = 0;
            int origindx = 0;
            NormalDistribution normdist = new NormalDistribution(0.00, 1.00);

            // Get initial glidepath provided, assign to both original GP array and prior GP array.
            var prevGP = new double[fxTD];
            var origGP = new double[fxTD];
            for (int y = 0; y < fxTD; ++y)
            {
                origGP[y] = prevGP[y] = gpath[y];
            }

            // Define the step size for climbing. The step size depends on the largest gradient element and grows exponentially.
            // This is a heuristic that has worked well and can be modified if desired. Better step sizes can reduce runtimes.
            for (int i = 1; i <= 10; ++i)
            {
                if (mxgrd >= 1.00 / (10.00 * Math.Pow(10.00, i)) && mxgrd < 1.00 / (10.00 * Math.Pow(10.00, i - 1)))
                {
                    iindx = i;
                    iter = Math.Pow(Math.Exp(Math.Log(5.00) / 4.00), iindx);
                }
            }

            // Output details for the current iteration.
            Trace.WriteLine(new String('=', 70));
            Trace.WriteLine($"Iteration step size = {iter:F10}");
            Trace.Write("Trying to improve on success probability = ");
            Trace.WriteLine(stdpnr);
            Trace.WriteLine(new String('=', 70));
            // Climb in the direction of the gradient.
            for (int t = 0; cont == 1 || fstimpr == 0; ++t) // Iterate until no more progress is made
            {
                if (cont == 0 && fstimpr == 0) // If no progress is made reduce step size and try again.
                {
                    // Set the original index value when entering this problematic scenario.
                    if (origindx == 0)
                        origindx = iindx;
                    // Adjust glidepath back one iteration since it failed to improve the probability.
                    for (int y = 0; y < fxTD; ++y)
                    {
                        gpath[y] = prevGP[y]; // Reverse final update, since no progress was made.
                    }
                    if (iindx != 0)
                    {
                    }
                    Trace.WriteLine("");
                    Trace.Write("No Progress Made: Iteration step size changed from ");
                    Trace.Write(iter);
                    if (iindx == 0)
                    {
                        Trace.WriteLine(" to 0.00. (No additional climbing attempts will be made.)");
                        Trace.WriteLine("");
                        Trace.WriteLine("ERROR: No progress can be made, the procedure is stuck. (Step size has been reduced to 0.)");
                        Trace.WriteLine(" You may be operating along the boundary where the process is not well defined or your");
                        Trace.WriteLine(" estimation/approximation precision level is not adequate for your epsilon level.");
                        Trace.WriteLine("");
                        Trace.WriteLine("");
                        Trace.Write("Current Glide-Path: ");
                        WrtAry.Run(maxpnr, gpath, "GP", fxTD);
                        Trace.WriteLine("");
                        Trace.WriteLine("EXITING...Climb()...");
                        Console.Read();
                        Environment.Exit(1);
                    }
                    else if (iindx == 1 && tryup == 5)
                    {
                        iindx = iindx - 1;
                        iter = iter / 2.00;
                    }
                    else if (iindx > 1 && tryup == 5)
                    {
                        if (iindx == origindx + 5)
                            iindx = origindx - 1;
                        else
                            iindx = iindx - 1;
                        iter = Math.Pow(Math.Exp(Math.Log(5.00) / 4.00), iindx);
                    }
                    else if (iindx > 1 && tryup < 5)
                    {
                        iindx = iindx + 1;
                        iter = Math.Pow(Math.Exp(Math.Log(5.00) / 4.00), iindx);
                        tryup = tryup + 1;
                    }
                    Trace.Write(" to ");
                    Trace.Write(iter);
                    Trace.WriteLine(". (Attempting to climb again.)");
                    cont = 1;
                }
                for (int y = 0; y < fxTD; ++y) // Iterate over glide-path and update it
                {
                    prevGP[y] = gpath[y]; // Reset the previous glide-path element
                    gpath[y] = gpath[y] + (iter) * grdnt[y]; // Update each individual glide-path element
                    if (gpath[y] < Funcs.mva(prms) + 0.0001)
                    {
                        gpath[y] = Funcs.mva(prms) + 0.0001; // Stay above MVA and consistent with ThrdPNRdyn() and ThrdPNRsim().
                    }
                    else if (gpath[y] > 1.00)
                    {
                        gpath[y] = 1.00; // Consistent with ThrdPNRdyn() and ThrdPNRsim().
                    }
                }
                double newpnr = GetPNR.Run(type, prms, gpath, fxTD, rf0, n, nbuckets, prec, prtls, plproc);
                Trace.WriteLine("");
                Trace.Write("Base Prob(NR) = ");
                Trace.Write(maxpnr);
                if (type == "sim")
                {
                    Trace.Write(" (N=");
                    Trace.Write(maxn);
                    Trace.Write(")");
                }
                Trace.WriteLine("");
                Trace.Write("New Prob(NR) = ");
                Trace.Write(newpnr);
                if (type == "sim")
                    Trace.Write($" (N={n})");
                else if (newpnr > maxpnr)
                    Trace.Write(" (Better, CONTINUE climbing ...)");
                else
                    Trace.Write(" (Worse, STOP climbing ...)");
                Trace.WriteLine("");

                // If using simulation, conduct a non-inferiority test of the new vs max base GP.
                // =====> Continue to climb if the new GP is at least as good as the max base GP.
                // Otherwise, compare new probability with old and climb while making progress.
                if (type == "sim")
                {
                    double cmbvar = maxpnr * (1.00 - maxpnr) / maxn + newpnr * (1.00 - newpnr) / n;
                    double ts = (newpnr - maxpnr) / Math.Sqrt(cmbvar);
                    double pval = normdist.DistributionFunction(ts);
                    Trace.Write("Test Statistic = ");
                    Trace.WriteLine(ts);
                    Trace.Write("P-Value = ");
                    Trace.Write(pval);
                    Trace.Write(" (Alpha=");
                    Trace.Write(alpha);
                    Trace.WriteLine(")");
                    if (pval > alpha)
                    {
                        fstimpr = 1;
                        Trace.WriteLine("=====> Accept Ho (non-inferiority), CONTINUE climbing ...");
                    }
                    else
                    {
                        cont = 0;
                        Trace.WriteLine("=====> Reject Ho (non-inferiority), STOP climbing ...");
                    }
                    // Update PNR and sample size for base GP.
                    if (newpnr > maxpnr)
                    {
                        maxpnr = newpnr;
                        maxn = n;
                    }
                }
                else if (type == "dp")
                {
                    if (newpnr > maxpnr)
                    {
                        fstimpr = 1;
                        maxpnr = newpnr;
                    }
                    else
                        cont = 0;
                }
                // For lengthy climbing, display the current glidepath at 100 iteration intervals.
                if ((t + 1) % 100 == 0)
                {
                    Trace.WriteLine("");
                    Trace.Write("Current Glide-Path at Iteration: ");
                    Trace.Write(t + 1);
                    WrtAry.Run(newpnr, gpath, "GP", fxTD);
                }
                // Stop when maximum number of iterations has been reached, if specified.
                //=========================================================================
                if (nitrs > 0 && t + 1 == nitrs && cont == 1)
                {
                    Trace.WriteLine("");
                    Trace.WriteLine($"Climbing limit reached at {nitrs} iterations.");
                    cont = 2;
                }
            }

            // Adjust glidepath back one iteration since it failed to improve the probability.
            // (This is only done when climbing failed to improve, not when limit is reached.)
            if (cont != 2)
            {
                for (int y = 0; y < fxTD; ++y)
                {
                    gpath[y] = prevGP[y];
                }
            }
            if (type == "sim")
            {
                Trace.WriteLine("");
                Trace.WriteLine("Resetting the probability (to remove any built-in upward sampling bias) ...");
                maxpnr = ThrdPNRsim.Run(prms, gpath, fxTD, rf0, 2 * n, prtls, plproc);
            }

            // Return the max success probability.
            return maxpnr;
        }
    }
}
