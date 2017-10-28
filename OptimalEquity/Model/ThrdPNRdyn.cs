/*
/ Copyright (c) 2015 Chris Rook
/
/ Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of
/ the License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by applicable law or agreed to in writing, software distributed under the License
/ is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language
/ governing permissions and limitations under the License.
/
/ Filename: ThrdPNRdyn.cpp
/
/ Function: ThrdPNRdyn()
/
/ Summary:
/
/ This function serves as a wrapper for the PNRdyn() function and it splits the job of computing a ruin probability with a dynamic program across the
/ various independent processing units on the computer running the application. This function is passed the number of independent processing units on
/ the computer and it divides the DP discretization by this number and launches simultaneous calls to PNRdyn() to process specific sets of buckets
/ concurrently populating a pre-defined array at these bucket indices. When invoking this function expect a single success probability to be returned.
/
/ Parameters:
/
/ 1.) Six member long double array of parameter settings for stock mean, stock variance, bond mean, bond variance, stock-bond covariance, expense ratio.
/ 2.) Array of long doubles holding the glide-path to use for this ruin calculation (as pointer). The corresponding array must have exactly TD elements.
/ (Note that we only consider alphas between the low-volatility portfolio and 1.00. If any alpha provided in this array is not between the low-
/ volatility portfolio and 1.00 it is forced to the nearest acceptable value.)
/ 3.) Number of time-points for this implementation. This represents a fixed number of years to consider for the retirement horizon (i.e., 30).
/ 4.) The fixed inflation-adjusted withdrawal rate that the retiree is using during decumulation (i.e., 0.04 for a 4% withdrawal rate).
/ 5.) The number of buckets to use when discretizing the ruin factor dimension. This value equals RFMax*(Discretization Precision), where these 2 values
/ are set by the user in the control file.
/ 6.) A 4-element indicator array that specifies which densities should represent inflation/expense-adjusted returns at the given time-point.
/ The settings for this array are:
/ a.) prtls[] = {-1,-1,-1,-1} for all standard real/expense adjusted return densities.
/ b.) prtls[] = {i,-1,-1,-1} for gradient element "i" which uses gradient density g() for the real/expense-adjusted return at the i-th
/ time-point.
/ c.) prtls[] = {i,j,-1,-1} for the off-diagonal Hessian elements which uses the gradient density g() at both time-points "i" and "j".
/ d.) prtls[] = {-1,-1,i,-1} for the diagonal Hessian elements which uses h1() as the real/expense-adjusted return at time-point "i".
/ e.) prtls[] = {-1,-1,-1,i} for the diagonal Hessian elements which uses h2() as the real/expense-adjusted return at time-point "i".
/ 7.) The number of independent processing units on the computer running the application or alternatively the number of parallel processes to use during
/ execution as specified by the user. (If user specified it is set as a parameter when invoking the application.)
/ 8.) The number of buckets to use when discretizing each ruin factor unit value. For example, if this value=1000 then each ruin factor unit is
/ represented by 1000 buckets. This value is specified by the user in the control file when type="dp".
/
/ Return Value:
/
/ This function returns the success probability derived using a DP.
/------------------------------------------------------------------------------------------------------------------------------------------------------------*/
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading.Tasks;
using Accord.Statistics.Distributions.Univariate;
using OptimalEquity.Helpers;

namespace OptimalEquity.Model
{
    public static class ThrdPNRdyn
    {
        public static double Run(double[] prms, double[] a, int fxTD, double rf0, int nbuckets, int[] partls, int plproc, int prec)
        {
            // Declare local variables.
            double[] vp = new double[nbuckets];
            double[] v = new double[nbuckets];
            double[] mts = new double[7];
            List<int> uBkts = new List<int> {1, nbuckets};
            List<double> c;
            List<GammaDistribution> g;
            // Initiate prior timepoint's probabilities with zeros.
            for (int b = 1; b <= nbuckets; ++b)
            {
                vp[b - 1] = 0.00;
            }
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
            // Check that all equity ratios are between MVA and 1.00.
            //=========================================================
            for (int y = 0; y < fxTD; ++y)
            {
                if (a[y] < Funcs.mva(prms) + 0.0001 || a[y] > 1.00)
                {
                    if (a[y] > 1.00)
                        a[y] = 1.00;
                    else if (a[y] < Funcs.mva(prms) + 0.0001)
                        a[y] = Funcs.mva(prms) + 0.0001;
                }
            }

            // Iterate over all time points, launching separate threads to process equal sized collections of
            // buckets concurrently within each time point.
            for (int y = fxTD - 1; y >= 1; y--)
            {
                // Populate moments array for this time point.
                mts[0] = Funcs.m(prms, a[y]);
                mts[1] = Funcs.mp(prms);
                mts[2] = Funcs.v(prms, a[y]);
                mts[3] = Funcs.vp(prms, a[y]);
                mts[4] = Funcs.vpp(prms);
                mts[5] = Funcs.mva(prms);
                mts[6] = Funcs.kh1(prms, a[y]);

                // Process collections of buckets concurrently.
                var t = new Task[plproc];
                c = GetConst.Run(mts, partls, y);

                // Define needed distributions.
                g = GetGamma.Run(partls, y);

                // Update prnbkt and derive the new # of buckets to process per run for the next iteration.
                // (Buckets with trivial assignments are handled separately.)
                var cont = 1;
                var pwr = 1;
                int prnbkt = nbuckets;
                if (y < fxTD - fxTD / 6)
                {
                    prnbkt = 1;
                    pwr = 2;
                }
                for (int b = prnbkt + (int) Math.Pow(-1.00, pwr) * 2 * plproc;
                    cont == 1 && b >= 1 && b <= nbuckets;
                    b = b + (int) Math.Pow(-1.00, pwr) * (2 * plproc))
                {
                    int[] prnbkts = {b, b, nbuckets};
                    PNRdyn.Run(mts, prec, prnbkts, vp, v, uBkts, c, g);
                    if (pwr == 1 && v[b - 1] >= 1.00)
                        prnbkt = b;
                    else if (pwr == 2 && v[b - 1] < 1.00)
                        prnbkt = Math.Min(b + (int) Math.Pow(-1.00, pwr) * 2 * plproc, nbuckets);
                    else
                        cont = 0;
                }

                // The value of prnbkt should be ge plproc and le nbuckets.
                if (prnbkt > nbuckets)
                    prnbkt = nbuckets;
                else if (prnbkt < plproc)
                    prnbkt = plproc;
                int bktsprun = prnbkt / plproc + 1;
                if (bktsprun * plproc > nbuckets)
                    bktsprun = nbuckets / plproc;

                // If # buckets/run is still not right, exit with an error and fix.
                if (bktsprun * plproc > nbuckets || bktsprun < 1)
                {
                    Trace.WriteLine("");
                    Trace.WriteLine("The # of buckets per run is not being derived correctly, must fix:");
                    Trace.Write("# of buckets per run = ");
                    Trace.WriteLine(bktsprun);
                    Trace.Write("Total # of buckets = ");
                    Trace.WriteLine(nbuckets);
                    Trace.Write("# of concurrent processes being used = ");
                    Trace.WriteLine(plproc);
                    Trace.WriteLine("EXITING...ThrdPNRdyn()...");
                    Console.Read();
                    Environment.Exit(1);
                }

                // An array of size 3 is used to specify the start/end/total buckets for each thread.
                var bktarys = new int[plproc][];
                for (int i = 0; i < plproc; ++i)
                {
                    bktarys[i] = new int[3];
                    bktarys[i][0] = bktsprun * i + 1;
                    bktarys[i][1] = Math.Min(bktsprun * (i + 1), nbuckets);
                    bktarys[i][2] = nbuckets;

                    var j = i;
                    var cLocal = c.ToList();
                    var gLocal = g.ToList();
                    t[i] = Task.Run(() =>
                    {
                        PNRdyn.Run(mts, prec, bktarys[j], vp, v, uBkts, cLocal, gLocal);
                    });
                }

                // Assign trivial (known) values.
                if (bktarys[plproc - 1][1] < nbuckets)
                {
                    for (int b = bktsprun * plproc + 1; b <= nbuckets; ++b)
                    {
                        v[b - 1] = 1.00;
                    }
                }
                // Wait for all threads to finish, then proceed.
                Task.WaitAll(t);

                // Free temporary memory allocations and reused containers.
                for (int i = 0; i < plproc; ++i)
                {
                    bktarys[i] = null;
                    bktarys[i] = null;
                }
                uBkts.Clear();
                c.Clear();
                g.Clear();

                // Update quantities needed for next iteration.
                // Reset the prior year's probabilities in Vp[] and derive the unique bucket quantities.
                cont = 1;
                var prevprob = 0.00;
                for (int b = 1; b <= nbuckets; ++b)
                {
                    vp[b - 1] = v[b - 1];

                    // Verify the integrity of the probabilities derived during this iteration.
                    if (vp[b - 1] < prevprob - 1e-15 || vp[b - 1] > 1.00 + 2.00 * Math.Pow(0.1, 16))
                    {
                        Trace.WriteLine("");
                        Trace.WriteLine($"There is an issue with the probabilities derived at this timepoint (t={y}), see below:");
                        if (b > 1)
                        {
                            Trace.Write("Vp[");
                            Trace.Write(b - 2);
                            Trace.Write("] = ");
                            Trace.WriteLine(vp[b - 2]);
                        }
                        Trace.Write("Vp[");
                        Trace.Write(b - 1);
                        Trace.Write("] = ");
                        Trace.WriteLine(vp[b - 1]);
                        Trace.WriteLine("EXITING...ThrdPNRdyn()...");
                        Console.Read();
                        Environment.Exit(1);
                    }
                    prevprob = vp[b - 1];

                    // Derive the new vector of bucket #'s with unique probabilities at the prior timepoint.
                    if (b == 1 || b != nbuckets && Math.Abs(v[b - 1] - v[b]) > 1e-15 && cont == 1 || (b == nbuckets))
                    {
                        uBkts.Add(b);
                    }

                    // Once PRuin=1 stop collecting unique buckets.
                    if (cont == 1 && vp[b - 1] >= 1.00)
                    {
                        cont = 0;
                    }
                }

                // Check that the last bucket at this time point has a PRuin of 1.00.
                // (Otherwise RFMax needs to be increased. This must hold when using special densities also.)
                if (vp[nbuckets - 1] < 1.00)
                {
                    Trace.Write("Timepoint (t=");
                    Trace.Write(y);
                    Trace.Write("), has V[");
                    Trace.Write(nbuckets - 1);
                    Trace.Write("]=");
                    Trace.Write(v[nbuckets - 1]);
                    Trace.WriteLine(", which is < 1.00.");
                    Trace.WriteLine("(Increase RFMax.)");
                    Trace.WriteLine("EXITING...ThrdPNRdyn()...");
                    Console.Read();
                    Environment.Exit(1);
                }
            }

            // Process final timepoint only for the given RF0.
            mts[0] = Funcs.m(prms, a[0]);
            mts[1] = Funcs.mp(prms);
            mts[2] = Funcs.v(prms, a[0]);
            mts[3] = Funcs.vp(prms, a[0]);
            mts[4] = Funcs.vpp(prms);
            mts[5] = Funcs.mva(prms);
            mts[6] = Funcs.kh1(prms, a[0]);
            int rf0Bkt = (int)(rf0 * prec + 0.5);
            c = GetConst.Run(mts, partls, 0);
            g = GetGamma.Run(partls, 0);
            int[] fnlbkts = { rf0Bkt, rf0Bkt, nbuckets };
            PNRdyn.Run(mts, prec, fnlbkts, vp, v, uBkts, c, g);

            // Retrieve the probability to return.
            double rtprob = 1.00 - v[rf0Bkt - 1];

            // The single PNR derived using a DP is returned.
            return rtprob;
        }
    }
}
