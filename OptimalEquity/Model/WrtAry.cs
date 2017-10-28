/*
/ Copyright (c) 2015 Chris Rook
/
/ Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except in compliance with the License. You may obtain a copy of
/ the License at http://www.apache.org/licenses/LICENSE-2.0. Unless required by applicable law or agreed to in writing, software distributed under the License
/ is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the License for the specific language
/ governing permissions and limitations under the License.
/
/ Filename: WrtAry.cpp
/
/ Function: WrtAry()
/
/ Summary:
/
/ This function will write the contents of an array as a block to either standard output or to standard output and a file specified by the user.
/
/ Parameters:
/
/ 1.) The success probability for the glide-path being printed, if it applies. (Set to something not between 0.00 and 1.00 when printing gradient.)
/ 2.) The array being printed to standard output or a file specified as a pointer. (Will be either the glide-path or the gradient.)
/ 3.) The text to use for labeling each element of the array. (For example, "GP" for GP[] or "Grd" for Grd[].)
/ 4.) The number of elements in the array to be printed. (This will be the # of time-points for the arrangement being analyzed.)
/ 5.) An optional filename. If specified the array is printed the file specified, which is placed in the directory specified by the user.
/
/ Return Value:
/
/ This function has no return value.
/------------------------------------------------------------------------------------------------------------------------------------------------------------*/
using System;
using System.Diagnostics;
using System.IO;
using System.Text;

namespace OptimalEquity.Model
{
    public static class WrtAry
    {
        public static void Run(double stdpnr, double[] a, string lbl, int fxTD, string filenm = " ")
        {
            // Declare and initialize local variables.
            int maxcols = 5;
            int nrows = fxTD % maxcols == 0 ? fxTD / maxcols : fxTD / maxcols + 1;

            // Write glide-path to output window by default.
            Trace.WriteLine("");
            if (stdpnr >= -0.000001 && stdpnr <= 1.000001)
            {
                Trace.Write("--> Success probability for this Glide-Path=");
                Trace.WriteLine(stdpnr);
            }
            for (int r = 0; r < nrows; ++r)
            {
                for (int c = 0; c < maxcols; ++c)
                {
                    if (r + nrows * c < fxTD)
                    {
                        //C++ TO C# CONVERTER TODO TASK: The cout 'setfill' manipulator is not converted by C++ to C# Converter:
                        //ORIGINAL LINE: cout << lbl + "[" << setfill('0') << setw(2) << noshowpos << r + nrows *c << "]=" << right << setw(10) << showpos << a[r + nrows *c] << " ";
                        //C++ TO C# CONVERTER TODO TASK: The cout 'noshowpos' manipulator is not converted by C++ to C# Converter:
                        //C++ TO C# CONVERTER TODO TASK: The cout 'showpos' manipulator is not converted by C++ to C# Converter:
                        Trace.Write(lbl + $"[{r + nrows * c:D2}]={(a[r + nrows * c] > 0 ? "+" : "")}{a[r + nrows * c]:F10}");
                        Trace.Write(" ");
                    }
                }
                Trace.WriteLine("");
            }

            // If a filename is given, write to file.
            if (filenm != " ")
            {
                var sb = new StringBuilder("\n");
                if (stdpnr >= -0.000001 && stdpnr <= 1.000001)
                    sb.AppendLine($"--> Success probability for this Glide-Path = {Math.Round(stdpnr, 12)}");
                for (int r = 0; r < nrows; ++r)
                {
                    for (int c = 0; c < maxcols; ++c)
                    {
                        if (r + nrows * c < fxTD)
                            sb.Append($"{lbl}[{r + nrows * c:D2}]={a[r + nrows * c]:F10} ");
                    }
                    sb.AppendLine();
                }
                File.WriteAllText(filenm, sb.ToString());
            }
        }
    }
}
