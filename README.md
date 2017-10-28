# Optimal Equity Glidepaths in Retirement

***NOTE: the project is done for Lucas Weatherill (OnTrack Retirement).***

***The app is based on CHRISTOPHER J. ROOK C++ implementation from [the doc](doc/1506.08400-1.pdf).***

## How to run

Compile and run 

    # OptimalEquity.exe [<optional parameter: number of threads>]
    
Input files are created automatically with default values.  
A single parameter is optionally accepted which is the number of tasks the user wants to process concurrently as it runs. If no value is supplied the code determines the maximum number of concurrent processing units on the machine running it and uses this value. The maximum number of concurrent processing units is generally the number of CPU cores.

## Configure

There are two input files that are created automatically in "./config/" directory if do not exist.

### control.txt

Example:

<pre>
0.082509 0.0402696529 0.021409 0.0069605649 0.0007344180 0.000
30 0.05 0.00000000001
(nr or ga) (dp or sim) <parameters>
</pre>

The 1st line of the control file contains the real stock and bond means, variances, covariance and the expense 
ratio in the following order: μ<sub>s</sub>, σ<sup>2</sup>
<sub>s</sub>, μ<sub>b</sub>, σ<sup>2</sup><sub>b</sub>, σ<sub>(s,b)</sub>, and, E<sub>R</sub>. This setup file uses 
historical returns. 

The 2nd line of the control file contains 
the horizon length, withdrawal rate, and epsilon convergence level as: T<sub>D</sub>, W<sub>R</sub>, and, ε. 
(Note: The first 2 lines will always contain the terms just described.) 

The 3rd line of the control file specifies the
optimization scheme (“nr” vs. “ga”), the estimation/approximation method (“dp” vs. “sim”),
followed by settings specific to the estimation/approximation method.  

When a DP is used (dp) the 2 settings are the discretization level (PR) and the maximum ruin factor used during discretization
(RF<sub>Max</sub>). 

When simulation is used (sim) there are 3 additional parameters to specify in the following order: N, alpha-level 1, alpha-level 2.

For details see page 18 and 19 from [the original doc](doc/1506.08400-1.pdf).

### gp.txt

From [the original doc](doc/1506.08400-1.pdf): The procedure requires an initial glidepath be specified and optimization begins at this
point. If the surface being optimized has multiple local optimums then the starting glidepath
determines where the procedure ends. In Section IV.E we found that all 5 starting glidepaths
converge to the same optimum glidepath in Scenarios 1-8. The same was found for Scenarios 9
and 10 in Section V, using different starting glidepaths. Contents of “gp.txt” using Random
Glidepath #1 from Figure 5 is shown next. We do not include blank lines or extra spaces in this
file, and we start on line #1. Therefore this file has exactly 30 lines.

```
0.636
0.214
0.193
0.637
0.626
0.597
0.943
0.877
0.254
0.823
0.903
0.294
0.444
0.513
0.529
0.160
0.564
0.293
0.698
0.228
0.311
0.776
0.689
0.764
0.596
0.793
0.911
0.624
0.709
0.205
```

## Output

2 files as result:
* *MinimizeRuinProbability.exe.log* - log file, all you have on console is saved there.

* *output.txt* - is written to the directory that contains the control
file and the initial glidepath. A sample of this file is shown below. GP[00] is α<sub>1</sub>, GP[01] is α<sub>2</sub>,
etc... The first equity ratio is set at time t=0 and the corresponding return is observed at time t=1.

Example:
```
--> Success probability for this Glide-Path = 0.527952155270
GP[00]=+0.8235272966 GP[06]=+0.8590020928 GP[12]=+0.7307821069 GP[18]=+0.6235703362 GP[24]=+0.5436024264
GP[01]=+0.8617503896 GP[07]=+0.8386754281 GP[13]=+0.7107993614 GP[19]=+0.6085468596 GP[25]=+0.5323875470
GP[02]=+0.8850732670 GP[08]=+0.8170174723 GP[14]=+0.6916598009 GP[20]=+0.5942528104 GP[26]=+0.5216839028
GP[03]=+0.8931568162 GP[09]=+0.7949412393 GP[15]=+0.6733802825 GP[21]=+0.5806493950 GP[27]=+0.5114608960
GP[04]=+0.8890061069 GP[10]=+0.7730074379 GP[16]=+0.6559544071 GP[22]=+0.5676980327 GP[28]=+0.5016899088
GP[05]=+0.8765350143 GP[11]=+0.7515547956 GP[17]=+0.6393611079 GP[23]=+0.5553611039 GP[29]=+0.4923442815
```

## How to use

The app calculates max success probability and the optimal static glidepath for it, for the specified input parameters.

First fill *control.txt*, run the app and check *output.txt*.

### Example

Input *control.txt*:

```
0.082509 0.0402696529 0.021409 0.0069605649 0.0007344180 0.000
20 0.05 0.00001
nr dp 500 2.75
```

Results *output.txt*:

```
--> Success probability for this Glide-Path = 0.950775818493
GP[00]=0.2974686786 GP[04]=0.3373336816 GP[08]=0.3915503321 GP[12]=0.4669838487 GP[16]=0.5772678161 
GP[01]=0.3063624255 GP[05]=0.3493274495 GP[09]=0.4080631847 GP[13]=0.4904722346 GP[17]=0.6134931396 
GP[02]=0.3159278463 GP[06]=0.3622919553 GP[10]=0.4260511472 GP[14]=0.5163658358 GP[18]=0.6562485613 
GP[03]=0.3262273557 GP[07]=0.3763271362 GP[11]=0.4455968640 GP[15]=0.5451994954 GP[19]=0.7047500770 
```

