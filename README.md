# random-loop-toolbox

This repository features code for experimenting with random monodromy loops for linearly-parametrized families of polynomial systems. To use the main package "MonodromyMixing', you should have [Macaulay2](http://www2.macaulay2.com/Macaulay2/) vesion 1.10 or higher. The plots are made using [R](https://www.r-project.org/).

To get started, clone this repo, start a Macaulay2 session from the cloned directory, install the package from the source directory:

        git clone https://github.com/timduff35/random-loop-toolbox.git
        cd random-loop-toolbox/
        M2
        installPackage "MonodromyMixing"

To use the package in a future Macaulay2 session, issue the following command:
        needsPackage "MonodromyMixing"

The experiments provided in "bezout.m2" and "powerFlow.m2" give a sense of what can be done. 

The main functions exported by "MonodromyMixing" can be used as follows:

1. simulate(P,p0,sols)
  * P should be a linearly-parametrized family of polynomial systems. In Macaulay2, this may be represented by an object of class "PolySystem" over a ring whose variables are the parameters defining the family, as in the examples provided in the MonodromySolver [documentation](http://www2.macaulay2.com/Macaulay2/doc/Macaulay2..*1.11/share/doc/Macaulay2/MonodromySolver/html/_solve__Family.html).
  * p0, an object of class ["Point"](https://faculty.math.illinois.edu/Macaulay2/doc/Macaulay2..*1.10/share/doc/Macaulay2/NAGtypes/html/___Point.html), gives an initial base point in the parameter space for seeding the solver. 
  * sols, a list of all solutions to P specialized at p0.
  * The most useful options are "MixLimit" and "Iterations." In each iteration, we start from a random base point in the parameter space and collect statistics obtained by concatenating "MixLimit" many loops based at this point. 
  * returns a mutable list with elements of the form (Number of loops concatenated, "OK"/"NA", permutation), where the permutation is represented as a mutable list and "NA" indicates a numerical error (ie. path..*jumping.)
2. simulate(P,p0,sols,filename)
  * same as above, but with a string "filename", with output to be written to "filename.csv"
3. simulateModel(n,iters,mixLimit)
  * Simulates a random walk on the symmetric group S_n for comparison purposes. Currently, either a stataionary (ie uniform) random walk (triggered by setting option RandMethod=>"Unif") or a multiStep random walk with binomial distribution B(n,p). 


Preferred (but not foolproof) way to get reasonable "p0" & corresponding "sols" from P:
        needsPackage "MonodromySolver"
        (V,npaths)=monodromySolve P
        p0=V.BasePoint
        sols=points V.PartialSols

Numerous other functions are currently exported, but undocumented. This may change as the project progresses.