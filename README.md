# random-loop-toolbox

This repository features code for experimenting with random monodromy loops for linearly-parametrized families of polynomial systems. To use the main package "MonodromyMixing', you should have [http://www2.macaulay2.com/Macaulay2/](Macaulay2) vesion 1.10 or higher. The plots are made using [https://www.r-project.org/](R).

The main functions exported by "MonodromyMixing" can be used as follows:

1. simulate(P,p0,sols)
--* P should be a family of polynomial systems---in Macaulay2, this is an object of class "PolySystem" over a ring whose variables are the parameters defining the family, as in the examples provided in the MOnodromySolver [http://www2.macaulay2.com/Macaulay2/doc/Macaulay2-1.11/share/doc/Macaulay2/MonodromySolver/html/_solve__Family.html](documentation).
--* p0, an object of class [https://faculty.math.illinois.edu/Macaulay2/doc/Macaulay2-1.10/share/doc/Macaulay2/NAGtypes/html/___Point.html]("Point"), gives an initial base point in the parameter space for seeding the solver. 
--* sols, a list of initial solutions.
--* The most useful options are "MixLimit" and "Iterations." In each iteration, we start from a random base point in the parameter space and collect statistics obtained by concatenating "MixLimit" many loops based at this point. 
--* returns a mutable list with elements of the form (Number of loops concatenated, "OK"/"NA", permutation), where the permutation is represented as a mutable list and "NA" indicates a numerical error (ie. path-jumping.)
2. simulate(P,p0,sols,filename)
--* same as above, but with a string "filename", with output to be written to "filename.csv"
3. simulateModel(n,iters,mixLimit)
--* Simulates a random walk on the symmetric group S_n for comparison purposes. Currently, either a stataionary (ie uniform) random walk (triggered by setting option RandMethod=>"Unif") or a multiStep random walk with binomial distribution B(n,p). 

Numerous other functions are currently exported, but un-documented. This may change as the project progresses.