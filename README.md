# MOCHA

MOCHA stands for _Matroid Optimization: Combinatorial Heuristics and Algorithms_.
MOCHA provides heuristics and algorithms to solve problems in the burgeoning field of multicriteria matroid optimization.
Besides specific algorithms and heuristics, matroid data structures are included in order to provide a foundation for old and new algorithms.
MOCHA was initially developed for the paper "Computation in Multicriteria Matroid Optimization", accepted by the [Journal of Experimental Algorithmics](http://www.jea.acm.org/), 2009, and authored by J. De Loera, D. Haws, J. Lee, and A. O'Hair.

## Project Manager

Dr. David Haws,
University of Kentucky,
Department of Statistics,
871 Patterson Office Tower,
Lexington, KY 40506-0027

## Installation
```
git clone --branch stable/1.0 https://github.com/coin-or/MOCHA.git
cd MOCHA
./configure
make
make install
```

### Third-party software dependencies

GMP, LAPACK, BLAS. MOCHA can compile without GMP, though some functionality
will be lost such as rank calculation using arbitrary precision arithmetic.

More information is available in the INSTALL file.


## Matroid Optimization

Matroids encapsulate the combinatorial notion of independence. They can be
found in graphs, matrices, point sets, hyperplane arrangements and many more.
There are many ways to define a matroid, but for MOCHA the most useful
axiomatization is that of bases. Given a finite set S (ground set) with n
elements, the collection of subsets (B) of S are the bases for a matroid M if

1) every element of (B) has the same cardinality, 
2) if B1 and B2 are in (B), for all x in B1 not in B2, there exists y in B2 
not in B1 such that B1 - x + y is in (B).

Perhaps the simplest matroid example is given by a graph. Given a graph G, the
bases of the matroid on G are all the spanning trees of G. It is not hard to
see that the collection of spanning trees of G all have the same cardinality
and the exchange property 2) holds.

Another intuitive example is given by a matrix A. The ground set is the labeled
columns of A, and a base of the matroid on A is any maximally linear
independent collection of columns of A. ( Mat[roid/rix] )

For more on matroids see: "Matroid Theory", Welsh, 1976, "Matroid Theory",
Oxley, 1992 or for a sufficient online coverage [http://en.wikipedia.org/wiki/Matroid].

Matroid optimization appears when weights are placed on the ground set S. For
example, by placing a single real valued weight on every element of S, every
base B1 has a weight given by the sum of the weights it contains. One could
then optimize over the bases of a matroid by asking for the base with maximal
(minimal) weight. This problem can be solved by the strongly polynomial greedy
algorithm (which MOCHA has a function for).

MOCHA can tackle an even larger family of optimization problems. Instead of
placing a single real value weight on each element of S, one places d real
valued weights on every element of S. In this case, instead of a base B1 having
a real valued weight, the weight of B1 is a vector in R^d^, given by the sum of
the weights it contains. To ask for an optimal base, given multiple weightings,
one needs a way to decide between two vectors in R^d^. There are multiple
methods to distinguish two points in R^d^, such as a functional f: R^d^ --> R, 
Pareto inequality, or min-max. We call these balancing functions, and they may
be linear, convex, or highly non-linear. 

The difficulty of matroid optimization (with multiple weightings) lays in the
fact that there are exponentially many bases with respect to n, the size of the
ground set. So, one can not simply iterate over all bases of a matroid and
optimize. Fortunately, under some assumptions on the weightings applied to the
ground set S, the number of points obtained by considering the weights of all
bases is polynomially bounded. There are deterministic algorithms to compute
the projected bases (projected by the weighting on S), but they are
computationally _heavy_. See "Nonlinear Matroid Optimization and Experimental
Design", Berstein et. al., 2008.

MOCHA is a software packaged developed in conjunction with the paper
"Computation in Multicriteria Matroid Optimization", submitted to the ACM
Journal of Experimental Algorithmics, 2009. Its goals are to provide fast
heuristics and algorithms aimed at solving multiple weight matroid optimization
problems. Typically, the choice of weightings on the ground set ensure that the
number of projected bases is much smaller than all bases. Many of MOCHA's
heuristics attempt to give all the projected bases (the most difficult part),
leaving the final optimization over a balancing function for later.

The software in MOCHA will take in as input a matroid (currently supported
formats are matrices, graphs, uniform) and the weightings on the ground
set.  It will output a subset of the projected bases (guaranteed to be all
projected bases in some cases).


## How MOCHA works
Data Format:
    The input format for MOCHA is a text file in the following format:
```
        <MATROID TYPE>
        <MATROID DESCRIPTION>
        ...
        <Number of weightings/criteria>
        <Number of matroid elements>
        <#> <#>  ... <#>
        <#> <#>  ... <#>
        ...
        <#> <#>  ... <#>
```

Example (Instances/Calibration/gn9e18d2w0w20.mo):
```
    GRAPH
    ADJACENCY
    9
    9
    1   0   0   0   0   1   0   1   0
    0   1   0   0   1   1   1   1   0 
    0   0   1   0   0   0   0   1   0 
    0   0   0   1   1   0   0   1   0 
    0   1   0   1   1   1   1   1   0 
    1   1   0   0   1   1   1   1   1 
    0   1   0   0   1   1   1   1   1 
    1   1   1   1   1   1   1   1   1 
    0   0   0   0   0   1   1   1   1 
    2
    18
    18   1   7   1   8   7  14  17   6   9  20   6  10   3   8   5   4  16 
    1  18   0   4  16  17  19   1  18   5   1  11   4  17  10   9  19  16 
```

MATROID TYPE can be either "GRAPH" or "VECTOR".
If the type is "GRAPH" then the next line should be "ADJACENCY".
Following this the adjacency representation of the graph should follow.
It should be in the format:
```
  <Number of rows>
  <Number of columns>
```
followed by an appropriate row and column of numbers.

If the type is "VECTOR" then the following lines should be a matrix given as
```
    <Number of rows>
    <Number of columns>
```
followed by an appropriate row and column of numbers.

After the matroid type is specified, the weight is given as a matroid given as
```
    <Number of rows>
    <Number of columns>
```
followed by an appropriate row and column of numbers.

An example of a vector matroid representation is:
```
    VECTOR
    3
    10
    0   6   8   7   0  10   6  10   9  10 
    3   2   0   6  10   3   5   1   5   5 
    4   5   4   9   1   0   1   7  10   6 
    2
    10
    1   7   1   5   7   8   3   7   5   2 
    3   8   7   7   0   6   5   3   8   5 
```

### What MOCHA does

The most important program is 'matroidtest' which contains heuristics/algorithms
which output subsets of the projected bases. The projected bases are 
outputted in matlab format. The file should also be easy to parse using GNU
utilities. We hope to add more output formats.

The three general methods to find subsets of the projected bases are:

#### DFBFS: Different Fiber BFS

This is a modified version of breadth first search which only pivots to new basis
as long as it is a new projected basis.

#### Pivot Test

This finds a tight rectangular box containing all projected basis then proceeds
to use local search or tabu search on each point (many times if specified), using
a specialized convex function. All the projected bases found are output.

#### Boundary

This outputs the vertices boundary of the convex hull of the projected basis.
Our heuristic also has the side effect of outputting some extreme points besides
the vertices.

For further description of our algorithms heuristics, see below.

## Instances from our Paper

The directory 'Instances/' contains all of our input data for the paper as well as
some simple examples demonstrating our software.

### Detailed How-To/Example

    To reproduce most of the tests in our paper, only three programs
    are needed 'matroidtest', 'localsearch', and 'tabusearch'. Below
    we give instructions how to reproduce experiments shown in our
    paper.

    The program 'matroidtest' has the capability of running the 
    following algorithms: Different Fiber BFS, Pivot Test, Projected 
    Boundary Calculation, Pareto Calculation, and Brute force 
    projected bases calculation (Graphs only). See our paper for
    a full description of these algorithms. 'matroidtest' has the 
    option to output all computed points in a MATLAB ready file
    which will plot the points.

    The program 'matroidtest' can be run with no arguments. It is
    interactive and will ask which algorithms to run. If more than
    one are selected, e.g. DFBFS and Boundary calculation, then
    all the points will be output to the same file, with different
    labels.

    'matroidtest' can also be invoked as:
```
        matroidtest <inputfile>
```
        or 
```
        matroidtest <inputfile> <outputfile>
```


### Different Fiber BFS

        Run 'matroidtest'. Specify the input file either in the command line or
        when prompted. Next it will ask "Different Fiber BFS Enumerations?
        (y/n)". Enter "y". It will ask for the number of searches. This is the
        number of DFBFS searches to perform on the interior and boundary.
        'matroidtest' will perform this number of searches using (pseudo)random
        projected boundary points. Next, 'matroidtest' starts at a random
        unvisited projected base a number of times specified above.  Next
        'matroidtest' will ask for the BFS depth. This is the recursion level
        that DFBFS will terminate at.  Following this, 'matroidtest' will ask
        for the "Interior Random retry limit?".  This is the limit of how many
        times DFBFS will try to find a new random projected bases before
        giving up. 'matroidtest' will then ask for the "Boundary random retry
        limit?". This is the number of times DFBFS will try to find a new
        random projected base on the boundary before giving up. When prompted
        for all other algorithms/test, enter "n". When prompted for the output
        file, give a name. Alternatively enter the output filename as the
        second command-line argument.

        Here is a quick example:

        Run the following command from the src/ directory: 
```
        ./matroidtest ../Instances/Calibration/gn11e20d2w0w100.mo tmpout.m
```

        Enter the following sequence of input:
```
        y
        100
        8
        10000
        100
        n
        n
        n
        n
        n
        n
        n
```

        This will run DFBFS 100 times, with a truncation level of 8, an
        interior retry limit of 10000, and a boundary retry limit of 100.
        Typically it takes MOCHA longer to find a random projected boundary
        point than an interior point.  It will write the outputed points to
        tmpout.m. In MATLAB execute 'tmpout.m'.
    
### Local Search/Tabu Search

        The objective for local/tabu search must be hard coded into the files
        localsearch.cpp and tabusearch.cpp and recompiled. By default the
        object function is set to x^2. We did not want to spend excessive time
        implementing or including a bulky format to read in arbitrary objective
        functions. To change the objective function, study the function "MyFct"
        in localsearch.cpp or tabusearch.cpp.

        localsearch and tabusearch are fully interactive and will ask for all
        the necessary arguments. Alternatively they can be invoked as 
```
localsearch <inputfile> <outputfile> <number of searches>

tabusearch <inputfile> <outputfile> <number of searches> <tabu searchlimit>
```
        Again, the output files are written in MATLAB format which can
        be run to display graphically the local and tabu searches. 

        Example from within /src/ directory:
```
./localsearch ../Instances/Calibration/gn11e20d2w0w20.mo tmp.m 100
```
        This will use gn11e20d2w0w20.mo as input, output to the file tmp.m
        and run 100 local searches.

### Auto Box Pivot Heuristic
   
        This will perform the Pivot Test described in our paper. Run 
        'matroidtest', specify input and output files and when asked
        "Run Auto Box Pivot Heuristic (Local Search)? (y/n)" or
        "Run Auto Box Pivot Heuristic (Tabu Search)? (y/n)" select
        y according to which pivot heuristic you want to use. 'matroidtest'
        will then ask for the number of times to run local search or
        tabu search on each point. For tabu search, 'matroidtest' will 
        also ask for the tabu search limit. The auto box pivot heuristic
        will then automatically calculate an appropriate bounding box of
        the projected bases using 2(number of weightings) local searches.

        An example run (for localsearch) would be as follows:

        Run the command from the src/ directory: 
```
./matroidtest ../Instances/Calibration/gn11e20d2w0w100.mo tmpout.m
```
        Enter the following sequence of input:
```
        n
        y
        20
        n
        n
        n
        n
        n
        n
```
        This will run the auto box pivot heuristic using local search using
        20 attempts for each point.

### Box Pivot Heuristic 

        This is much the same as the Autobox heuristic above in 'matroidtest'
        except that the user is prompted for coordinates of a box of points
        to test using Pivot Test.

### Boundary Calculation

        When there are only two criteria/weightings 'matroidtest' can 
        calculate the boundary (and more) of the projected bases. When asked
        "Run Boundary Calculation (dim=2 only)? (y/n)" enter "y".

### Pareto Optimum Test

There are four options for finding Pareto optima:
1. Boundary Only,
2. Boundary Triangular Region Search,
3. BFS Pareto Search, and
4. BFS and Boundary Triangular Region Pareto Search.

#### Boundary Only: This will compute the boundary of the projected bases.
            'matroidtest' will then output the Pareto optima that are contained
            in the boundary.

#### Boundary Triangular Region Search: This follows the Boundary and 
            Triangular Region Pareto Test heuristic described in our paper.
            'matroidtest' will compute the boundary, find regions where the
            remaining Pareto optima will be, then run the Pivot Heuristic using
            either local search or tabu search to enumerate a subset of 
            potential Pareto optima.

#### BFS Pareto Search: This will perform Different Fiber BFS and from the
            returned set, compute the Pareto optima. Note that these are only 
            guaranteed to be Pareto optima of the points found by DFBFS, which
            may not necessarily be Pareto optima of the original problem.

#### BFS and Boundary Triangular Region Pareto Search: This will perform
            DFBFS and Boundary Triangular Region Search and combine the 
            results.

### Brute Force Enumeration:

        'matroidtest' can also enumerate all projected bases using reverse search
        when the matroid is graphical. This uses the algorithm of Matsui.


## Executables contained in the MOCHA package

- `alltoproj`
    This file is meant to convert the output of findChildren printing all
    spanning trees and project them by some given weighting.  Expects input to
    be sets given by unsigned ints on each line.  First argument should be the
    matrix that is our weighting.  The second argument should be the number of
    elements in each set

- `designtomatrix`
    This program will read in two matrices; a n x k exponent vector matrix B, 
    a m x k design point matrix P

    It then computes A_ij := \prod_{h=1}^k P_{i,h}^B_{j,h} 
    It reads in the matrices from standard input and prints the matrix to 
    standard output

- `estimatebases`
    This program estimates the number of bases of a matrix by assigning
    random weights to each column of the matrix and finding the maximal
    weight basis. This is done m times and the resulting average estimates
    the function GAMMA. Then upper and lower bounds for the number of bases
    are calculated.

- `genmatrix`
    Generates a random matrix. Useful for creating random examples.
    run genmatrix for input format.

- `genrandmo`
   This is an interactive program used to generate random Vector and Graph 
   matroids along with random weightings.

- `graphtest`
    Program to check validity of some internal graph routines.

- `localsearch`
    Program to perform localsearch heuristic as described in our paper.
    See explanation above for usage.

- `matroidtest`
    The program 'matroidtest' has the capability of running the 
    following algorithms: Different Fiber BFS, Pivot Test, Projected 
    Boundary Calculation, Pareto Calculation, and Brute force 
    projected bases calculation (Graphs only). See our paper for
    a full description of these algorithms. 

- `mvbalclust`
    Non-linear matroid optimization can solve the problem of min variance
    balanced clustering. The input is a file that is a matrix of
    2k points in R^n^ which are to be partitioned into two sets of size
    k such that the variance of their euclidean distances is minimized.
    Input is a file with the following format.
    ```
    <Number of rows>
    <Number of columns>
    ```
    followed by an appropriate row and column of numbers.

- `nagibatest`
    This program was written to test the our implementation of Matsui
    and Nagamochi-Ibariki algorithms. The input is 
    ```
    ADJACENCY
    <Number of rows>
    <Number of columns>
    ```
    followed by an appropriate row and column of numbers.
    The previous matrix is the adjacency representation of the graph.

    nagibatest will enumerate all spanning trees and only print out
    the total number of trees. The code can be modified easily
    to print out all trees to stdout. Set "printTrees = 1;".

- point2po
    This program reads in from stdin points, and finds the pareto optimum using 
    a straightforward search

- tabusearch
    Program to perform localsearch heuristic as described in our paper.
    See explanation above for usage.

- testmatrix
    Program to internally check matrix class.


## Data structures and code:
    Effort was made to make the code legible and intuitive. The authors
    hope to continue to add functionality and subroutines to many of
    the objects. Others are welcome and encouraged to contribute! 
