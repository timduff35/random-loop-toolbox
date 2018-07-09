-*
TODO:
-) separate residual, correspondence, and other tolerances + defaults
-) refine  
*-

newPackage(
        "MonodromyMixing",
        Version => "0.1", 
        Date => "June 23, 2018",
        Authors => {{Name => "Tim Duff", 
                  Email => "tduff3@gatech.edu", 
                  HomePage => "http://people.math.gatech.edu/~tduff3/"}},
        Headline => "Package useful for running experiments related to monodromy solving of polynomial systems",
	PackageImports => {"PHCpack"},
	PackageExports => {"NumericalAlgebraicGeometry","MonodromySolver"},
        DebuggingMode => true
        )

export {"ResidualTolerance","solsToPerm","MixLimit","simulate","trans","cycles",
        "sparseFamily","denseFamily","viewPerms","Refine","writeSimulations","simulateModel",
	"unifRList","RandMethod"}

---GLOBAL VARIABLES
numericalErrorMessage="failed correspondence--is the residual tolerance set realistically?"
defaultTol = 0.1--  should it be 

-- PHCpack configuration will depend on your machine
--loadPackage ("PHCpack", Configuration=>{"path"=>"~/","PHCexe"=>"./phc"},Reload=>true)


-*
Methods for constructing random points
*-

-- is this the current random type implemented in M2? if so, not needed
old'random'Type = lookup(random,Type)
random Type := o -> R -> (
    if class R === ComplexField then (
	exp(2 * pi * random RR * ii)
	) 
    else (old'random'Type o) R
    ) 

-- box mueller transform for standard normal RVs
gaussCC = () -> (
    (u1,u2):=(random RR,random RR);
    sqrt(-2*log(u1))*cos(2*pi*u2)+ii*sqrt(-2*log(u1))*sin(2*pi*u2)
    )

-- samples from the unit sphere in C^n
sphere = n -> (
	v:=apply(n,i->gaussCC());
	(1/norm(v))*v
	)

--sphere with standard normal noise
noiseySphere = n -> (
	v:=sphere(n);
	w:=apply(n,i->gaussCC());
	v+w
	)

-- random number on complexu nit circle returned 
gam = () -> (
    exp(ii * random RR)
    )

-*
Methods for constructing random systems
*-
	    
-- random dense polynomial
randDensePoly = (R,d) -> (
    t:=symbol t;
    S:=CC[gens R|{t}];
    sub(sub(random(d,S),{t=>1}),R)
    )

-- random dense polynomial
denseFamily = (R,d) -> sparseFamily polySystem apply(numgens R,i->randDensePoly(R,d))

-- B: a matrix encoding the sparsity pattern of an unmixed system
randSparsePoly = B -> B*random(CC^(numcols B),CC^1)

-- sparse
sparseFamily = method(Options=>{})
sparseFamily PolySystem := o -> PS -> (
    polys := flatten entries PS.PolyMap;
    ind := flatten apply(#polys,i-> -- indices for parameters                                              
    	apply(exponents polys#i, t->(i,t))                                                                 
    	);                                                                   
    R := PS.PolyMap.ring;  
    W:=symbol W;-- what if the user wants to use a different symbol?                 
    AR := CC[apply(ind,i->W_i)][gens R];                                                                   
    polysP := for i to #polys-1 list -- system with parameteric coefficients and same support              
    sum(exponents polys#i, t->(W_(i,t))_AR*AR_(t));	
    polySystem transpose matrix {polysP}
    )


-*
--generic unmixed system over random sparse basis
createBaseSystem = method(Options=>{Sparsity=>1.0,Hypersurface=>false})
createBaseSystem (Ring,ZZ) := o -> (R,d) -> (    
    B:=first coefficients(randDensePoly(R,d));
    s:=ceiling((numcols B)*o.Sparsity);
    B=B_(take(random toList(0..numcols B-1),{0,s-1}));
    if o.Hypersurface then return polySystem(apply(n-1,i->random CC+ random(1,S))|{randSparsePoly B});
    return polySystem apply(n,i->randSparsePoly B)
    )

-- symmetric distance between permuted items in a transposition
dist = Q -> (
    T:=positions(apply(Q,i->(i-Q#i)),k->k!=0);
	if length(T) != 2 then error("not a transposition");
    min((first T-last T)%diam,(last T - first T)%diam)
    )
*-


-*
Methods for manipulating permutations and computing various statistics
*-

-- extract cycle decomposition of permutaiton represented as a list
-- returns mutable hash table
cycles = L -> (
    --if not(ancestor(class L,List)) then return "NA";
    copyPerm:=new MutableList from L;
    cycleTable:=new MutableHashTable from {};
    for i from 0 to #L-1 do (
	k:=i;
	counted:=false;
	if copyPerm#k == -1 then counted = true;
	if not counted then (
	    newCycleLength:=0;
	    if copyPerm#k == k then (
	    	copyPerm#k=-1;
	    	newCycleLength=1;
	    	);
	    while copyPerm#k != -1 do (
	    	copyPerm#k=-1;
	    	k=L#k;
	    	newCycleLength=newCycleLength+1;
	    	);
	    copyPerm#i=-1;
	    if not member(newCycleLength,keys cycleTable) then cycleTable#newCycleLength = 0;
	    cycleTable#newCycleLength = cycleTable#newCycleLength +1;
	    );
	);
    cycleTable
    )

-- average cycle length
avgCyc = method()
avgCyc MutableHashTable := H -> sub(sum apply(keys(H),i->i*H#i) / sum(values(H)),RR)

-- min cycle length
minCyc = method()
minCyc MutableHashTable := H -> min keys H

-- max cycle length
maxCyc = method()
maxCyc MutableHashTable := H -> max keys H

-- median cycle length
medCyc = method()
medCyc MutableHashTable := H -> (
    (L,n,i,j):=(sort keys H,sum values H,0,0);
    if length L == 1 then return L#0;
--    if length L = 2 then return L#0
    while j<n/2 do (
	j=j+H#(L#i);
	i=i+1;
	);
    if i==#L then return last L;
    if (sum values H)%2==1 then (
    	if j > n/2 then return L#(i-1)
    	else return L#i
	);
    sub((L#(i-1)+L#i)/2,RR)
    )
    

-- compute the length given cycle hash table
trans=method()
trans MutableHashTable := H -> (
    if not(ancestor(class H,MutableHashTable)) then return "NA";
    sub(sum(apply(keys(H),i->(i-1)*H#i)),RR)
    )
trans String := s -> s

-- compose permutations
pCompose = method()
pCompose (MutableHashTable, MutableHashTable) := (H1, H2) -> (
    H:=new MutableHashTable from {};
    apply(keys H1,k-> if H2#?(H1#k) then k=> H2#(H1#k))
    )
pCompose (List, List) := (L1,L2) -> (
	apply(length(L1),i-> L2#(L1#i))
		)
pCompose (MutableList, MutableList) := (L1,L2) -> (
	new MutableList from apply(#L1,i-> L2#(L1#i))
		)
pCompose (String,MutableList) := (s,L) -> s	    
pCompose (MutableList,String) := (L,s) -> s	    
pCompose (String,String) := (s,s') -> s

--- methods for random permutations

coin = p -> (
    r:=random RR;
    if r <= sub(p,RR) then return 1
    else return 0
    )

bin = (n,p) -> sum apply(n,i->coin(p))

unifRPerm = n -> (
    L:=new MutableList from 0..(n-1);
    for i from 0 to (n-2) do (
	j:=random(i,n-1);
	temp:=L#j;
	L#j=L#i;
	L#i=temp;
	);
    L
    )    


--array representing a random transposition
randomTransposition = n -> (
    L':=new MutableList from 0..(n-1);
    i:=random(0,n-1);
    j:=random(0,n-1);
    L'#i=j;
    L'#j=i;
    L'
    )

binTrans = (n,p) -> (
    newPerm:=new MutableList from 0..(n-1);
    for k from 0 to bin(n,p)-1 do (
	newPerm=pCompose(randomTransposition n,newPerm);
	);
    newPerm
    )


-*
initializeChain = method(Options=>{ResidualTolerence=>0.01,RootCount=>false,SparseBasis=>false,MSOptions=>new OptionTable from {}})
initializeChain PolySystem := o -> P -> (
    if instance(coefficientRing ring p,PolynomialRing) then L:=solveFamily(P,MSOptions)
    else if SprarseBasis then L=sparseMonodromySolve(P,MSOptions)
    else L=solveSystem P,x->x.Coordinates;
    L=apply(L,x->x.Coordinates);
    if not checkResidual(L,P,o.ResidualTolerance,rc) or ((o.RootCount != false) and (#L == o.RootCount)) then return false;
    L
    )
initializeChain (PolySystem, OptionTable) := o -> (MSOptions, P) -> (
    L:=solveFamily(P,MSOptions)
    )

-- we also need methods for families
singleStep = method(Options=>{ResidualTolerence=>0.01,SparsityLower=>0.0,SparsityUpper=>1.0,ResidualTolerance=>0.01,Family=>false})
singleStep (PolySystem, PolySystem, List) := o -> (Pbase,Pstep,L) -> (
    track((random CC)*Pbase,(random CC)*Pstep,L)
    )
*-

checkSolutions = (L,P,eps,rc) -> (
    residualCheck:=all(apply(L,p->max flatten entries evaluate(P,p)<eps),s->s);
    residualCheck and (#L==rc)
    )

-- In: L, the original solutions list
-- S', solutions obtained by following
-- Out: array representation of a random monodromy loop, status
solsToPerm = method(Options=>{ResidualTolerance=>defaultTol})
solsToPerm (List,List)  := o -> (Lhash,Shash) -> (
    l:=length Lhash;
    perm:=new MutableList from toList(l:-1);
    for k from 0 to  l-1 do perm#k=minPosition apply(Shash,s->norm(s - Lhash#k));
    if sort(toList perm)==toList(0..l-1) then return (perm,"OK") else return (perm,"NA")
    )

TEST ///
setRandomSeed 1
R=CC[x,y]
F={x^5+y^2-1,x+y-2}
L=solveSystem F
G={random(5,R)+random(4,R)+random CC,random(1,R)+random CC}
S'=track((random CC)*G,(random CC)*F,track((random CC)*F,(random CC)*G,L))
assert(toList first solsToPerm(L,S')=={0,2,1,3,4})
///

-*
Methods for simulation
*-


simulateModel = method(Options=>{RandMethod=>"mix"})
simulateModel (ZZ,ZZ,ZZ,RR) := o -> (n,iters,mixLimit,p) -> (
    L:=new MutableList from 0..(n-1);
    M:=new MutableList from  (iters*mixLimit):(0,"OK",L);
    for i from 0 to iters-1 do (
    	for j from 0 to mixLimit-1 do (
	    if o.RandMethod=="mix" then newPerm := binTrans(n,p) else newPerm = unifRPerm n;
	    if j==0 then M#(i+j*iters)=(j,"OK",newPerm) else M#(i+j*iters)=(j,"OK",pCompose(newPerm,last M#(i+(j-1)*iters)));
	    );
	);
    M
    )


simulate = method(Options=>{TargetSolutionCount=>null,ResidualTolerance=>defaultTol,
	                    MixLimit=>1,Iterations=>2,Refine=>true})
simulate (PolySystem,Point,List) := o -> (P,p0,sols) -> (
    --assumes dominant family!
    rc:=#sols;
    M:=new MutableList from (o.Iterations*(o.MixLimit)):(o.Iterations,"OK",toList(0..#sols-1));
    m:=numgens coefficientRing ring P;
    initialSystem:=specializeSystem(p0,P);
    if o.Refine then refine(initialSystem,sols);
    for i from 0 to o.Iterations-1 do (
	print("iteration " | i | ", initialize!\n");
	-- TODO: this "goodInit, goodLoop" stuff can be abstracted
--	print(toExternalString "bad Initialization?" | "\n");
	baseSystem:=specializeSystem(point random(CC^m,CC^1),P);
	L:=track((random CC)*initialSystem,(random CC)*baseSystem,sols);
	if o.Refine then refine(baseSystem,L);
	hashMatrix:=random(CC^(length (first L).Coordinates),CC^1);
    	Lhash:=apply(L,p->((matrix p)*hashMatrix)_(0,0));
	for j from 0 to (o.MixLimit-1) do (
	    print("mix step!\n");
	    G:=specializeSystem(point random(CC^m,CC^1),P);   
	    S'':=track((random CC)*baseSystem,(random CC)*G,L);
	    if o.Refine then refine(G,S'');
	    S':=track((random CC)*G,(random CC)*baseSystem,S'');
	    if o.Refine then refine(baseSystem,S');
	    Shash:=apply(S',p->((matrix p)*hashMatrix)_(0,0));
	    (perm,isNa):=solsToPerm(Lhash,Shash);
	    if j==0 then perm = new MutableList from perm else perm=pCompose(perm,last M#(i+(j-1)*o.Iterations));
	    M#(i+j*o.Iterations)=(j,isNa,perm);
	    );
	);
    M
    )
simulate (PolySystem,Point,List,String) := o -> (P,p0,sols,filename) -> (
    M:=simulate(P,p0,sols,o);
    writeSimulations(M,filename);
    M
    )

writeSimulations = (M,filename) -> (
    file := openOut (currentFileDirectory | filename | ".csv");--maybe we should create a data directory?
    file << "isNA,Steps,Length,AvgLen,MedLen,MaxLen,MinLen" << endl;
    for i from 0 to (#M-1) do (
	C:=cycles toList last M#i;
	file << M#i#1 << "," << first M#i << "," << trans C << ","  << avgCyc C << "," << medCyc C << "," << maxCyc C << "," << minCyc C << endl;
	);
    close file;
    )


    

viewPerms = M -> (
    if #(M#0) > 1 then M:=apply(M,x->toList last x); 
    netList apply(toList M,x->x)
    )

TEST ///
setRandomSeed 0
R=CC[a,b,c,d][x,y]
P=polySystem {a*x^4+b,c*y^4+d}
V = first monodromySolve P
iters=10
mixlim=3
M=simulate(P,V.BasePoint,points V.PartialSols,Iterations=>iters,MixLimit=>mixlim)
toList String := s -> s
toList apply(M,x->if instance(last x,BasicList) then trans cycles toList last x else last x)
assert(#M == (iters+1)*(mixlim+1))
x///

end

uninstallPackage packageName
restart
packageName = "MonodromyMixing"
installPackage packageName

-* quadric experiment from before
setRandomSeed 0
(n,d)=(2,5)
(iters,steps)=(1000,8)
R=CC[x_1..x_n]
P=denseFamily(R,d)
m=numgens coefficientRing ring P
p0=point random(CC^m,CC^1)
sols=solveSystem specializeSystem(p0,P)
-- small example
elapsedTime M=simulate(P,p0,sols)
netList toList M
-- generate data file
setRandomSeed 0
elapsedTime M= simulate(P,p0,sols,"bezout"|n|"vars"|d|"deg"|iters|"iters"|steps|"steps",Iterations=>iters,MixLimit=>steps)
N=simulateModel(25,iters,steps,14.69/25)
writeSimulations(N,"bezoutModel")
N=simulateModel(25,iters,steps,14.69/25,RandMethod=>"Unif")
writeSimulations(N,"bezoutUnif")

*-

needs "~/Workshop-2018-Leipzig/NumericalAG/examples/powerNetwork.m2"
n=5
G=completeGraph n
P=powerEquations G
setRandomSeed 1
(V,npaths)=monodromySolve(P)
(iters,steps)=(1000,8)
elapsedTime N=simulate(P,V.BasePoint,sols,"k5",Iterations=>iters,MixLimit=>steps)

N=simulateModel(54,iters,steps,46.7/54)
writeSimulations(N,"k5Model")

N=simulateModel(54,iters,steps,46.5/54,RandMethod=>"Unif")
writeSimulations(N,"k5Unif")






-- Wnt experiment: should probably just cache the results!
needs("../dev/M2/M2/Macaulay2/packages/MonodromySolver/paper-examples/example-traceCRN.m2")

P=mSys
--product apply(equations P,e->first degree e) Bezout bound is about 1Bil
p0=W.BasePoint
sols=points W.PartialSols
--Pp0=polySystem specializeSystem(p0,P)
--max apply(sols,x->norm evaluate(polySystem PtraceWitness,x))




setRandomSeed 0
elapsedTime N=simulate(mSys,p0,sols,Iterations=>10,MixLimit=>1,Refine=>true)
apply(toList N, x->length unique toList last x)
viewPerms N

