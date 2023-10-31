
AttachSpec("/mnt/data/smp19lmc_code/ArtinAlgebras/ArtinAlgebras.spec");
ChangeDirectory("/mnt/data/smp19lmc_code/thesis/H1/");

/*
ChangeDirectory("/home/smp19lmc/Documents/work/maths/code/modular forms/H1");
AttachSpec("/home/smp19lmc/Documents/work/maths/code/modular forms/ArtinAlgebras/ArtinAlgebras.spec");
*/

load "ProjAction.m";
load "ModuleAction.m";
load "RecursiveMatrix.m";
load "Utility.m";


load "LevelData_pgl.m";
load "space_mat.m";
load "Space_pgl.m";

load "Hecke_poly.m";
load "Hecke.m";
load "Utility_Hecke.m";
load "GetFormVals.m";

load "Periods.m";



//
// Change things here
//


//field of definition
d:=1;
K:=QuadraticField(-d);
O<z>:=MaximalOrder(K);

//level, weight, twists (TBA), number of eigenvalues, dirichlet character
N:=-4*z + 3;
level:=N*O;
wt:=[0,0];
tw:=[0,0];
chi:=1;
HB:=40;
char:=0;


//
//
//



LD:=LEVELDATA(level,d,char,chi);
space:=DIM(LD,wt[1],wt[2],tw[1],tw[2]);

print "Full dimension:", space`dim;
print "Cocycle dimension:", Dimension(space`der);
print "Coboundary dimension:", Dimension(space`inn);


print "Eisenstein dimension:", CuspNumber_SL(K,level);


der:=space`der;
inn:=space`inn;

mforms,pi:=quo<der|inn>;
g:=Inverse(pi);


if char eq 0 then
  HP:=[TP : TP in PrimesUpTo(HB,K) | GCD(TP,level) eq 1*O];
  //HP:=[TP : TP in PrimesUpTo(HB,K)];
else
  HP:=[TP : TP in PrimesUpTo(HB,K) | GCD(TP,level*char) eq 1*O];
end if;
HNF:=[HNF_basis(J): J in HP];
ParallelSort(~HNF,~HP);


time HH:= GetHeckeMatrices(space,HP);


EV_systems,form_vals:=GetFormVals(space,HH,HP,g);


MFPeriods:=function(n)
  for j in [1..#form_vals[n,3]] do
    if GCD(form_vals[n,3][j]*O,level) eq 1*O then
      PeriodAverage(space,form_vals[n],j);
    end if;
  end for;
  return n;
end function;


EV_systems;

//
