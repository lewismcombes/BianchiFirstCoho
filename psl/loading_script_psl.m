
AttachSpec("/mnt/data/smp19lmc_code/ArtinAlgebras/ArtinAlgebras.spec");
ChangeDirectory("/mnt/data/smp19lmc_code/thesis/H1/");

load "ProjAction.m";
load "ModuleAction.m";
load "RecursiveMatrix.m";
load "Utility.m";

load "psl/LevelData_psl.m";
load "space_mat.m";
load "psl/Space_psl.m";

load "psl/Hecke_poly_psl.m";
load "psl/Hecke_psl.m";
load "Utility_Hecke.m";
load "GetFormVals.m";

load "psl/Periods_psl.m";



//
// Change things here
//


//field of definition
d:=2;
K:=QuadraticField(-d);
O<z>:=MaximalOrder(K);

//level, weight, twists (TBA), number of eigenvalues
N:=2+5*z;
level:=N*O;
wt:=[0,0];
tw:=[0,0];
chi:=1;
HB:=30;
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



HP:=[TP : TP in PrimesUpTo(HB,K) | GCD(TP,level) eq 1*O];
//HP:=[TP : TP in PrimesUpTo(HeckeBound,K)];
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
