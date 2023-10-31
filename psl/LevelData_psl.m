ToSpace:= recformat <
  level             : RngOrdIdl,
  projmatA          : List,
  projmatAi         : List,
  projmatB          : List,
  projmatU          : List,
  projmatUi         : List,
  projmatJ          : List,
  d                 : RngIntElt,
  char              : RngQuadIdl,
  PL                : SetIndx,
  r                 : UserProgram,
  chi               : GrpDrchNFElt >;

/* PGL(2,Z[-d]) level data  */


LEVELDATA:=function(level,d,char,chi)

K:=QuadraticField(-d);
O<z>:=MaximalOrder(K);

char:=char*O;

if char eq 0*O then
  F:=K;
else
  char:=Factorization(char)[1,1];
  F<t>:=ResidueClassField(char);
end if;

PL,r:=ProjectiveLine(quo<O|level>);

A,Ai,B,U,Ui,J:=StandardMats(d);

DG:=DirichletGroup(level);
chi:=Elements(DG)[chi];


if level ne 1*O then
PA:=ProjMat(F,level,A,PL,r,chi);
PAi:=ProjMat(F,level,Ai,PL,r,chi);
PB:=ProjMat(F,level,B,PL,r,chi);
PU:=ProjMat(F,level,U,PL,r,chi);
PUi:=ProjMat(F,level,Ui,PL,r,chi);
else
PA:=[**];
PAi:=[**];
PB:=[**];
PU:=[**];
PUi:=[**];
end if;



Data:=rec< ToSpace | level:=level, projmatA:=PA, projmatAi:=PAi, projmatB:=PB, projmatU:=PU, projmatUi:=PUi, d:=d, char:=char, PL:=PL, r:=r, chi:=chi  >;

return Data;
end function;
