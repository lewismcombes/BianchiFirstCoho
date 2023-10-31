ToHecke:= recformat <
  der               : ModTupFld,
  inn               : ModTupFld,
  dim               : RngIntElt,
  coeff_dim         : RngIntElt,
  ta                : Mtrx,
  tb                : Mtrx,
  tu                : Mtrx,
  weight            : SeqEnum,
  level             : RngOrdIdl,
  field             : FldQuad,
  id_index          : RngIntElt,
  d                 : RngIntElt,
  char              : RngQuadIdl,
  det_twists        : SeqEnum,
  PL                : SetIndx,
  r                 : UserProgram,
  chi               : GrpDrchNFElt  >;



/* given level ideal J and weight (k,l), we compute the first cohomology for PSL(2,Z_K) with coefficients in
coinduced-module tensor V_(k,l) */


/* PSL(2,Z[-2])  */

DIM:=function(level_data,k,l,a,b)



level:=level_data`level;
PA:=level_data`projmatA;
PAi:=level_data`projmatAi;
PB:=level_data`projmatB;
PU:=level_data`projmatU;
PUi:=level_data`projmatUi;
weight:=[k,l,a,b];
d:=level_data`d;
char:=level_data`char;
PL:=level_data`PL;
r:=level_data`r;

K<z>:=QuadraticField(-d);
O:=Integers(K);


A,Ai,B,U,Ui,J:=StandardMats(d);


if char eq 0*O then
  F:=K;
else
  F<t>:=ResidueClassField(char);
end if;


if level ne 1*O then
TA:=RecursiveMatrix(PA,ModuleMat(char,weight,A));
TAi:=RecursiveMatrix(PAi,ModuleMat(char,weight,Ai));
TB:=RecursiveMatrix(PB,ModuleMat(char,weight,B));
TU:=RecursiveMatrix(PU,ModuleMat(char,weight,U));
TUi:=RecursiveMatrix(PUi,ModuleMat(char,weight,Ui));
else
TA:=ModuleMat(char,weight,A);
TAi:=ModuleMat(char,weight,Ai);
TB:=ModuleMat(char,weight,B);
TU:=ModuleMat(char,weight,U);
TUi:=ModuleMat(char,weight,Ui);
end if;


E:=space_mat_psl(d,TA,TAi,TB,TU,TUi);

DER:=Kernel(E);

ID:=TA^0;
F:=HorizontalJoin([ID-TA,ID-TB,ID-TU]);
INN:=Image(F) meet DER;

dimension:=Dimension(DER)-Dimension(INN);

t,N:=IsPrincipal(level);

Data:=rec< ToHecke |
der:=DER,
inn:=INN,
dim:=dimension,
coeff_dim:=#Rows(TA),
ta:=TA,
tb:=TB,
tu:=TU,
level:=level,
weight:=[k,l,a,b],
field:=K,
id_index:=PLIndex(K,N),
d:=d,
char:=char,
det_twists:=[a,b],
PL:=PL,
r:=r,
chi:=level_data`chi
>;

return Data;
end function;
