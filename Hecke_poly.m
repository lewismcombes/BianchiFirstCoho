
/* M is the matrix; A,B or U
   TP is the prime ideal of the Hecke operator
   i is the index
   TD is the matrix associated to the commensurator element associated to the Hecke opertor

*/

POLY:=function(M,permM,reps,TP,i,DIMdata,TD,d)

  JJ:=DIMdata`level;

  TA:=DIMdata`ta;
  TB:=DIMdata`tb;
  TU:=DIMdata`tu;
  TJ:=DIMdata`tj;
  coeff_size:=DIMdata`coeff_dim;
  char:=DIMdata`char;
  weight:=DIMdata`weight cat DIMdata`det_twists;


  K:=QuadraticField(-d);
  O<z>:=MaximalOrder(K);
  Z:=Integers();

  A,Ai,B,U,Ui,J:=StandardMats(d);

  if char eq 0*O then
    F:=K;
  else
    F<t>:=ResidueClassField(char);
  end if;

  t,s:=IsPrincipal(TP);
  D:=Matrix(K,2,2,[s,0,0,1]);      /* this is D itself, we will invert it, so it is over K */


  SSS:=MatrixAlgebra(F,coeff_size);
  N:=SSS!0;
  ID:=SSS!1;


  PRR<X,Y,Z,W>:=PolynomialRing(SSS,4);
  IDD:=Matrix(O,2,2,[1,0,0,1]);

  hM:=D*reps[i]*M*reps[permM[i]]^(-1)*D^(-1);

  g:=PRR!0;

  if hM eq IDD then
    return(g);
  end if;
  if hM eq -IDD then
    return(g);
  end if;

  decomp:=WORD(hM,d);

  t_decomp:=[];

  for i in [1..#decomp] do //creates the same decomposition, but with the T matrices
    if decomp[i,1] eq A then
      Append(~t_decomp,MatPow(TA,decomp[i,2]));
    elif decomp[i,1] eq B then
      Append(~t_decomp,MatPow(TB,decomp[i,2]));
    elif decomp[i,1] eq U then
      Append(~t_decomp,MatPow(TU,decomp[i,2]));
    elif decomp[i,1] eq J then
      Append(~t_decomp,MatPow(TJ,decomp[i,2]));
    else
      Append(~t_decomp,"you stop that");
    end if;
  end for;

  for i in [1..#decomp] do //does the fiddly cocycle calculating
    if decomp[i][1] eq A then
      g+:=HeckeMatP(N,TA,decomp[i][2])*&*t_decomp[i+1..#t_decomp]*X;
    elif decomp[i][1] eq B then
      g+:=HeckeMatP(N,TB,decomp[i][2])*&*t_decomp[i+1..#t_decomp]*Y;
    elif decomp[i][1] eq U then
      g+:=HeckeMatP(N,TU,decomp[i][2])*&*t_decomp[i+1..#t_decomp]*Z;
    elif decomp[i][1] eq J then
      g+:=HeckeMatP(N,TJ,decomp[i][2])*&*t_decomp[i+1..#t_decomp]*W;
    else
      g:="you stop that";
    end if;
  end for;

  Rj:=reps[permM[i]];
  PRj:=ProjMat(F,JJ,Rj,DIMdata`PL,DIMdata`r,DIMdata`chi);
  TRj:=RecursiveMatrix(PRj,ModuleMat(char,weight,Rj));

  return MonomialCoefficient(g,X)*TD*TRj*X+
         MonomialCoefficient(g,Y)*TD*TRj*Y+
         MonomialCoefficient(g,Z)*TD*TRj*Z+
         MonomialCoefficient(g,W)*TD*TRj*W;


end function;
