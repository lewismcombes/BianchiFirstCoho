/* **********************************************************************
   ** hecke.m                                                          **
   **                                                                  **
   ** Computing Bianchi modular forms and the Hecke action             **
   **                                                                  **
   ** Mehmet Haluk Sengun & Lewis Combes                               **
   **                                                                  **
   ********************************************************************** */

/* TP is the prime ideal for the Hecke operator */

HECKE:=function(TP,DIMdata)

  JJ:=DIMdata`level;
  coeff_size:=DIMdata`coeff_dim;
  char:=DIMdata`char;
  weight:=DIMdata`weight cat DIMdata`det_twists;


  DER:=DIMdata`der;
  INN:=DIMdata`inn;

  d:=DIMdata`d;
  K:=QuadraticField(-d);
  O<z>:=MaximalOrder(K);

  A,Ai,B,U,Ui,J:=StandardMats(d);
  t,s:=IsPrincipal(TP);
  D:=Matrix(O,2,2,[s,0,0,1]);

  if char eq 0*O then
    F:=K;
  else
    F<t>:=ResidueClassField(char);
  end if;

  PD:=ProjMat(F,JJ,D,DIMdata`PL,DIMdata`r,DIMdata`chi);
  TD:=RecursiveMatrix(PD,ModuleMat(char,weight,D));

  // time saved here by computing coset reps and permutations once
  cc:=CosetMats(K,TP);

  CosetPerm:=function(reps,mat,I) //I is the level defining reps
    perm:=[];
    for R in reps do
      M:=R*mat;
      for i in [1..#reps] do
        conj:=M*reps[i]^(-1);
        if conj[2,1] mod I eq 0 and conj in Parent(M) then
          Append(~perm,i);
        end if;
      end for;
    end for;
    return perm;
  end function;

  permA:=CosetPerm(cc,A,TP);
  permB:=CosetPerm(cc,B,TP);
  permU:=CosetPerm(cc,U,TP);
  permJ:=CosetPerm(cc,J,TP);

  SSS:=MatrixAlgebra(F,coeff_size);   /* this algebra is used for the polynomials */

  PRR<X,Y,Z,W>:=PolynomialRing(SSS,4);  /* X,Y,Z,W represent f(A), f(B), f(U), f(J) coefficients are matrices */

  poly_A:=PRR!0;
  poly_B:=PRR!0;
  poly_U:=PRR!0;
  poly_J:=PRR!0;

  for i in [1..Norm(TP)+1] do
  /* we will get the summation that defines the Hecke operator's image on A,B, U and J */
  poly_A:= PRR!poly_A + PRR!POLY(A,permA,cc,TP,i,DIMdata,TD,d);
  poly_B:= PRR!poly_B + PRR!POLY(B,permB,cc,TP,i,DIMdata,TD,d);
  poly_U:= PRR!poly_U + PRR!POLY(U,permU,cc,TP,i,DIMdata,TD,d);
  poly_J:= PRR!poly_J + PRR!POLY(J,permJ,cc,TP,i,DIMdata,TD,d);
  end for;

  H11:=MonomialCoefficient(poly_A,X);
  H21:=MonomialCoefficient(poly_A,Y);
  H31:=MonomialCoefficient(poly_A,Z);
  H41:=MonomialCoefficient(poly_A,W);

  C1:=VerticalJoin([H11,H21,H31,H41]);

  H12:=MonomialCoefficient(poly_B,X);
  H22:=MonomialCoefficient(poly_B,Y);
  H32:=MonomialCoefficient(poly_B,Z);
  H42:=MonomialCoefficient(poly_B,W);

  C2:=VerticalJoin([H12,H22,H32,H42]);

  H13:=MonomialCoefficient(poly_U,X);
  H23:=MonomialCoefficient(poly_U,Y);
  H33:=MonomialCoefficient(poly_U,Z);
  H43:=MonomialCoefficient(poly_U,W);

  C3:=VerticalJoin([H13,H23,H33,H43]);

  H14:=MonomialCoefficient(poly_J,X);
  H24:=MonomialCoefficient(poly_J,Y);
  H34:=MonomialCoefficient(poly_J,Z);
  H44:=MonomialCoefficient(poly_J,W);

  C4:=VerticalJoin([H14,H24,H34,H44]);


  H0:=HorizontalJoin([C1,C2,C3,C4]);

  R,f:=quo<DER|INN>;

  B:=Basis(R);

  pB:=[];
  g:=Inverse(f);

  for b in B do
    Append(~pB,g(b));
  end for;

  Mat:=[];

  for b in pB do
    Append(~Mat,ElementToSequence(f(b*H0)));
  end for;

  return Matrix(F,Mat);
end function;
