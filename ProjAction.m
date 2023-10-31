///
// For integers x & y with GCD(x,y) = d, gives integers [a,b] such that ax+by = d
//
Bezout:=function(x,y) /* a Bezout algoirthm*/
     Q:=Rationals();
     if y eq 0 then
        return [Sign(Q!x),0];
     end if;
     q1:=x div y;
      r:= x - q1*y;
     if r eq 0 then
          if Norm(q1) ge 0 then
         return [1 , 1 - q1];
     else
         return [1 , -1 - q1];
     end if;
     else
     q2:=$$(y,r);
          return [q2[2],q2[1]-q2[2]*q1];
     end if;
end function;



ProjMat:=function(K,J,M,PL,r,chi)

     O:=CoefficientRing(Parent(M));

     if Type(O) ne RngQuad then // sometimes M is defined over K rather than O, this deals with that case
       O:=MaximalOrder(O);
     end if;

     SS:=MatrixAlgebra(O,2);
     M:=SS!M;
     FacMatrix:=[**];

      proj_action:=function(M,i)
         t,im,sca :=r(PL[i]*M,true,true);
         if t then
           return [*Index(PL,im),sca*];
         else
           return [*-1,0*];
         end if;
      end function;

      perm:=[];
      scalars:=[];
         for i in [1..#PL] do
           PA:=proj_action(M,i);
         Append(~perm,PA[1]);
         Append(~scalars,PA[2]);
      end for;

      l:=#PL;
      Mat:=[];
      for i in [1..l] do
        new_row:=[];
        for j in [1..l] do
          if j eq perm[i] then
            Append(~new_row,chi(scalars[i]));
          else
            Append(~new_row,0);
          end if;
        end for;
        Append(~Mat,new_row);
      end for;

      Append(~FacMatrix,ChangeRing(Matrix(Mat),K));

     Insert(~FacMatrix,1,#PL);

     return FacMatrix;
end function;





//
// Computes coset representative matrices for Gamma0(J) in PSL(Z_K)
//
CosetMats:=function(K,J)

  O<z>:=MaximalOrder(K);
  PL,r:=ProjectiveLine(quo<O|J>);
  t,j:=IsPrincipal(J);

  coset_reps:=[];

  for i in [1..#PL] do
    bottomRow:=[PL[i][1],PL[i][2]];
    while Abs(Norm(GCD(bottomRow[1],bottomRow[2]))) ne 1 do /* moving to a rep with gcd = 1 */
      bottomRow[1]+:=j;
    end while;
    topRow:=Bezout(bottomRow[1],bottomRow[2]);
    newMatrix:=Matrix(O,2,2,[topRow[2],-topRow[1],bottomRow[1],bottomRow[2]]);
    u:=Determinant(newMatrix);
    if u ne 1 then /*this makes the determinant 1 always */
       newMatrix:=Matrix(O,2,2,[newMatrix[1][1]/u,newMatrix[1][2]/u,newMatrix[2][1],newMatrix[2][2]]);
    end if;
    assert Determinant(newMatrix) eq 1;
    Append(~coset_reps,newMatrix);
  end for;

  return coset_reps;
end function;



//
// Computes the matrices t_i(M) and coset representatives for Gamma_0(P), where P is a prime ideal giving the Hecke operator
//
CoresMat:=function(K,P,M)
     O<z>:=MaximalOrder(K);
     SS:=MatrixAlgebra(O,2);
     PL,r:=ProjectiveLine(quo<O|P>);
     level:=#PL;
     coset_reps:=[**];
     permuted_reps:=[**];

     proj_action:=function(M,i)
         t,im :=r(PL[i]*(SS!M),true,false);
         return Index(PL,im);
     end function;


    coset_reps:=CosetMats(K,P);

    for i in [1..#PL] do
       t_iM:=coset_reps[proj_action(M,i)];
       Append(~permuted_reps,t_iM);
    end for;

     return [*coset_reps, permuted_reps *];
end function;





CoresMatAlt:=function(K,P,M,space)
     O<z>:=MaximalOrder(K);
     SS:=MatrixAlgebra(O,2);
     PL:=space`PLlevel;
     r:=space`rpl;
     level:=#PL;
     coset_reps:=[**];
     permuted_reps:=[**];

     proj_action:=function(M,i)
         t,im :=r(PL[i]*(SS!M),true,false);
         return Index(PL,im);
     end function;


    coset_reps:=CosetMats(K,P);

    for i in [1..#PL] do
       t_iM:=coset_reps[proj_action(M,i)];
       Append(~permuted_reps,t_iM);
    end for;

     return [*coset_reps, permuted_reps *];
end function;



//
// Returns the index of (0,1) in Proj Line
//
PLIndex:=function(K,level)
  O<z>:=MaximalOrder(K);
  PL,r:=ProjectiveLine(quo<O|level*O>);

  pos:=-1;

  for i in [1..#PL] do
    if PL[i][1] eq 0 and PL[i][2] eq 1 then
      pos:=i;
      break i;
    end if;
  end for;

  return pos;
end function;
