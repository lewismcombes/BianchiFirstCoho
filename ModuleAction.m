

 ModuleMat:=function(char,weight,M)

   R:=CoefficientRing(Parent(M));
   K:=NumberField(R);

   if char eq 0*R then
     F:=K;
     down:=IdentityHomomorphism(R); //compatibility
   else
     F<t>,down:=ResidueClassField(char);
   end if;

   P<x,y>:=PolynomialRing(F,2);
   p:=Characteristic(F);

   k:=weight[1];
   l:=weight[2];
   d1:=weight[3];
   d2:=weight[4];

  Symm:=function(k,d,T)
     ST:=ZeroMatrix(F,k+1,k+1);
     for i in [0..k] do
       Q:=(down(T[1,1])*x+down(T[1,2])*y)^(k-i)*(down(T[2,1])*x+down(T[2,2])*y)^(i);
       for j in [0..k] do
         ST[i+1,j+1]:=MonomialCoefficient(Q,x^(k-j)*y^(j));
       end for;
     end for;
     return down(Determinant(T))^d*ST;
   end function;

   Mc:=Matrix(R,2,2,[R!Conjugate(M[1][1]),R!Conjugate(M[1][2]),R!Conjugate(M[2][1]),R!Conjugate(M[2][2])]);

   TM:=TensorProduct(Symm(k,d1,M),Symm(l,d2,Mc));

   return TM;

end function;
