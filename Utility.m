

StandardMats:=function(d)

  K:=QuadraticField(-d);
  O<z>:=MaximalOrder(K);

  A:=Matrix(O,2,2,[1,1,0,1]);
  Ai:=Matrix(O,2,2,[1,-1,0,1]);
  B:=Matrix(O,2,2,[0,-1,1,0]);    /* order 2 */
  U:=Matrix(O,2,2,[1,z,0,1]);
  Ui:=Matrix(O,2,2,[1,-z,0,1]);

  if d in {2,7,11} then
    J:=Matrix(O,2,2,[-1,0,0,1]);
  elif d in {1,3} then
    J:=Matrix(O,2,2,[z,0,0,1]);
  else
    J:="stop that";
  end if;

  return A,Ai,B,U,Ui,J;
end function;


DEWORD:=function(list)
  M:=list[1,1]^0;
  Z:=Integers();
  for i in list do
    M*:=i[1]^Z!i[2];
  end for;
  return M;
end function;



//
// Write W as a word in [A,i], [B,1], [U, j], [J,k]
//
WORD:=function(W,d)

  K:=QuadraticField(-d);
  O<z>:=MaximalOrder(K);

  A,Ai,B,U,Ui,J:=StandardMats(d);

  if d eq 1 then

    w:=W;
    S:=[];

    if Determinant(W) eq z then
      Append(~S,[*J,1*]);
      W:=J^(-1)*W;
    elif Determinant(W) eq -z then
      Append(~S,[*J,-1*]);
      W:=J*W;
    elif Determinant(W) eq -1 then
      S cat:= [[*J,1*],[*J,1*]];
      W:=J^2*W;
    end if;

    if Abs(Norm(W[2,1])) ge Abs(Norm(W[1,1])) then
      W:=B*W;
      Append(~S,[*B,1*]);
    end if;

    while Norm(W[2,1]) ne 0 do
      q:=O!W[1,1] div O!W[2,1];
      seq:=Eltseq(O!q);
      Append(~S,[*A,seq[1]*]);
      Append(~S,[*U,seq[2]*]);
      Append(~S,[*B,1*]);
      Q:=Matrix(O,2,2,[1,-q,0,1]);
      W:=B*Q*W;
    end while;

    if W[1,1] eq z then
      S cat:= [[*B,1*],[*U,-1*],[*B,1*],[*U,1*],[*B,1*],[*U,-1*]];
      W:=U^(-1)*B*U*B*U^(-1)*B*W;
    elif W[1,1] eq -z then
      S cat:= [[*U,-1*],[*B,1*],[*U,1*],[*B,1*],[*U,-1*],[*B,1*]];
      W:=B*U^(-1)*B*U*B*U^(-1)*W;
    end if;

    if W[1,1] eq 1 then
      seq:=Eltseq(O!W[1,2]);
      Append(~S,[*A,seq[1]*]);
      Append(~S,[*U,seq[2]*]);
    else
      seq:=Eltseq(O!W[1,2]);
      Append(~S,[*A,-seq[1]*]);
      Append(~S,[*U,-seq[2]*]);
    end if;

    assert (DEWORD(S) eq w or DEWORD(S) eq -w or DEWORD(S) eq z*w or DEWORD(S) eq -z*w);


  elif d eq 3 then

    w:=W;

    S:=[];

    while Determinant(W) ne 1 do
      Append(~S,[*J,1*]);
      W:=J^(-1)*W;
    end while;

    if Abs(Norm(W[2,1])) ge Abs(Norm(W[1,1])) then
      W:=B*W;
      Append(~S,[*B,1*]);
    end if;

    while Norm(W[2,1]) ne 0 do
      q:=O!W[1,1] div O!W[2,1];
      seq:=Eltseq(O!q);
      Append(~S,[*A,seq[1]*]);
      Append(~S,[*U,seq[2]*]);
      Append(~S,[*B,1*]);
      Q:=Matrix(O,2,2,[1,-q,0,1]);
      W:=B*Q*W;
    end while;

    R:=Matrix(O,2,2,[z,0,0,z^5]);

    while W[1,1] ne 1 do
      S cat:=[[*B,1*],[*A,1*],[*B,1*],[*U,1*],[*B,1*],[*U,-1*],[*B,1*],[*A,-1*]];
      W:=R^(-1)*W;
    end while;


      seq:=Eltseq(O!W[1,2]);
      Append(~S,[*A,seq[1]*]);
      Append(~S,[*U,seq[2]*]);


    M:=DEWORD(S);
    unitmats:=[z^k*w : k in [1..6]];

    assert M in unitmats;

  elif d in {2,7,11} then

    w:=W;
    S:=[];

    if Determinant(W) eq -1 then
      Append(~S,[*J,1*]);
      W:=J*W;
    end if;

    if Abs(Norm(W[2,1])) ge Abs(Norm(W[1,1])) then
      W:=B*W;
      Append(~S,[*B,1*]);
    end if;

    while Norm(W[2,1]) ne 0 do
      q:=O!W[1,1] div O!W[2,1];
      seq:=Eltseq(O!q);
      Append(~S,[*A,seq[1]*]);
      Append(~S,[*U,seq[2]*]);
      Append(~S,[*B,1*]);
      Q:=Matrix(O,2,2,[1,-q,0,1]);
      W:=B*Q*W;
    end while;


    if W[1,1] eq 1 then
      seq:=Eltseq(O!W[1,2]);
      Append(~S,[*A,seq[1]*]);
      Append(~S,[*U,seq[2]*]);
    else
      seq:=Eltseq(O!W[1,2]);
      Append(~S,[*A,-seq[1]*]);
      Append(~S,[*U,-seq[2]*]);
    end if;

    assert (DEWORD(S) eq w or DEWORD(S) eq -w);

  else
    return "you stop that";
  end if;

  return S;

end function;






// ================================================================================
// IDENTIFYING IDEALS VIA THEIR MAGMA PRIME FACTORIZATIONS
// =================================================================================

/*  builds ideal from our standard code  [<<p_1,a_1>,exp_1>, ... , <<p_n,a_n>,exp_n> */

get_ideal:=function(list,K)
  final:=1*Integers(K);
  for item in list do
      I:=item[1][1]*Integers(K);
      J:=Factorization(I)[item[1][2]][1];
      final:=final*J^item[2];
  end for;
  return final;
end function;


identify_prime:=function(J)
   p:=PrimeDivisors(Norm(J))[1];
   fac:=Factorization(p*Order(J));
   if J eq fac[1][1] then
      return <p,1>;
   else
      return <p,2>;
   end if;

end function;

// factorizes and identifies factors with exponents

identify:=function(I)
    fac:=Factorization(I);
    list:=[];
    for factor in fac do
        Append(~list,<identify_prime(factor[1]),factor[2]>);
    end for;

    return <Norm(I),list>;
end function;

//======================================================================
// IDENTIFYING IDEALS VIA THEIR HNF BASES
// ====================================================================

detect_hnf:=function(J,M)
      N:=Norm(J);
      a:=M[1,1];  d:=M[1,2];
	  b:=M[2,1];  c:=M[2,2];
	  if (d eq 0) and (N eq a*c) and (b in [0..a-1]) then
	     return 1;
	  else
         return 0;
      end if;
end function;

HNF_basis:=function(J)
      N:=Norm(J);
      M:=BasisMatrix(J);
	  Mt:=Matrix(Integers(),2,2,[M[1,2],M[1,1],M[2,2],M[2,1]]);
      HN:=HermiteForm(Mt);
	  H:=Matrix(Integers(),2,2,[HN[2,2],HN[2,1],HN[1,2],HN[1,1]]);
      c:=H[2,2];  b:=H[2,1];   a:=H[1,1];
	  if c lt 0 then
	     c:=-c;   b:=-b;
	  end if;

	  b:=(b mod a);
	  assert 1 eq detect_hnf(J,Matrix(Integers(),2,2,[a,0,b,c]));
	  return [Norm(J), b,c];
end function;

get_ideal_hnf:=function(K,list)
      O<w>:=Integers(K);
	  N:=list[1];   /* norm */
      gen1:=O!(N/list[3]);
	  gen2:=O!(list[2]+w*list[3]);
	  return ideal<O| gen1, gen2>;
end function;




//======================================================================================

    /*  K gives field, I is the ideal in Gamma0(I)  */

CuspNumber_GL:=function(K,I)
      O<z>:=Integers(K);

      if Discriminant(O) eq -4 then U:={1,z,-1,-z};
      elif Discriminant(O) eq -3 then U:={z,z^2,z^3,z^4,z^5,z^6};
      else U:={1,-1};
      end if;

      sum:=0;
      for d in Divisors(I) do
		 R,h:=quo<O | d+I*d^-1>;
         sum:=sum+(#UnitGroup(R)/#h(U));
     end for;
      return ClassNumber(K)*sum;
end function;



CuspNumber_SL:=function(K,I)
      O:=Integers(K);
      sum:=0;
      for d in Divisors(I) do
		 R:=quo<O | d+I*d^-1>;
         sum:=sum+#UnitGroup(R);
     end for;
     return ClassNumber(K)*sum;
end function;




//////////////////////////////
// WAS PREVIOUSLY HECKE_AUX //
//////////////////////////////

/* =================================================================
given a sequence G of algebraic integers a+bw and index r, returns the first coordinate "a" of the r'th element of G
============================================================= */
R:=function(G,r)

h:=Eltseq(G[r])[1];
h:=Integers()!h;
return h;
end function;

/* =================================================================
returns the second coordinate "b"         Be careful when the field is not Q(-1)
============================================================== */
I:=function(G,r)

h:=Eltseq(G[r])[2];
h:=Integers()!h;
return h;
end function;


// needs to be here for strange compatibility reasons
MatPow:=function(matrix,pow)
  mat_list:=[];
  if pow eq 0 then
    return matrix^0;
  end if;
  for i in [1..Abs(pow)] do
    Append(~mat_list,matrix^Sign(pow));
  end for;
  return &*mat_list;
end function;

/*=================================================================
Given a square matrix M and m, it will return
 1+M+..+M^(m-1) or -M^m*(1+M+...+M^(-m-1)) if k is negative and 0 if k=0
=====================================*/

HeckeMatP:=function(N,M,m)

W:=N;

if m lt 0 then
for j in [0..(-m-1)] do
W:=W+MatPow(M,j);
end for;
W:= -W*MatPow(M,m);
end if;

if m gt 0 then
for j in [0..(m-1)] do
W:=W+MatPow(M,j);
end for;
end if;


return W;
end function;
