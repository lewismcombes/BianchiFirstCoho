//
// Returns 1+M+M^2+...+M^(m-1) (or the correct version for an inverse)
//
MatP:=function(M,m)
  Z:=Integers();
  W:=M-M;
  if m lt 0 then
    for j in [0..(-m-1)] do
      W:=W+M^j;
    end for;
    W:= -W*M^Z!(m);
  end if;
  if m gt 0 then
    for j in [0..Z!(m-1)] do
      W:=W+M^j;
    end for;
  end if;
  return W;
end function;



//
// Returns a matrix in Gamma0(N) which takes 0 to k/p
//
PathMat:=function(K,N,k,p)

  O<z>:=MaximalOrder(K);
  J:=N*O;
  ord:=#UnitGroup(quo<O|J>);

  a:=p^(ord-1) mod N;

  t:=(a*p-1)/N;

  coeffs:=Bezout(p,k);

  u:=-t*coeffs[1]/GCD(p,k);
  v:=t*coeffs[2]/GCD(p,k);

  M:=Matrix(O,2,2,[a+u*N,k,v*N,p]);
  assert Determinant(M) eq 1;
  assert M[2,1] mod N eq 0;

  return M;
end function;



//
// Given f a cocycle and W a matrix, returns f(w)
//
EvaluateCocycle:=function(space,f,W)

  K:=space`field;
  O<z>:=MaximalOrder(K);

  d:=space`d;
  A,Ai,B,U,Ui,J:=StandardMats(d);

  der:=space`der;
  inn:=space`inn;

  TA:=space`ta;
  TB:=space`tb;
  TU:=space`tu;
  TJ:=space`tj;

  n:=space`coeff_dim;

  f_A:=[];
  f_B:=[];
  f_U:=[];
  f_J:=[];

  for i in [1..n] do // fills the f_A, f_B, f_U vectors with their respective parts of the full vector of f
    Append(~f_A,f[i]);
    Append(~f_B,f[n+i]);
    Append(~f_U,f[2*n+i]);
    Append(~f_J,f[3*n+i]);
  end for;

  f_A:=Vector(n,f_A);
  f_B:=Vector(n,f_B);
  f_U:=Vector(n,f_U);
  f_J:=Vector(n,f_J);

  Z:=Integers();

  decomp:=WORD(W,d);
  t_decomp:=[]; // creating the decomposition e.g. [A,3],[B,1],[U,-4],[B,1] -> TA^3,TB,TU^-4,TB
  for i in [1..#decomp] do
    if decomp[i][1] eq A then
      Append(~t_decomp,TA^(Z!decomp[i][2]));
    elif decomp[i][1] eq B then
      Append(~t_decomp,TB^(Z!decomp[i][2]));
    elif decomp[i][1] eq U then
      Append(~t_decomp,TU^(Z!decomp[i][2]));
    elif decomp[i][1] eq J then
      Append(~t_decomp,TJ^(Z!decomp[i][2]));
    else
      print "you stop that", "t_decomp";
    end if;
  end for;

  image:=f_A-f_A; //zero vector
  N:=Parent(TA)!0; // zero matrix

  for i in [1..#decomp] do //adding on each term in the cocycle decomposition
    if decomp[i][1] eq A then
      image+:=f_A*MatP(TA,decomp[i][2])*&*t_decomp[i+1..#t_decomp];
    elif decomp[i][1] eq B then
      image+:=f_B*MatP(TB,decomp[i][2])*&*t_decomp[i+1..#t_decomp];
    elif decomp[i][1] eq U then
      image+:=f_U*MatP(TU,decomp[i][2])*&*t_decomp[i+1..#t_decomp];
    elif decomp[i][1] eq J then
      image+:=f_J*MatP(TJ,decomp[i][2])*&*t_decomp[i+1..#t_decomp];
    else
      print "you stop that","decomp";
    end if;
  end for;

  return image;
end function;



//
// Calculates the period sum f(M_a) for a in Z_K/P
//
PeriodAverage:=function(space,form,r)

  d:=space`d;
  K:=QuadraticField(-d);
  O<z>:=MaximalOrder(K);
  level:=space`level;
  t,N:=IsPrincipal(level);
  p:=form[3][r];
  ModPoints:=Exclude(Elements(quo<O|p*O>),0); // 0 borks the code, and integral from 0 to 0 doesn't contribute, so we ignore it

  mats:=[];
  images:=[];
  image:=EvaluateCocycle(space,form[1],Matrix(O,2,2,[1,0,0,1]));

  for k in ModPoints do
    Append(~mats,PathMat(K,N,k,p));
    Append(~images,EvaluateCocycle(space,form[1],mats[#mats]));
    image+:=images[#images];
  end for;
  R:=Parent(form[1][1]);
  if R!(1+Norm(O!p)-form[2][r]) ne 0 then
    return form[3][r],"good", image[space`id_index]/(1+Norm(O!p)-form[2][r]), [images[i][space`id_index] : i in [1..#images]];
    //return form[3][r],"good", image[space`id_index], [images[i][space`id_index] : i in [1..#images]];
  else
    return form[3][r], " bad", image[space`id_index], [images[i][space`id_index] : i in [1..#images]];
  end if;
end function;
