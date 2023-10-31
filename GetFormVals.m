
GetFormVals:=function(space,HECKEMATS,HECKEPRIMES,g)


  //gathers the eigenvectors

  der:=space`der;
  inn:=space`inn;

  mforms,pi:=quo<der|inn>;
  g:=Inverse(pi);


  EVs,list:=GET_EV(HECKEMATS);

  //gathers the generators of the hecke operators
  gens:=[];
  for i in HECKEPRIMES do
    t,gen:=IsPrincipal(i);
    Append(~gens,gen);
  end for;

  EV_systems:=[*<EVs[i,2],EVs[i,3]> : i in [1..#EVs]*];


  form_vals:=[];
  for i in [1..#EVs] do
    e:=EVs[i];
    if IsRationalSystem(e[3]) then
      if i gt 1 then
        offset:=&+[EVs[j,2] : j in [1..i-1]];
      else
        offset:=0;
      end if;
      for j in [1..e[2]] do
        Append(~form_vals,[*g(list[offset+j]),e[3],gens*]);
      end for;
    end if;
  end for;

  return EV_systems,form_vals;
end function;
