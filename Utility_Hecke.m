


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// This function gather ONLY RATIONAL eigenvalue systems from a given
// set of commuting matrices.
// needs to be loaded after Hecke functionalities and Artin Algebras


GET_EV_RAT:=function(list)
   A:=Parent(list[1]);
   D:=Decomposition(A);
   EV_List:=<>;

   for j in [1..#D] do
       TT:=[BaseChange(T,D[j]): T in list];
       m:=NumberOfRows(TT[1]);
       ev_size:=[];
	   for M in TT do
	       total:=0;
		   for item in Eigenvalues(M) do
		       total:=total+item[2];
			end for;
	       Append(~ev_size,total);
	   end for;

       if m eq Min(ev_size) then
         Zx<x>:=PolynomialRing(Integers());
         fx:=Zx!DefiningPolynomial(BaseRing(TT[1]));
         mult:=SetToSequence(Eigenvalues(TT[1]))[1][2];
         ev:=[SetToSequence(Eigenvalues(M))[1][1] : M in TT];
         Append(~EV_List,<fx,mult,ev>);
       end if;
	end for;

	 return EV_List;
end function;

//////
// This is the function to use. all others are buggy and incomplete as far as I can tell
//////

GET_EV:=function(list)

  char:=Characteristic(CoefficientRing(list[1]));
  if char eq 0 then
    F<w>:=OptimizedRepresentation(SplittingField(&*[MinimalPolynomial(m) : m in list]));
  else
    F<w>:=SplittingField(&*[MinimalPolynomial(m) : m in list]);
  end if;

  list:=[ChangeRing(h,F) : h in list];
  A:=MatrixAlgebra(list);

  if Dimension(A) eq 0 then
    EV_list:=[<DefiningPolynomial(F),Ncols(list[1]),[0 : i in [1..#list]]>];
    basis:=Basis(Kernel(list[1]));
  else
    D:=Decomposition(A);
    EV_list:=<>;

    for j in [1..#D] do
        TT:=[BaseChange(T,D[j]): T in list];
        m:=NumberOfRows(TT[1]);
        ev_size:=[];
      for M in TT do
          total:=0;
        for item in Eigenvalues(M) do
            total:=total+item[2];
       end for;
          Append(~ev_size,total);
      end for;

        if m eq Min(ev_size) then
          Zx<x>:=PolynomialRing(Integers());
          if char eq 0 then
            fx:=Zx!DefiningPolynomial(BaseRing(TT[1]));
          else
            fx:=DefiningPolynomial(BaseRing(TT[1]));
          end if;
          mult:=SetToSequence(Eigenvalues(TT[1]))[1][2];
          ev:=[SetToSequence(Eigenvalues(M))[1][1] : M in TT];
          Append(~EV_list,<fx,mult,ev>);
        end if;
    end for;
    basis:=&cat[Rows(D[i,1]) : i in [1..#D]];
  end if;

  return EV_list,basis;
end function;






//
// Takes a list of ideals and returns the relevant Hecke matrices
//
GetHeckeMatrices:=function(space,list)
  if space`dim ne 0 then // doing hecke stuff
    heckes:=[];
    for TP in list do
        time Append(~heckes,HECKE(TP,space));
    end for;
    return heckes;
  else
    return [];
  end if;
end function;




  IsRationalSystem:=function(evs)
    isitrational:=true;
    for e in evs do
      if not e in Integers() then
        isitrational:=false;
      end if;
    end for;
    return isitrational;
  end function;






















//
