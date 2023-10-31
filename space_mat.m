
space_mat_pgl:=function(d,TA,TAi,TB,TU,TUi,TJ)

  ID:=TA^0;
  N:=ID-ID;

  TAB2:=TA*TB*TA*TB;
  TAB:=TA*TB;


  if d eq 1 then

    TU2:=TU^2;
    TJ2:=TJ^2;
    TJi:=TJ^(-1);

    E1:=VerticalJoin([N,TB+ID,N,N]);
    E2:=VerticalJoin([TB*(TAB2+TAB+ID),TAB2+TAB+ID,N,N]);
    E3:=VerticalJoin([TAi*(ID-TUi),N,(TAi-ID)*TUi,N]);
    E4:=VerticalJoin([N, TU * TB * TUi * TB * TU * TB * TUi * TB * TU * TB * TUi + TUi * TB
    * TU * TB * TUi * TB * TU * TB * TUi + TU * TB * TUi * TB * TU * TB *
    TUi + TUi * TB * TU * TB * TUi + TU * TB * TUi + TUi, TB * TUi * TB
    * TU * TB * TUi * TB * TU * TB * TUi - TUi * TB * TU * TB * TUi * TB *
    TU * TB * TUi + TB * TUi * TB * TU * TB * TUi - TUi * TB * TU * TB *
    TUi + TB * TUi - TUi, N ]);
    E5:=VerticalJoin([N, TU2 * TB * TUi * TB * TU2 * TB * TUi + TUi * TB * TU2 * TB * TUi
    + TU2 * TB * TUi + TUi, TU * TB * TUi * TB * TU2 * TB * TUi + TB *
    TUi * TB * TU2 * TB * TUi - TUi * TB * TU2 * TB * TUi + TU * TB *
    TUi + TB * TUi - TUi, N]);
    E6:=VerticalJoin([ TU * TB * TA * TUi * TB * TA * TU * TB * TA * TUi * TB + TUi * TB * TA *
    TU * TB * TA * TUi * TB + TU * TB * TA * TUi * TB + TUi * TB, TA * TUi *
    TB * TA * TU * TB * TA * TUi * TB + TA * TU * TB * TA * TUi * TB + TA *
    TUi * TB + ID, TB * TA * TUi * TB * TA * TU * TB * TA * TUi * TB -
    TUi * TB * TA * TU * TB * TA * TUi * TB + TB * TA * TUi * TB - TUi * TB,
    N]);

    J0:=VerticalJoin([N, N, N, TJ*TJ2 + TJ2 + TJ + ID]);
    J1:=VerticalJoin([ -TAi * TJi * TU, N, ID, TAi * TJi * TU - TJi * TU  ]);
    J2:=VerticalJoin([ ID, N, TJi * TA, TU * TJi * TA - TJi * TA  ]);
    J3:=VerticalJoin([ N, TJ * TB * TJ * TB * TJ * TB * TJ + TJ * TB * TJ * TB * TJ + TJ * TB * TJ +
    TJ, N, (TB * TJ)^3 + (TB * TJ)^2 + TB * TJ + ID ]);

    E:=HorizontalJoin([E1,E2,E3,E4,E5,E6,J0,J1,J2,J3]);

  elif d eq 2 then

    E1:=VerticalJoin([N,TB+ID,N,N]);
    E2:=VerticalJoin([TB*(TAB2+TAB+ID),TAB2+TAB+ID,N,N]);
    E3:=VerticalJoin([TAi*(ID-TUi),N,(TAi-ID)*TUi,N]);
    E4:=VerticalJoin([ N, TUi * TB * TU * TB * TUi * TB * TU + TU * TB * TUi * TB * TU + TUi* TB * TU
    + TU, -TUi * TB * TU * TB * TUi * TB * TU + TB * TUi * TB * TU - TUi* TB * TU +
    ID, N ]);

    J0:=VerticalJoin([N,N,N,ID+TJ]);
    J1:=VerticalJoin([ TJ * TA + ID, N, N, TA * TJ * TA + TA ]);
    J2:=VerticalJoin([ N, TJ * TB + ID, N, TB * TJ * TB + TB ]);
    J3:=VerticalJoin([ N, N, TJ * TU + ID, TU * TJ * TU + TU ]);

    E:=HorizontalJoin([E1,E2,E3,E4,J0,J1,J2,J3]);

  elif d eq 3 then

    TJi:=TJ^(-1);

    E1:=VerticalJoin([N,TB+ID,N,N]);
    E2:=VerticalJoin([TB*(TAB2+TAB+ID),TAB2+TAB+ID,N,N]);
    E3:=VerticalJoin([TAi*(ID-TUi),N,(TAi-ID)*TUi,N]);
    E4:=VerticalJoin([  TUi * TA * TB * TU * TB * TUi * TA * TUi * TA * TB + TB * TU * TB *
    TUi * TA * TUi * TA * TB + TUi * TA * TB + TB, TUi * TA * TUi * TA *
    TB * TU * TB * TUi * TA * TUi * TA * TB + TU * TB * TUi * TA * TUi * TA
    * TB + TUi * TA * TUi * TA * TB + ID, TB * TUi * TA * TUi * TA *
    TB * TU * TB * TUi * TA * TUi * TA * TB - TUi * TA * TUi * TA * TB * TU
    * TB * TUi * TA * TUi * TA * TB - TUi * TA * TB * TU * TB * TUi * TA *
    TUi * TA * TB + TB * TUi * TA * TUi * TA * TB - TUi * TA * TUi * TA *
    TB - TUi * TA * TB, N  ]);
    E5:=VerticalJoin([TB * TU * TB * TUi * TA * TB * TU * TB * TUi * TA * TB + TB * TU * TB *
    TUi * TA * TB + TB, TUi * TA * TB * TU * TB * TUi * TA * TB * TU * TB *
    TUi * TA * TB + (TU * TB * TUi * TA * TB)^2 + TUi * TA * TB * TU * TB *
    TUi * TA * TB + TU * TB * TUi * TA * TB + TUi * TA * TB + ID, TB *
    TUi * TA * TB * TU * TB * TUi * TA * TB * TU * TB * TUi * TA * TB - TUi
    * TA * TB * TU * TB * TUi * TA * TB * TU * TB * TUi * TA * TB + TB * TUi *
    TA * TB * TU * TB * TUi * TA * TB - TUi * TA * TB * TU * TB * TUi * TA *
    TB + TB * TUi * TA * TB - TUi * TA * TB, N ]);
    E6:=VerticalJoin([ TU * TB * TUi * TA * TB * TAi * TU * TB * TAi * TU * TB * TUi * TA *
    TB + TB * TAi * TU * TB * TAi * TU * TB * TUi * TA * TB - TAi * TU * TB
    * TAi * TU * TB * TUi * TA * TB - TAi * TU * TB * TUi * TA * TB + TB,
    TUi * TA * TB * TAi * TU * TB * TAi * TU * TB * TUi * TA * TB + TAi *
    TU * TB * TAi * TU * TB * TUi * TA * TB + TAi * TU * TB * TUi * TA * TB
    + TUi * TA * TB + ID, TB * TUi * TA * TB * TAi * TU * TB * TAi *
    TU * TB * TUi * TA * TB - TUi * TA * TB * TAi * TU * TB * TAi * TU * TB
    * TUi * TA * TB + TB * TAi * TU * TB * TUi * TA * TB + TB * TUi * TA *
    TB - TUi * TA * TB, N ]);

    J0:=VerticalJoin([N, N, N, TJ^5 + TJ^4 + TJ^3 + TJ^2 + TJ + ID]);
    J1:=VerticalJoin([  TJi * TUi, N, -TUi, TA * TJi * TUi - TJi * TUi  ]);
    J2:=VerticalJoin([ ID, N, TJi * TUi * TA - TUi * TA, TU * TJi * TUi * TA - TJi
    * TUi * TA  ]);
    J3:=VerticalJoin([ TB * TU * TB * TUi * TB * TAi * TB - TAi * TB, TJi * TB * TA * TB * TU
    * TB * TUi * TB * TAi * TB + TA * TB * TU * TB * TUi * TB * TAi * TB +
    TU * TB * TUi * TB * TAi * TB + TUi * TB * TAi * TB + TAi * TB +
    ID, TB * TUi * TB * TAi * TB - TUi * TB * TAi * TB, TB * TJi *
    TB * TA * TB * TU * TB * TUi * TB * TAi * TB - TJi * TB * TA * TB * TU *
    TB * TUi * TB * TAi * TB  ]);

    E:=HorizontalJoin([E1,E2,E3,E4,E5,E6,J0,J1,J2,J3]);

  elif d eq 7 then

    E1:=VerticalJoin([N,TB+ID,N,N]);
    E2:=VerticalJoin([TB*(TAB2+TAB+ID),(TAB2+TAB+ID),N,N]);
    E3:=VerticalJoin([(TAi*(ID-TUi)),N,(TAi-ID)*TUi,N]);
    E4:=VerticalJoin([ TUi * TB * TU * TB * TA * TUi * TB * TU + TUi*TB*TU, TA * TUi * TB * TU * TB
    * TA * TUi * TB * TU + TU * TB * TA * TUi * TB * TU + TA * TUi * TB * TU +
    TU, -TUi * TB * TU * TB * TA * TUi * TB * TU + TB * TA * TUi * TB * TU -
    TUi*TB*TU + ID,N ]);

    J0:=VerticalJoin([N,N,N,ID+TJ]);
    J1:=VerticalJoin([ TJ * TA + ID, N, N, TA * TJ * TA + TA ]);
    J2:=VerticalJoin([ N, TJ * TB + ID, N, TB * TJ * TB + TB ]);
    J3:=VerticalJoin([ N, N, TJ * TU + ID, TU * TJ * TU + TU ]);

    E:=HorizontalJoin([E1,E2,E3,E4,J0,J1,J2,J3]);

  elif d eq 11 then

    E1:=VerticalJoin([N,TB+ID,N,N]);
    E2:=VerticalJoin([TB*(TAB2+TAB+ID),TAB2+TAB+ID,N,N]);
    E3:=VerticalJoin([TAi*(ID-TUi),N,(TAi-ID)*TUi,N]);
    E4:=VerticalJoin([ TUi * TB * TU * TB * TA * TUi * TB * TU * TB * TA * TUi * TB * TU +
    TUi * TB * TU * TB * TA * TUi * TB * TU + TUi*TB*TU, TA * TUi * TB * TU * TB *
    TA * TUi * TB * TU * TB * TA * TUi * TB * TU + TU * TB * TA * TUi * TB *
    TU * TB * TA * TUi * TB * TU + TA * TUi * TB * TU * TB * TA * TUi * TB *
    TU + TU * TB * TA * TUi * TB * TU + TA * TUi * TB * TU + TU, -TUi * TB *
    TU * TB * TA * TUi * TB * TU * TB * TA * TUi * TB * TU + (TB * TA * TUi *
    TB * TU)^2 - TUi * TB * TU * TB * TA * TUi * TB * TU + TB * TA * TUi * TB
    * TU - TUi*TB*TU + ID, N ]);

    J0:=VerticalJoin([N,N,N,ID+TJ]);
    J1:=VerticalJoin([ TJ * TA + ID, N, N, TA * TJ * TA + TA ]);
    J2:=VerticalJoin([ N, TJ * TB + ID, N, TB * TJ * TB + TB ]);
    J3:=VerticalJoin([ N, N, TJ * TU + ID, TU * TJ * TU + TU ]);

    E:=HorizontalJoin([E1,E2,E3,E4,J0,J1,J2,J3]);

  else
    return "you stop that", "space_mat_pgl";
  end if;

  return E;

end function;





space_mat_psl:=function(d,TA,TAi,TB,TU,TUi)

  ID:=TA^0;
  N:=ID-ID;

  TAB2:=TA*TB*TA*TB;
  TAB:=TA*TB;


  if d eq 1 then

    TU2:=TU^2;

    E1:=VerticalJoin([N,TB+ID,N]);
    E2:=VerticalJoin([TB*(TAB2+TAB+ID),TAB2+TAB+ID,N]);
    E3:=VerticalJoin([TAi*(ID-TUi),N,(TAi-ID)*TUi]);
    E4:=VerticalJoin([N, TU * TB * TUi * TB * TU * TB * TUi * TB * TU * TB * TUi + TUi * TB
    * TU * TB * TUi * TB * TU * TB * TUi + TU * TB * TUi * TB * TU * TB *
    TUi + TUi * TB * TU * TB * TUi + TU * TB * TUi + TUi, TB * TUi * TB
    * TU * TB * TUi * TB * TU * TB * TUi - TUi * TB * TU * TB * TUi * TB *
    TU * TB * TUi + TB * TUi * TB * TU * TB * TUi - TUi * TB * TU * TB *
    TUi + TB * TUi - TUi ]);
    E5:=VerticalJoin([N, TU2 * TB * TUi * TB * TU2 * TB * TUi + TUi * TB * TU2 * TB * TUi
    + TU2 * TB * TUi + TUi, TU * TB * TUi * TB * TU2 * TB * TUi + TB *
    TUi * TB * TU2 * TB * TUi - TUi * TB * TU2 * TB * TUi + TU * TB *
    TUi + TB * TUi - TUi ]);
    E6:=VerticalJoin([ TU * TB * TA * TUi * TB * TA * TU * TB * TA * TUi * TB + TUi * TB * TA *
    TU * TB * TA * TUi * TB + TU * TB * TA * TUi * TB + TUi * TB, TA * TUi *
    TB * TA * TU * TB * TA * TUi * TB + TA * TU * TB * TA * TUi * TB + TA *
    TUi * TB + ID, TB * TA * TUi * TB * TA * TU * TB * TA * TUi * TB -
    TUi * TB * TA * TU * TB * TA * TUi * TB + TB * TA * TUi * TB - TUi * TB ]);

    E:=HorizontalJoin([E1,E2,E3,E4,E5,E6]);

  elif d eq 2 then

    E1:=VerticalJoin([N,TB+ID,N]);
    E2:=VerticalJoin([TB*(TAB2+TAB+ID),TAB2+TAB+ID,N]);
    E3:=VerticalJoin([TAi*(ID-TUi),N,(TAi-ID)*TUi]);
    E4:=VerticalJoin([ N, TUi * TB * TU * TB * TUi * TB * TU + TU * TB * TUi * TB * TU + TUi* TB * TU
    + TU, -TUi * TB * TU * TB * TUi * TB * TU + TB * TUi * TB * TU - TUi* TB * TU +
    ID ]);

    E:=HorizontalJoin([E1,E2,E3,E4]);

  elif d eq 3 then

    E1:=VerticalJoin([N,TB+ID,N]);
    E2:=VerticalJoin([TB*(TAB2+TAB+ID),TAB2+TAB+ID,N]);
    E3:=VerticalJoin([TAi*(ID-TUi),N,(TAi-ID)*TUi]);
    E4:=VerticalJoin([  TUi * TA * TB * TU * TB * TUi * TA * TUi * TA * TB + TB * TU * TB *
    TUi * TA * TUi * TA * TB + TUi * TA * TB + TB, TUi * TA * TUi * TA *
    TB * TU * TB * TUi * TA * TUi * TA * TB + TU * TB * TUi * TA * TUi * TA
    * TB + TUi * TA * TUi * TA * TB + ID, TB * TUi * TA * TUi * TA *
    TB * TU * TB * TUi * TA * TUi * TA * TB - TUi * TA * TUi * TA * TB * TU
    * TB * TUi * TA * TUi * TA * TB - TUi * TA * TB * TU * TB * TUi * TA *
    TUi * TA * TB + TB * TUi * TA * TUi * TA * TB - TUi * TA * TUi * TA *
    TB - TUi * TA * TB ]);
    E5:=VerticalJoin([TB * TU * TB * TUi * TA * TB * TU * TB * TUi * TA * TB + TB * TU * TB *
    TUi * TA * TB + TB, TUi * TA * TB * TU * TB * TUi * TA * TB * TU * TB *
    TUi * TA * TB + (TU * TB * TUi * TA * TB)^2 + TUi * TA * TB * TU * TB *
    TUi * TA * TB + TU * TB * TUi * TA * TB + TUi * TA * TB + ID, TB *
    TUi * TA * TB * TU * TB * TUi * TA * TB * TU * TB * TUi * TA * TB - TUi
    * TA * TB * TU * TB * TUi * TA * TB * TU * TB * TUi * TA * TB + TB * TUi *
    TA * TB * TU * TB * TUi * TA * TB - TUi * TA * TB * TU * TB * TUi * TA *
    TB + TB * TUi * TA * TB - TUi * TA * TB ]);
    E6:=VerticalJoin([ TU * TB * TUi * TA * TB * TAi * TU * TB * TAi * TU * TB * TUi * TA *
    TB + TB * TAi * TU * TB * TAi * TU * TB * TUi * TA * TB - TAi * TU * TB
    * TAi * TU * TB * TUi * TA * TB - TAi * TU * TB * TUi * TA * TB + TB,
    TUi * TA * TB * TAi * TU * TB * TAi * TU * TB * TUi * TA * TB + TAi *
    TU * TB * TAi * TU * TB * TUi * TA * TB + TAi * TU * TB * TUi * TA * TB
    + TUi * TA * TB + ID, TB * TUi * TA * TB * TAi * TU * TB * TAi *
    TU * TB * TUi * TA * TB - TUi * TA * TB * TAi * TU * TB * TAi * TU * TB
    * TUi * TA * TB + TB * TAi * TU * TB * TUi * TA * TB + TB * TUi * TA *
    TB - TUi * TA * TB ]);

    E:=HorizontalJoin([E1,E2,E3,E4,E5,E6]);

  elif d eq 7 then

    E1:=VerticalJoin([N,TB+ID,N]);
    E2:=VerticalJoin([TB*(TAB2+TAB+ID),(TAB2+TAB+ID),N]);
    E3:=VerticalJoin([(TAi*(ID-TUi)),N,(TAi-ID)*TUi]);
    E4:=VerticalJoin([ TUi * TB * TU * TB * TA * TUi * TB * TU + TUi*TB*TU, TA * TUi * TB * TU * TB
    * TA * TUi * TB * TU + TU * TB * TA * TUi * TB * TU + TA * TUi * TB * TU +
    TU, -TUi * TB * TU * TB * TA * TUi * TB * TU + TB * TA * TUi * TB * TU -
    TUi*TB*TU + ID ]);

    E:=HorizontalJoin([E1,E2,E3,E4]);

  elif d eq 11 then

    E1:=VerticalJoin([N,TB+ID,N,N]);
    E2:=VerticalJoin([TB*(TAB2+TAB+ID),TAB2+TAB+ID,N]);
    E3:=VerticalJoin([TAi*(ID-TUi),N,(TAi-ID)*TUi]);
    E4:=VerticalJoin([ TUi * TB * TU * TB * TA * TUi * TB * TU * TB * TA * TUi * TB * TU +
    TUi * TB * TU * TB * TA * TUi * TB * TU + TUi*TB*TU, TA * TUi * TB * TU * TB *
    TA * TUi * TB * TU * TB * TA * TUi * TB * TU + TU * TB * TA * TUi * TB *
    TU * TB * TA * TUi * TB * TU + TA * TUi * TB * TU * TB * TA * TUi * TB *
    TU + TU * TB * TA * TUi * TB * TU + TA * TUi * TB * TU + TU, -TUi * TB *
    TU * TB * TA * TUi * TB * TU * TB * TA * TUi * TB * TU + (TB * TA * TUi *
    TB * TU)^2 - TUi * TB * TU * TB * TA * TUi * TB * TU + TB * TA * TUi * TB
    * TU - TUi*TB*TU + ID ]);

    E:=HorizontalJoin([E1,E2,E3,E4]);

    else
      return "you stop that","space_mat_psl";
    end if;

  return E;

end function;
