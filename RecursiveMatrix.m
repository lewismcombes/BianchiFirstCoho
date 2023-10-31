

 RecursiveMatrix_Sp:=function(MatrixList,A)


   /* the first entry of the MatrixList is supposed to be the size of the full
      projective line  */

   M:=TensorProduct(MatrixList[2],A);
   for j in [3..#MatrixList] do
     M:=TensorProduct(MatrixList[j],M);
   end for;

   return SparseMatrix(M);
 end function;


 // non-sparse version. somehow sparse matrices cannot made to act on spaces ???

 RecursiveMatrix:=function(MatrixList,A)


   /* the first entry of the MatrixList is supposed to be the size of the full
      projective line  */

   M:=TensorProduct(MatrixList[2],A);
   for j in [3..#MatrixList] do
     M:=TensorProduct(MatrixList[j],M);
   end for;

   return Matrix(M);
 end function;




 
