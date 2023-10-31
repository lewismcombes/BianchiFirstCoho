# BianchiFirstCoho
Magma code to compute H^1 of Bianchi modular forms of general level, weight and character. Computes the Hecke action. Can compute all of this in characteristic 0 and p. Adapted from various pieces of code originally written by Haluk Şengün.

Contains 

 + loading_script_psl.m, the main way to use this code. Change the various quantities like generator of the level, weights, determinant twists, character and characteristic. Computes H^1 using Fox calculus per e.g. Haluk Şengün's thesis _"Serre's conjecture over imaginary quadratic fields"_. Returns cocycle classes in the list form_vals, which contains cocycles corresponding to Hecke eigenvalues also contained in the list, and the generators of Hecke matrices. Only works for Euclidean imagianry quadratic fields.
 + GetFormVals.m, which packages all the eigenvalue systems up into the form_vals list
 + Hecke.m, which computes the Hecke action.
 + Hecke_poly.m, a utility function that helps compute the Hecke action by treating evaluations of the cocycle as variables in a polynomial ring, making computation much more straightforward to code.
 + LevelData_pgl.m, packages the projection line action into one record.
 + ModuleAction.m, computes the action of PSL_2(Z_K) on the weight module, allowing Bianchi forms with non-trivial weight to be computed.
 + Periods.m, computes the "period average" of a Bianchi form. Only works with trivial weight, uses an averaging trick in e.g. Section 2.8 of Cremona's book _"Algorithms for modular elliptic curves"_ to compute a special L-value (up to some scaling we do not address).
 + ProjAction.m, computes the projective line action of PSL_2(Z_K) to allow Bianchi forms with non-trivial level to be computed.
 + space_mat.m, which is where most of what makes the code tick is contained. Contains relations between matrices obtained from Fox calculus.
 + Space_pgl.m, combines the level and weight actions to produce the actual spaces of cocycles and coboundaries.
 + Utility.m, some useful functions; particularly useful is WORD, which computes the word decomposition of a matrix in terms of the generators used throughout.
 + Utility_Hecke.m, some useful functions that aid in the Hecke computation. 
