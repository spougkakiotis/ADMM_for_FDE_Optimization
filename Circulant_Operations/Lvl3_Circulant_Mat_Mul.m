function C1 = Lvl3_Circulant_Mat_Mul(C1,C2)
% ============================================================================================================================================================ %
% C1 = Lvl3_Circulant_Mat_Mul(C1,C2)
% ------------------------------------------------------------------------------------------------------------------------------------------------------------ %
% This function takes as an input Level-3 Circulant matrix, represented by their eigenvalues, 
% and performs the matrix multiplication, by performing a component-wise product of their eigenvalues.
% ============================================================================================================================================================ %
    C1 = C1 .* C2;
end

