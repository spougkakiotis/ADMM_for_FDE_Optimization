function cmat = Multilevel_Circulant_Extrapolation(T, lvls, f_lvl)
% ================================================================================ %
% [cmat_eigs] = Multilevel_Circulant_Extrapolation(T, lvls, f_lvl)
% -------------------------------------------------------------------------------- %
% This function takes as an input the Toeplitz matrix of our FDE, a matrix 
% containing the dimension size of each level and the final level (i.e.
% the number of levels). Then, it recursively computes the 
% Level-k Circulant extrapolation of the given Level-k Toeplitz matrix.
% It returns the eigenvalues of the circulant matrix.
% ________________________________________________________________________________ %
    
    current_lvl = size(lvls,1);
    if (current_lvl == 1)
        cmat = Lvl1_Circulant_Extrapolation(T);
        return
    end
    if (f_lvl ~= 3)
        errror("Current implementation does not support other lvl options.")
    end
    nc = lvls(end,1);
    n = 2*nc;
    for i = 1:current_lvl-1
        k = current_lvl - i;
        n = n*(2*lvls(k,1));
    end
    nrC = n/(2*nc);
    nrT = nrC/(2^(current_lvl-1));
    cmat = zeros(n,1);
    k = current_lvl - 1;
    for i = 1:(nc)
        cmat((i-1)*nrC + 1: i*nrC,1) = Multilevel_Circulant_Extrapolation(T((i-1)*nrT+1:i*nrT,1),lvls(1:k,1),f_lvl);
    end
    if (current_lvl == f_lvl) % This implementation is correct only for f_lvl = 3.
        k = current_lvl - 1;
        cmat(nc*nrC+1:(nc+1)*nrC,1) = zeros(nrC,1);       
        for i = (nc+2):(2*nc)
            cmat((i-1)*nrC + 1: i*nrC,1) = Multilevel_Circulant_Extrapolation(T((3*nc-i)*nrT+1:(3*nc-i+1)*nrT,1),lvls(1:k,1),f_lvl);
        end
        cmat = reshape(cmat,nrC,2*nc);
        k = current_lvl-1;
        for i = 1:(2*nc)
            tmp_mat = cmat(:,i);
            tmp_mat = reshape(tmp_mat,(2.*lvls(1:k,1))');
            tmp_mat = fft(tmp_mat,[],1);
            tmp_mat = tmp_mat';
            tmp_mat = fft(tmp_mat,[],1);
            cmat(:,i) = reshape(tmp_mat,nrC,1);
        end
        cmat = cmat';
        cmat = fft(cmat,[],1);
        cmat = reshape(cmat,n,1);
    else
        cmat(nc*nrC+1:(nc+1)*nrC,1) = zeros(nrC,1);
        for i = (nc+2):(2*nc)
             cmat((i-1)*nrC + 1: i*nrC,1) = cmat((2*nc-i + 1)*nrC+1:(2*nc-i + 2)*nrC,1);
        end
    end
end

function circ_vals = Lvl1_Circulant_Extrapolation(T)
    circ_vals = [T;0;T(end:-1:2,1)];
end
