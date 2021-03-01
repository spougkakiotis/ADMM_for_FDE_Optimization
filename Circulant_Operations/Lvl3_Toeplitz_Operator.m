function x = Lvl3_Toeplitz_Operator(circ_eigenvals,x,n_1,n_2,n_3)
    w = zeros(8*n_1*n_2*n_3,1);
    k = 1;
    for i = 1:2:(2*n_3-1)
        for j = 1:2:(2*n_2-1)
            w((i-1)*(2*n_2*n_1)+(j-1)*n_1+1:(i-1)*(2*n_2*n_1)+j*n_1,1) = x((k-1)*n_1+1:k*n_1,1);
            k = k+1;
        end
    end
    w = Lvl3_Circulant_Operator(circ_eigenvals,w,2*n_1,2*n_2,2*n_3,"mul");
    k = 1;
    for i = 1:2:(2*n_3-1)
        for j = 1:2:(2*n_2-1)
            x((k-1)*n_1+1:k*n_1,1) = w((i-1)*(2*n_2*n_1)+(j-1)*n_1+1:(i-1)*(2*n_2*n_1)+j*n_1,1);
            k = k+1;
        end
    end
end

