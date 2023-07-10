% QR FACTORIZATION VIA ROTATION MATRICES: 
% the qr factorization via rotation matrices is useful for sparse matrices
% because it avoids the n^2 computations of the classic Graham Schmidt. 
% This is generally written for square matrices of dim n
% Note: for large, dense matrices the multiplication errors compound,
% making the computation unreliable! The method should be used only for
% SPARSE matrices.

function [Q, r] = qr_reflections(matrix)
    
    format long

    n = length(matrix);
    
    W = matrix; 
    Q = eye(n); 

    for i = 1:n
        x = W(i,i); 

        for j = i+1:n
            y = W(j,i);

            if y ~= 0 
                % this check is fundamental, otherwise the algorithm 
                % does not work!
                G = eye(n); 
    
                hold = sqrt(x^2 + y^2); 
    
                c = abs(x)/hold; 
                s = abs(y)/hold; 
    
                G(i,i) = c; 
                G(i,j) = -s; 
                G(j,i) = s; 
                G(j,j) = c; 
                 
                Q = G*Q; 
                W = G*W; 
            end
        end
    end
    
    Q = Q'; 
    r = W; 