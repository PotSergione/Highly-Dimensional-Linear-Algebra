% Q R FACTORIZATION VIA GRAHM-SCHMIDT STABLE ALGORITHM: 

% returns the Q r factorization of a given matrix, where q is orthogonal
% and r is uooer triangular such that V = Q*r. The stability of the
% algorithm comes from the use of w as placeholder for the colum of v,
% which, if used, would result in an unaxceptable cumulation of consecutive
% errors. 

function [Q,r] = qr_factorization(matrix)
    
    [m, n] = size(matrix);  
    
    % initializing the output matrices
    Q = zeros(m,n); 
    r = zeros(n,n); 
    % here we initialize the first vector of the Q matrix
    %Q(:,1) = matrix(:,1)/norm(matrix(:,1)); 
    W = matrix; % guarantees the stability

    for j=1:n
%         w = matrix(:,j); % guarantees the stability.
%         for i=1:j-1
%             
%             r(i,j) = w'*Q(:,i); % graham schmidt method to find orthogonal
%             w = w - r(i,j)*Q(:,i); % basis
%             
%         end
%         r(j,j) = norm(w); % normalizing to obtain normality
%         Q(:,j) = w/r(j,j); 
        
        % without inner loop: 
        % we take advantage of the matrix product capabilities of matlab
        r(1:j-1,j) = Q(:,1:j-1)'*W(:,j); 
        W(:,j) = W(:,j)-Q(:,1:j-1)*r(1:j-1,j); 

        r(j,j) = norm(W(:,j),2); % int he loop free we normalize at the end
        Q(:,j) = W(:,j)/r(j,j); 

    end
    
return