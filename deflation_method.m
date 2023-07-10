% facoltativo 5

% deflation method for a generic matrix A

function eigenvalues = deflation_method(A)
        
% Deflation Method: gives an approximation of the matrix' eigenvalues
%                    
%                   INPUTS: A: square matrix, 
%
%                   OUTPUTS: b: estimated eigenvalues. 

        [n, m] = size(A); 

        if n ~= m
            error('not a square matrix')
        end

        % now for the real algorithm
        
        W = A; 
        eigenvalues = zeros(n,1); 
        for i=1:n-1
            
            m = length(W); 

            sigma = sign(W(1,1))*norm(W(:,1),2);
            utilde = W(:,1) + sigma*eye(m,1); 
            utilde = utilde/norm(utilde); 

            P = eye(m,m) - 2*(utilde*utilde');

            B = P*W%*P' 
            W = B(2:end, 2:end);  
            eigenvalues(i) = B(1,1); 

        end
