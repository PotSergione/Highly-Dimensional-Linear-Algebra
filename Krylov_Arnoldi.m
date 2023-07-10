% KRYLOV SPACE METHOD. This script finds a base of the Krylov space
% generated by a matrix A and a related vector b. 
% This implementation uses Arnoldi's method, which is a clever adaptation 
% of Graham/Shmidt's algorithm. 

function [V, H] = Krylov_Arnoldi(A, v, tol)
        % Takes as inputs A and v, where A is a square matrix and v is a
        % unit norm vector
        
        n = length(A);  
        H = zeros(n,n);
        V = zeros(n,n); 
        W = zeros(n,n); 
        
        m = 1; 
        V(:,1) = v; 
        

        for j = 1:n
            W(:,j) = A*V(:,j); 

            for i = 1:j
                H(i,j) = V(:,i)'*W(:,j); 
                W(:,j) = W(:,j) - H(i,j).*V(:,i);
            end
            
            H(j+1,j) = norm(W(:,j),2); 
            
            if H(j+1,j) <= tol 
                m = j; 
                V = V(:, 1:m); 
                H = H(:, 1:m);
            return 
            end

            V(:,j+1) = W(:,j)/H(j+1,j); 
        end
        
       V = V(:, 1:m); 
       H = H(:, 1:m); 
end