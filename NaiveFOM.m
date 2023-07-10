function x = NaiveFOM(A, b, x0, tol)
% NaiveKrylov: gives the Krylov approximating solution x. 
%
%
%              INPUTS: A square matrix n by n
%                      b image vector of size n by 1
%                      x0 n by 1 initial guess to be improved by the method.
%                      tol specifies the mazimum allowed  tolerance in the 
%                          algorithm.
%
%              The method is Naive because it does not use Householder
%              reflectors, but instead it relies on a straightforward 
%              implementation of the modified Graham-Shmidt algorithm. 

    n = length(A);  
    H = zeros(n,n);
    V = zeros(n,n); 
    W = zeros(n,n); 

    r0 = (b - A*x0); 
    beta = norm(r0,2); 
    V(:,1) = r0/beta; 
    
    for j = 1:n
        W(:,j) = A*V(:,j); 

        for i = 1:j
            H(i,j) = V(:,i)'*W(:,j); 
            W(:,j) = W(:,j) - H(i,j).*V(:,i);
        end
     
        
        if norm(W(:,j),2) <= tol 
            break
        end
        
        H(j+1,j) = norm(W(:,j),2);
        V(:,j+1) = W(:,j)/H(j+1,j); 
    end
   
   y = H\(beta.*eye(n,1));
   x = x0 + V*y; 
end