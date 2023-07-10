% QR FACTORIZATION VIA HOUSEHOLDER METHOD
% Is the most stable method, based on reflections along a prespecified
% axis, the algorithm constructs at every step the projection matrices. 

function [Q, W] = qr_Householder(matrix)

        n = length(matrix); 
        Q = eye(n); 
        W = matrix; 

        for i=1:n

            x = W(i:n, i);  % selecting the column to modify
            sigma = sign(x(1))*norm(x); 
            u = x + sigma*eye(n-i+1, 1); 
            
            P = eye(n); % building the new matrix for iteration i; 
            P(i:n, i:n) = eye(n-i+1) - u*u'/(sigma*(sigma + x(1)));

            Q = P*Q;
            W = P*W; 
        end
        Q = -Q'; 
        W = -W; 