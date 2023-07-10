% qr test

N = 500; 
n_0 = zeros(1,N); 
n_1 = zeros(1,N); 

for i = 2:N
    a = random('Exponential', 5, i,i); 
    tic
    [q, r] = qr_factorization(a); 
    n_0(i) = toc; 
end

for i = 2:N
    a = random('Exponential', 5, i,i); 
    tic
    [q, r] = qr_Householder(a); 
    n_1(i) = toc; 
end

scatter(linspace(1,N,N),n_0,'k') % cresce in modo quadratico
hold on
scatter(linspace(1,N,N), n_1, 'green') % cresce peggio perché più costoso 
% computazionalmente ma anche più stabile
    