function [W] = pentadiag_inverseodd(a,b,c,d,e,k)
%% The fucntion compute the odd dimension cass of the pentadiagonal matrices 
%%inverse problem.
M = pentadiag_inversetest(a,b,c,d,e,k-1)
%%compute k = a-v^{T}*M*u
v = zeros(k-1,1);
v(k-2) = e;
v(k-1) = c;
u = zeros(k-1,1);
u(k-2) = d;
u(k-1) = b;
A = zeros(1,k-1);
for j = 1:k-1
    A(j)= M(k-2,j)*e + M(k-1,j)*c;
end
K = a - A*u;
%% compute the column(k)
AU = transpose(zeros(1,k-1));
for i = 1:k-1
    AU(i) = -(M(i,k-2)*d + M(i,k-1)*b)/K;
    W(i,k) = AU(i);
end

%% compute the row(k)
VA = zeros(1,k-1);
for i = 1:k-1
    VA(i) = -(M(k-2,i)*e + M(k-1,i)*c)/K;
    W(k,i) = VA(i);
end

%% compute the W(k-1,k-1)
W(k,k) = 1/K;
%% compute the (k-1,k-1)part

AUF = zeros(k-1,1);
for i = 1:k-1
   AUF(i) = M(i,k-2) * d + M(i,k-1) * b;
end
VTA = zeros(1,k-1)
for i = 1:k-1
    VTA(i) = e * M(k-2,i) + c * M(k-1,i);
end

for i = 1:k-1
    for j = 1:k-1
        AUVA(i,j) = AUF(i) * VTA(j);
    end
end

VA = zeros(1,k-1);
for i = 1:k-1
    VA(i) = e*M(k-2,i) + c*M(k-1,i);
end
VAU = VA * u;
AUVA = AUVA/(a-VAU);
for i = 1:k-1
    for j =1:k-1
        W(i,j) = M(i,j) + AUVA(i,j) ;
    end
end
end