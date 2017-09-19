function [M] = pentadiag_inverse(a,b,c,d,e,k)
%
%  This function computes the inverse of the pentadiagonal matrix
%  M=(d,b,a,c,e) when the dimension k of the M is even
%  It returns the element of the T(i.e. M^{-1})  
%  It uses Matrix Mobius Transformations. 
A = [a b;c a];
C = [d 0;b d];
B = [e c;0 e];
m = k/2;
X = {};
Y = {};
T = {{}};
X{m} = zeros(2,2);
Y{1} = zeros(2,2);
M = zeros(k,k);
ZM = zeros(2,2);
M1 = 2;
%% Construct MMT for X_k 
Cinv= inv(C); 
Binv = inv(B);
MMX = [ZM  C; -Binv Binv*A]; 

% Compute the eigendecompsition of MMX
% MMX = EVX*evalX*inv(EV); where 'evalX' is a diagonal matix
[EVX, evalX]=eig(MMX); 

%% Construct MMT for Y_k
MMY = [ZM B; -Cinv Cinv*A ];

% Compute the eigendecomposition of MMY
[EVY, evalY]=eig(MMY);

%% Compute the inverse block of X_K
for i =1:m-1
    InvEVX = inv(EVX); 
    FY1 = InvEVX(1:M1,M1+1:end)*inv( InvEVX(M1+1:end,M1+1:end) );
    FY2 = ( evalX(1:M1,1:M1).^(m-i) ) * FY1; 
    Lambda2 = 1./ diag(evalX(M1+1:end,M1+1:end)); 
    FY2 = FY2 *( diag( Lambda2.^(m-i) ) ); 

    FY3 = EVX(1:M1,1:M1)*FY2+EVX(1:M1,M1+1:end); 
    FY4 = EVX(M1+1:end,1:M1)*FY2+EVX(M1+1:end,M1+1:end); 
    X{i} = FY3*inv(FY4);
end

%%compute the Y_k
for i = 2:m
    InvEVY = inv(EVY); 
    Y1 = InvEVY(1:M1,M1+1:end)*inv( InvEVY(M1+1:end,M1+1:end) );
    Y2 = ( evalY(1:M1,1:M1).^(i-1) ) * Y1; 
    Lambda1 = 1./ diag(evalY(M1+1:end,M1+1:end)); 
    Y2 = Y2 *( diag( Lambda1.^(i-1) ) ); 

    Y3 = EVY(1:M1,1:M1)*Y2+EVY(1:M1,M1+1:end); 
    Y4 = EVY(M1+1:end,1:M1)*Y2+EVY(M1+1:end,M1+1:end); 
    Y{i} = Y3*inv(Y4);
end
%%compute the T[1,1],T[m,m]
T{1,1}= inv(A-X{1});
T{m,m}= inv(A-Y{m});

%%compute the T[i,i]
for i = 2:m-1
    T{i,i} = inv(A-X{i}-Y{i});
end
%% compute the T(i,j) (i>j)
for j = 1:m-1
    for i = j+1:m
        T{i,j} = (X{i}-A)\B * T{i-1,j};
    end
end

%%compute the T(i,j)(i<j)
for i = m-1:-1:1
    for j = i+1:m
        T{i,j} = (Y{i}-A)\C * T{i+1,j};
    end
end
%% compute the (s,t) in the T
for  i = 1 : m
        M(2*i-1,2*i-1) = T{i,i}(1,1);
        M(2*i-1,2*i) = T{i,i}(1,2);
        M(2*i,2*i-1) = T{i,i}(2,1);
        M(2*i,2*i)= T{i,i}(2,2);
end

for i = 1 : m-1
    for j = i+1 :m
        M(2*i-1,2*j-1) = T{i,j}(1,1);
        M(2*i-1,2*j) = T{i,j}(1,2);
        M(2*i,2*j-1) = T{i,j}(2,1);
        M(2*i,2*j) = T{i,j}(2,2);    
    end
end
 
for i = 2 : m
    for j = 1:i-1
        M(2*i-1,2*j-1) = T{i,j}(1,1);
        M(2*i-1,2*j) = T{i,j}(1,2);
        M(2*i,2*j-1) = T{i,j}(2,1); 
        M(2*i,2*j) = T{i,j}(2,2);
    end
end


end

