function [D] = pentadiagonalodd(a,b,c,d,e,k)
%
%  This function computes the inverse of the pentadiagonal matrix
%  M=(d,b,a,c,e) when the dimension k of the M is even
%  It returns the element of the T(i.e. M^{-1})  
%  It uses Matrix Mobius Transformations. 
invs2 = @(a) [a(2,2),-a(1,2);-a(2,1),a(1,1)]/(a(1,1)*a(2,2)-a(1,2)*a(2,1));

A = [a b;c a];
C = [d 0;b d];
B = [e c;0 e];

X = zeros(2,k);
X(:,k-1:k) = 0;
Y = zeros(2,k);
Y(:,1:2) = 0;
 M = zeros(2,k);
 D = zeros(1,k+1);
 T = zeros(2,2);
% ZM = zeros(2,2);
% M1 = 2;
for i = k-3:-2:1
    X(:,i:i+1) = C * invs2(A-X(:,i+2:i+3)) * B;
end

for i = 3:2:k-1
    Y(:,i:i+1) = B * invs2(A-Y(:,i-2:i-1))* C;
end

M(:,1:2)= invs2(A-X(:,1:2));
M(:,k-1:k)= invs2(A-Y(:,k-1:k));
% 
% %%compute the T[i,i]
for i = 3:2:k-3
     M(:,i:i+1) = invs2(A-X(:,i:i+1)-Y(:,i:i+1));
%     T{i,i} = inv(A-X{i}-Y{i});
end

for i = 1:2:k-1
    D(1,i) = M(1,i);
    D(1,i+1) = M(2,i+1);
end

T = M(:,k-1:k);
n = a -(c*e*T(1,1) + c*d*T(2,1)+ b*e*T(1,2)+b*d*T(2,2));
MM = [(e*T(1,1)+c*T(2,1))*(d*T(1,1)+b*T(1,2)),(e*T(1,2)+c*T(2,2))*(d*T(1,1)+b*T(1,2));
    (e*T(1,2)+c*T(2,2))*(d*T(1,1)+b*T(1,2)),(e*T(1,2)+c*T(2,2))*(d*T(2,1)+b*T(2,2))];
T = T+ 1/n * MM;
D(1,k-1) = T(1,1);
D(1,k)=T(2,2);
D(1,k+1)= 1/n;
end
