function [x,y] = luinverse(l,u,k)
x = zeros(k,k);
y = zeros(k,k);
    for i = 1:k
        x(i,1) = 1/u(i,1);
        y(i,3) = 1;
    end
    
    for i = 1:k-1
        x(i,2) = -(u(i,2)*x(i,1))/u(i+1,1);
        y(i+1,2) = - l(i+1,2);
    end
    for i = 1:k-2
        x(i,3) = -(x(i,1)*u(i,3)+x(i,2)*u(i+1,2))/u(i+2,1);
        y(i+2,1) = -(y(i+2,2)*l(i+1,2)+y(i+2,3)*l(i+2,1));
    end
end