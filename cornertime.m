repeat=10;
numscale=10;
timetable3=zeros(numscale,repeat);
timetable4 = zeros(numscale,repeat);
scaletable = linspace(1100,11000,numscale);
for dimension = 1:numscale
    for t = 1:repeat
        dimension
        t
        %% This compute the pentadiagonal inverse problem by using the algorithm of solving the Toplizi matrices.
        %%first you need to enter the element of the pentadiagonal matrices.
       
         a = 1;
         b = 2;
         c = 2;
         d = 3;
         e = 3;
        %%ENTER THE dimension of the matrix to choose function
          k = scaletable(dimension);
          tic;
          W = corner(a,b,c,d,e,k);
             % W;
         timetable3(dimension,t)=toc;
         S = zeros(k,k);
        for  i = 1:k
            for j = 1:k
                if j == i
                    S(i,j) = a;
                end
                if j == i + 1
                    S(i,j) = b;
                end
                if j == i+2
                    S(i,j) = d;
                end
                if j == i-1
                    S(i,j)= c;
                end
                if j == i-2
                    S(i,j) = e;
                end
            end
        end
       tic;
       [l,u] = pentaLU(S,k);
%        d1 = [-2;-1;0];
%        d2 = [0;1;2];
%        L = spdiags(l,d1,k,k);
%        LF = full(L);
%        U = spdiags(u,d2,k,k);
%        UF = full(U);
%        IS = inv(UF)*inv(LF);
       %IS=inv(u)*inv(l);
       e1 = [1;zeros(k-2,1);1];
%        e2 = [0;1;zeros(k-2,1)];
%        en1 = [zeros(k-2,1);1;0];
%        en = [zeros(k-1,1);1];
       x1 = LUsolve(l,u,e1,k);
%        x2 = LUsovlve(l,u,e2,k);
%        x3 = LUsovlve(l,u,en1,k);
%        x4 = LUsovlve(l,u,en,k);
        %IS;
        timetable4(dimension,t) = toc;
         
        
    end
end
