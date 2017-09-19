repeat=5;
numscale=3;
timetable=zeros(numscale,repeat);
errortable= zeros(numscale,repeat);
scaletable = logspace(3,5,numscale);
for cccc = 1:numscale
    for cc = 1:repeat
        cccc
        cc
        %% This compute the pentadiagonal inverse problem by using the algorithm of solving the Toplizi matrices.
        %%first you need to enter the element of the pentadiagonal matrices.
        tic
         a = 1;
         b = 2;
         c = 2;
         d = 3;
         e = 3;
        %%ENTER THE dimension of the matrix to choose function
          k = scaletable(cccc);
          I = diag(ones(k,1));
          if rem(k,2) == 1
              W = pentadiag_inverseodd(a,b,c,d,e,k);
          else
              W = pentadiag_inversetest(a,b,c,d,e,k);
          end
             % W;
         timetable(cccc,cc)=toc;
         %% this is a test to check our codes
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
        S;
        IS=inv(S);
        IS;

        err = W*S-I;
        error=norm(err,'fro');
        errortable(cccc,cc)=error;

        %rel=abs(err)./abs(IS);
        %max(max(abs(rel)))
    end
end