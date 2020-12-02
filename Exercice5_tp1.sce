function [x] = Gauss(A,b,n)

    for k=1:n-1
        for i=(k+1):n
            m=A(i,k)/A(k,k)
            b(i)=b(i)-m*b(k)
            A(i,k+1:n)=A(i,k+1:n)-m*A(k,k+1:n);
         end
    end
    x=usolve(A,b,n)
endfunction













