
function [x] = usolve(U,b,n)
    x=zeros(n,1)
    j=n-1
    x(n)=b(n)/U(n,n);
    for i=1:n-1 
        U(j,j+1:n)=U(j,j+1:n)/U(j,j)
        x(j)=b(j)/U(j,j)-U(j,j+1:n)*x(j+1:n);
        j=j-1
    end
endfunction
function [x] = lsolve(L,b,n)
    x=zeros(n,1)
    j=n-1
    x(1)=b(1)/L(1,1)
    for i=2:n
        x(i)=b(i);
        L(i,1:i-1)=L(i,1:i-1)/L(i,i)
        x(i)=b(i)/L(i,i)-L(i,1:i-1)*x(1:i-1)
        
    end
endfunction
