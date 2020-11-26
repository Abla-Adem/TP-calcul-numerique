function [L,U]=LU(a,n)
    L=zeros(n,n);
    U=zeros(n,n);
    for i=1:n
        U(1,i)=a(1,i);
        L(i,i)=1;
    end
    for i=2:n
        L(i,1)=a(i,1)/U(1,1)
        for j=2:(i-1)
            L(i,j)=a(i,j);
            for k=1:(j-1)
            L(i,j)=L(i,j)-L(i,k)*U(k,j)
            end
            L(i,j)=L(i,j)/U(j,j);
        end
        U(i,i)=a(i,i);
        for k=1:(i-1)
                U(i,i)=U(i,i)-L(i,k)*U(k,i);
           end 
        for j=(i+1):n
            U(i,j)=a(i,j);
           for k=1:(i-1)
                U(i,j)=U(i,j)-L(i,k)*U(k,j);
           end 
        end
    end
endfunction
function [x] = usolve(U,b,n)
    x=zeros(1,n)
    j=n-1
    x(n)=b(n)/U(n,n);
    for i=1:n-1
        x(j)=b(j);
        for k=(j+1):n
            x(j)=x(j)-x(k)*U(j,k); 
        end
        x(j)=x(j)/U(j,j);
        j=j-1
    end
endfunction
function [x] = lsolve(L,b,n)
    x=zeros(1,n)
    j=n-1
    x(1)=b(1);
    for i=2:n
        x(i)=b(i);
        for k=1:(i-1)
            x(i)=x(i)-x(k)*L(i,k); 
        end
        x(i)=x(i)/L(i,i);
    end
endfunction
