function [L,U] = Mylu2(A,n)

    for k=1:n-1
        for i=(k+1):n
            A(i,k)=A(i,k)/A(k,k)
            A(i,k+1:n)=A(i,k+1:n)-A(i,k)*A(k,k+1:n);
         end
    end
   L=tril(A,-1)+eye(n,n)
   U=triu(A)
endfunction
function [P,L,U] = Mylu1(A,n)
    [P,A]=pivot(A,n)
    for k=1:n-1
            A(k+1:n,k)=A(k+1:n,k)/A(k,k)
            A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-A(k+1:n,k)*A(k,k+1:n);
         
    end
   L=tril(A,-1)+eye(n,n)
   U=triu(A)
endfunction

function [P,A] = pivot(A,n)
    P=zeros(n,n)
    t=[1:n]
    tem=0
    z=zeros(1,n)
    for k=1:n-1
        [val,i]=max(abs(A(k:n,k)))
         i=i+k-1
         
         P(t(i),k)=P(t(i),k)+1;
         z(t(i))=z(t(i))-1   
            if(k~=i)
               temp=A(k,1:n)
               tem=t(i)        
               t(i)=t(k)
               t(k)=tem
               A(k,1:n)=A(i,1:n)
               A(i,1:n)=temp;
            end
            
    end
    [val,i]=max(z)
    P(i,n)=P(i,n)+1
endfunction
function [L,U]=LU_crout(a,n)
    L=zeros(n,n);
    U=zeros(n,n);
    for i=1:n
        U(1,i)=a(1,i);
        L(i,i)=1;
    end
    for i=2:n
        L(i,1)=a(i,1)/U(1,1);
        for j=2:(i-1)
            
            L(i,j)=a(i,j);
            L(i,j)=L(i,j)-L(i,1:(j-1))*U(1:(j-1),j)
            L(i,j)=L(i,j)/U(j,j);
        end
        U(i,i)=a(i,i);
        U(i,i)=U(i,i)-L(i,1:(i-1))*U(1:(i-1),i);   
        for j=(i+1):n
            U(i,j)=a(i,j);
            U(i,j)=U(i,j)-L(i,1:(i-1))*U(1:(i-1),j);
        end
    end
endfunction

