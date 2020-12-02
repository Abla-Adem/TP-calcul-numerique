

function [x] = Gauss_tri(A,b,n)
    m=zeros(2,1)
    m1=A(2,1)/A(1,1)
    b(2)=b(2)-m1*b(1)
    A(2,2)=A(2,2)-m1*A(1,2);
    for k=2:n-2
        /*m(1:2)=m(1:2)+A(k+1:k+2,k)/A(k,k)
        b(k+1:k+2)=b(k+1:k+2)-m(1:2)*b(k)
        A(k+1:k+2)=A(k+1:k+2)-A(k,k+1:k+2)*m(1:2)
        m(1:2)=m(1:2)-m(1:2)*/
        //m=A(k+1,k)/A(k,k)
        A(k+1,k)=A(k+1,k)/A(k,k)
        
        b(k+1)=b(k+1)-A(k+1,k)*b(k)
        A(k+1,k+1)=A(k+1,k+1)-A(k+1,k)*A(k,k+1)
        A(k+1,k+2)=A(k+1,k+2)-A(k+1,k)*A(k,k+2)
        
        //m1=A(k+2,k)/A(k,k)
        A(k+2,k)=A(k+2,k)/A(k,k)
        b(k+2)=b(k+2)-A(k+2,k)*b(k)
        A(k+2,k+1)=A(k+2,k+1)-A(k+2,k)*A(k,k+1)
        A(k+2,k+2)=A(k+2,k+2)-A(k+2,k)*A(k,k+2)
        
           
    end
    A(n,n-1)=A(n,n-1)/A(n-1,n-1)
    b(n)=b(n)-A(n,n-1)*b(n-1)
    A(n,n)=A(n,n)-A(n,n-1)*A(n-1,n)
            
    x=usolve(A,b,n)
endfunction


function [L,U] = lu_tri(A,n)
    m=zeros(2,1)
    m1=A(2,1)/A(1,1)
    b(2)=b(2)-m1*b(1)
    A(2,2)=A(2,2)-m1*A(1,2);
    for k=2:n-2
        /*m(1:2)=m(1:2)+A(k+1:k+2,k)/A(k,k)
        b(k+1:k+2)=b(k+1:k+2)-m(1:2)*b(k)
        A(k+1:k+2)=A(k+1:k+2)-A(k,k+1:k+2)*m(1:2)
        m(1:2)=m(1:2)-m(1:2)*/
        //m=A(k+1,k)/A(k,k)
        A(k+1,k)=A(k+1,k)/A(k,k)
        
        b(k+1)=b(k+1)-A(k+1,k)*b(k)
        A(k+1,k+1)=A(k+1,k+1)-A(k+1,k)*A(k,k+1)
        A(k+1,k+2)=A(k+1,k+2)-A(k+1,k)*A(k,k+2)
        
        //m1=A(k+2,k)/A(k,k)
        A(k+2,k)=A(k+2,k)/A(k,k)
        b(k+2)=b(k+2)-A(k+2,k)*b(k)
        A(k+2,k+1)=A(k+2,k+1)-A(k+2,k)*A(k,k+1)
        A(k+2,k+2)=A(k+2,k+2)-A(k+2,k)*A(k,k+2)
        
           
    end
    A(n,n-1)=A(n,n-1)/A(n-1,n-1)
    b(n)=b(n)-A(n,n-1)*b(n-1)
    A(n,n)=A(n,n)-A(n,n-1)*A(n-1,n)
    
    L=tril(A,-1)+eye(n,n)
    U=triu(A)
endfunction
