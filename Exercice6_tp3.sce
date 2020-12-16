function [A]=trimat(n)
    A=zeros(n,n)
    A(1,1)=2
    A(1,2)=-1
    for i=2:n-1
    A(i,i-1)=-1
    A(i,i)=2
    A(i,i+1)=-1
    end
    A(n,n-1)=-1
    A(n,n)=2
endfunction
function [x,test,r]=jacobi(A,n,b,eps)
    
    
    
    M=-(triu(a,1)+tril(a,-1))
    temp=diag(A)
    N=zeros(n,n)
    for i=1:n
        N(i,i)=temp(i)
    end
    
    
    
    x=zeros(n,1);
    B=inv(N)*M
    v=inv(N)*b
    temp=zeros(n,1)
    ev=zeros(n,1)
    r=1
    test=0
    //on considere que l'intervalle est uniforme ,ainsi on le calcul une seul fois
    while(r>eps)
        ev=zeros(n,1)
       temp=zeros(n,1)
       temp=B*x+v                
       
        
          
        r=norm(b-A*x)
        // e=norm(ev)
             x=temp
          test=test+1
          
         
    end
endfunction

function [x,test,r]=jacobi_tri(A,n,b,eps)
    
    D=(1./diag(A))
    F=-diag(A,1)(2:n-1)
    E=-diag(A,-1)(1:n-2)
    F=D(2:n-1).*F
    E=D(2:n-1).*E
    x=zeros(n,1);
    temp=zeros(n,1)
    ev=zeros(n,1)
    v=D.*b
    r=1
    test=0
    //on considere que l'intervalle est uniforme ,ainsi on le calcul une seul fois
    while(r>eps)
       //temp=B()*x+v
       for i=2:n-1
       temp(i)=F(i-1)*x(i-1)+E(i-1)*x(i+1)+v(i)
       end
       temp(1)=D(1)*-A(1,2)*x(2)+v(1)
       temp(n)=D(n)*-A(n,n-1)*x(n-1)+v(n)
                  
        r=norm(b-A*x)
        // e=norm(ev)
             x=temp
          test=test+1
          
         
    end
endfunction



function [x,test,e]=gaus_seidel(A,n,b,eps)
    F=-triu(a,1)
    temp=diag(A)
    D=zeros(n,n)
    for i=1:n
        D(i,i)=temp(i)
    end
    
    
    
    x=zeros(n,1);
    B=inv(D+tril(a,-1))
    v=B*b
    temp=zeros(n,1)
    ev=zeros(n,1)
    e=1
    test=0
    //on considere que l'intervalle est uniforme ,ainsi on le calcul une seul fois
    while(e>eps)
      
       temp=B*F*x+v                
       
        
          ev=b-temp 
        
         
             x=temp
         e=norm(b-A*x) 
         test=test+1
    end
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
