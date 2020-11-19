function [c]=matmat3b(a,b,n)
    c=zeros(n,n);
    
    for i=1:n
        for j=1:n
           for k=1:n
              c(i,j)=a(i,k)*b(k,j)+c(i,j);
           end 
        end
    end
    
endfunction
function [c]=deux_boucle(a,b,n)
    c=zeros(n,n);
    
    for i=1:n
        for j=1:n
           
              c(i,j)=a(i,:)*b(:,j)+c(i,j);
            
        end
    end
    
endfunction
function [c]=une_boucle(a,b,n)
    c=zeros(n,n);    
    
    for i=1:n
           
              c(i,:)=a(i,:)*b+c(i,:);
            
        
    end
    
endfunction
n=3;
a=rand(n,n);
b=rand(n,n);
t=zeros(1,10);
t1=zeros(1,20);
t2=zeros(1,10);
t3=zeros(1,10);
step=1;
for i=1:10
    n=n+10
    a=rand(n,n);
    b=rand(n,n);

    for j=1:10
        tic();
        matmat3b(a,b,n);
        t(step)=t(step)+toc();
        tic();
        deux_boucle(a,b,n);
        t1(step)=t1(step)+toc();
        tic();
        une_boucle(a,b,n);
        t2(step)=t2(step)+toc(); 
        tic()
        a*b;
        t3(step)=t3(step)+toc();  
    end
    t(step)=t(step)/10.;
    t1(step)=t1(step)/10;
    t2(step)=t2(step)/10;
step=step+1;
end
t
t1
t2
t3
// triu() matrice triangulaire superier
