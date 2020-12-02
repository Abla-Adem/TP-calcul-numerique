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









