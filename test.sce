function [err_a,err_ar,err,cnd] = TestExo2(a,c)
    n=c
    t=zeros(1,a)
    err_a=zeros(1,(a));
    err_ar=zeros(1,(a));
    cnd=zeros(1,a);
    err=zeros(1,(a));
    for i=1:a
        t(i)=t(i)+i
        a=rand(n,n);
        a1=a;
        xex=rand(n,1);
        b=a*xex;
        x=inv(a)*b;        
        
        cnd(i)=cond(a);
        err_a(i)=norm(x-xex)/norm(xex);
                
        err_ar(i)=norm(b-a*x)/norm(b);
        err(i)=cnd(i)*err_ar(i);
        n=n+10
     end  
    err_a=err_a';
    err_ar=err_ar';
    err=err';
    plot2d(t,[err_a err_ar  ],style=[2,3]);
    legends(['erreur avant ';'erreur arriere '],[2 3],opt='lr')

endfunction




function [t,t1,t2,t3]=TestExo3(a,n)
a=rand(n,n);
x=zeros(n);
b=rand(n,n);
t=zeros(1,10);
t1=zeros(1,10);
t2=zeros(1,10);
t3=zeros(1,10);
step=1;
for i=1:10
    n=n+10
    x(i)=x(i)+n;
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
    t1(step)=t1(step)/10.;
    t2(step)=t2(step)/10.;
    t3(step)=t3(step)/10.;
step=step+1;
end
t=t';t1=t1';t2=t2';t3=t3';
plot2d(x,[t t1 t2 t3],style=[-1,2,3,4]);
legends(['3 boucle';'2 boucle';'boucle';'direct'],[-1,2 3,4],opt='lr')

endfunction
function [err_a,err_ar,err,cnd] = TestExo4(a,c,m)
    n=c
    t=zeros(1,a)
    err_a1=zeros(1,(a));
    err_ar1=zeros(1,(a));
    cnd1=zeros(1,(a));
    err1=zeros(1,(a));
    err_a=zeros(1,(a));
    err_ar=zeros(1,(a));
    cnd=zeros(1,a);
    err=zeros(1,(a));
    for i=1:a
        t(i)=t(i)+i
        if(m==0)
        a=tril(rand(n,n));
        xex=rand(n,1);
        b=a*xex;
        temp=lsolve(a,b,n);
        cnd(i)=cond(a);
        err_a(i)=norm(temp-xex)/norm(xex);
        err_ar(i)=norm(b-a*temp)/norm(b);
        err(i)=cnd(i)*err_ar(i);
        else
        a=triu(rand(n,n));
        xex=rand(n,1);
        b1=a*xex;
        x=usolve(a,b1,n);        
        
        cnd(i)=cond(a);
        err_a(i)=norm(x-xex)/norm(xex);
                
        err_ar(i)=norm(b1-a*x)/norm(b1);
        err(i)=cnd(i)*err_ar(i);
        end
        n=n+10
     end  
        
    
    
    err_a=err_a';
    err_ar=err_ar';
    err=err';
    plot2d(t,[err_a err_ar err ],style=[-1,2,3]);
    legends(['erreur avant ';'erreur arriere ';'erreur '],[-1,2 3],opt='lr')

endfunction



function [err_a,err_ar,err,cnd] = TestExo5(a,c)
    n=c
    t=zeros(1,a)
    err_a=zeros(1,(a));
    err_ar=zeros(1,(a));
    cnd=zeros(1,a);
    err=zeros(1,(a));
    for i=1:a
        
        a=rand(n,n);
        a1=a;
        xex=rand(n,1);
        b=a*xex;
        b1=a*xex;
        x=Gauss(a1,b1,n);        
        
        cnd(i)=cond(a);
        err_a(i)=norm(x-xex)/norm(xex);
                
        err_ar(i)=norm(b-a*x)/norm(b);
        err(i)=cnd(i)*err_ar(i);
        n=n+10
     end  
    err_a=err_a';
    err_ar=err_ar';
    err=err';
    plot2d(t,[err_a err_ar err ],style=[-1,2,3]);
    legends(['erreur avant ';'erreur arriere ';'erreur '],[-1,2 3],opt='lr')

endfunction

function [err_a,err_ar,err,cnd] = TestExo6(a,c)
    n=c
    t=zeros(1,a)
    err_a1=zeros(1,(a));
    err_ar=zeros(1,(a));
    cnd=zeros(1,a);
    err=zeros(1,(a));
    err1=zeros(1,(a));
    
    for i=1:a
        t(i)=t(i)+i
        a=rand(n,n);
        cnd(i)=cond(a);
        [L1,U1]=LU_crout(a,n)
        [L,U]=Mylu1(a,n)        
        
        
        err_ar1(i)=norm(L1*U1-a)/norm(a);
        err1(i)=cnd(i)*err_ar1(i);
        
        err_ar(i)=norm(L*U-a)/norm(a);
        err(i)=cnd(i)*err_ar(i);
        n=n+10
     end  
    err_ar1=err_ar1';
    err_ar=err_ar';
    err=err';
    err1=err1';
    plot2d(t,[ err  err1 ],style=[2,4]);
    legends(['erreur LU';'erreur LU_crout'],[2 ,4],opt='lr')

endfunction

function [err_ar1,err_ar,err,err1] = TestExo6_2(a,c)
    n=c
    t=zeros(1,a)
    err_ar1=zeros(1,(a));
    err_ar=zeros(1,(a));
    cnd=zeros(1,a);
    err=zeros(1,(a));
    err1=zeros(1,(a));
    
    for i=1:a
        t(i)=t(i)+i
        a=rand(n,n);
        cnd(i)=cond(a);
        [L1,U1]=lu(a)
        [L,U]=Mylu1(a,n)        
        
        
        err_ar1(i)=norm(L1*U1-a)/norm(a);
        err1(i)=cnd(i)*err_ar1(i);
        
        err_ar(i)=norm(L*U-a)/norm(a);
        err(i)=cnd(i)*err_ar(i);
        n=n+10
     end  
    err_ar1=err_ar1';
    err_ar=err_ar';
    err=err';
    err1=err1';
    plot2d(t,[ err  err1 ],style=[2,4]);
    legends(['erreur LU';'erreur LU_scilab'],[2 ,4],opt='lr')

endfunction
