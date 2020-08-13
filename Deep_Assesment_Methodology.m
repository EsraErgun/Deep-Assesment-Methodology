clc;
clear all;
P0=[
3007.123445
3066.562869
3243.843078
3374.515171
3573.941185
3827.52711
4146.316646
4336.426587
4695.92339
5032.144743
5234.296666
5609.3826
6094.01799
6726.358956
7225.69136
7801.456664
8592.253537
9452.576519
10564.94822
11674.18631
12574.79151
13976.10975
14433.78773
15543.89372
17121.22548
18236.82773
19071.22719
20038.9411
21417.01193
22857.15443
23888.60001
24342.2589
25418.99078
26387.29373
27694.85342
28690.8757
29967.71272
31459.139
32853.67698
34513.5615
36334.90878
37133.24281
38023.16111
39496.48588
41712.80107
44114.74778
46298.73144
47975.96768
48382.55845
47099.98047
48466.82338
49883.11398
51603.49726
53106.90977
55032.958
56803.47243
57904.20196
59927.92983
62794.58565

]






m=58;
M=3;
P=zeros(m,1);
 

x=linspace(1,m,m);
for extr=1:1
    if (-1)^m==1
    t=m/2;
    else
      t=(m-1)/2;  
    end
    t=10;
    l=t-1;



if extr==1

 for i=1:m
     P(i)=P0(i);
 end





finalF=zeros(m,1);

else
    for iq=1:m-1
  %P(iq)= x(iq)+2*sin(x(iq)*pi/5);
  
    end
    P(m)=F11(m);
   %P(m)= 32903.53;
 % x=linspace(1,NOD,num1);
    %P(iq)= Psm(iq); 
     %P(iq)=x(iq)+2*sin(x(iq)*pi/5);
finalF=zeros(m,1);
end







MD=1+l*M+l*M;
%MD=10;

A=zeros(MD,MD);
C=zeros(MD,1);
%niu=0.01;
NONIU=100;
niu1=linspace(0.01,1,NONIU);

%fileID = fopen('Error.dat','w');
%fprintf(fileID,'%6s\n','fError');

%fprintf(fileID,'%6s\n','x');
%fprintf(fileID,'%6s\n','AError');

%for i=1:Nxx

    
     
min=100000;
for i1=1:NONIU
    niu=niu1(i1);
 
    
    %%%%%First ROW%%%%%%%%%%%%%%%%
    A(1,1)=m-t+1;
   
    icol=1;
      for k=1:l
    for n=1:M
      
            icol=icol+1;
            A(1,icol)=0;
            for i=t:m
               
              A(1,icol)=A(1,icol)+(i-k).^(niu+n-1)  ;
               
            end
            A(1,icol)=A(1,icol)*gamma(n+1)/gamma(n+niu);
        end
      end
    
      
      icol=1+l*M;
       for k=1:l
    for n=1:M
      
            icol=icol+1;
            A(1,icol)=0;
            for i=t:m
               
              A(1,icol)=A(1,icol)+(i-k).^(niu+n-2)  ;
               
            end
            A(1,icol)=A(1,icol)*gamma(n+1)/gamma(n+niu-1);
        end
      end
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    

    
    %%%%%%%%% first scolomn%%%%%%%%%%%%%%%%
       irow=1;
        for k=1:l
    for n=1:M
       
            irow=irow+1;
            A(irow,1)=0;
            for i=t:m
              
              A(irow,1)=A(irow,1)+(i-k).^(niu+n-1)  ;
               
            end
            A(irow,1)=A(irow,1)*gamma(n+1)/gamma(n+niu);
        end
    end
    
       irow=1+l*M;
        for k=1:l
    for n=1:M
       
            irow=irow+1;
            A(irow,1)=0;
            for i=t:m
              
              A(irow,1)=A(irow,1)+(i-k).^(niu+n-2)  ;
               
            end
            A(irow,1)=A(irow,1)*gamma(n+1)/gamma(n+niu-1);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  
  
    
    %%%%%%%%%%%%%%%%%%%%%first kvadrant%%%%%%%%%%
    
        irow=1;
         for kpr=1:l
       for npr=1:M
        
            irow=irow+1;
            icol=1;
         for k=1:l    
    for n=1:M
      
            icol=icol+1;
            A(irow,icol)=0;
            for i=t:m
              
              A(irow,icol)=A(irow,icol)+(i-k).^(niu+n-1)*(i-kpr).^(niu+npr-1) ;
              
            end
            A(irow,icol)=A(irow,icol)*gamma(n+1)/gamma(n+niu)*gamma(npr+1)/gamma(npr+niu);
        end
    end
        end
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    
    
    
    
     %%%%%%%%%%%%%%%%%%%%%second kvadrant%%%%%%%%%%
    
        irow=1;
         for kpr=1:l
       for npr=1:M
        
            irow=irow+1;
            icol=1+l*M;
         for k=1:l    
    for n=1:M
      
            icol=icol+1;
            A(irow,icol)=0;
            for i=t:m
              
              A(irow,icol)=A(irow,icol)+(i-k).^(niu+n-2)*(i-kpr).^(niu+npr-1) ;
              
            end
            A(irow,icol)=A(irow,icol)*gamma(n+1)/gamma(n+niu-1)*gamma(npr+1)/gamma(npr+niu);
        end
    end
        end
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    
      %%%%%%%%%%%%%%%%%%%%%third kvadrant%%%%%%%%%%
    
        irow=1+l*M;
         for kpr=1:l
       for npr=1:M
        
            irow=irow+1;
            icol=1;
         for k=1:l    
    for n=1:M
      
            icol=icol+1;
            A(irow,icol)=0;
            for i=t:m
              
              A(irow,icol)=A(irow,icol)+(i-k).^(niu+n-1)*(i-kpr).^(niu+npr-2) ;
              
            end
            A(irow,icol)=A(irow,icol)*gamma(n+1)/gamma(n+niu)*gamma(npr+1)/gamma(npr+niu-1);
        end
    end
        end
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    
     %%%%%%%%%%%%%%%%%%%%% forth kvadrant%%%%%%%%%%
    
        irow=1+l*M;
         for kpr=1:l
       for npr=1:M
        
            irow=irow+1;
            icol=1+l*M;
         for k=1:l    
    for n=1:M
      
            icol=icol+1;
            A(irow,icol)=0;
            for i=t:m
              
              A(irow,icol)=A(irow,icol)+(i-k).^(niu+n-2)*(i-kpr).^(niu+npr-2) ;
              
            end
            A(irow,icol)=A(irow,icol)*gamma(n+1)/gamma(n+niu-1)*gamma(npr+1)/gamma(npr+niu-1);
        end
    end
        end
       end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    
    
    
   %%%%%%%%%%%%%%%%%% C matrica%%%%%%%%%%%%%%%%%%%%%
   C(1)=0;
   for i=t:m
   C(1)=C(1)+P(i);    
   end
   
  
   
   icol=1;
    for k=1:l
   for n=1:M
      
           icol=icol+1;
          C(icol)=0;
   for i=t:m
      
   C(icol)=C(icol)+P(i)*(i-k).^(niu+n-1);  
     
   end  
       C(icol)=C(icol)*gamma(n+1)/gamma(n+niu)  ;
       end
    end
   
   
      icol=1+l*M;
    for k=1:l
   for n=1:M
      
           icol=icol+1;
          C(icol)=0;
   for i=t:m
      
   C(icol)=C(icol)+P(i)*(i-k).^(niu+n-2);  
     
   end  
        C(icol)=C(icol)*gamma(n+1)/gamma(n+niu-1);     
       end
    end
    
    
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
   
%B=zeros(1,MD);
%Bfinal=zeros(1,MD);
B=A\C;
% q=0;
% for i=1:MD
%     q=q+A(1,i)*B(i);
% end
% q=q-C(1)
F=zeros(m,1);
eps1=zeros(m,1);
for i=t:m
   F(i)=B(1);
  
   irow=1;
      for k=1:l
   for n=1:M
    
   irow=irow+1;
  
    F(i)=F(i)+gamma(n+1)/gamma(n+niu)*(i-k).^(niu+n-1)*B(irow);
 
       end
      end
   
      
        irow=l*M+1;
      for k=1:l
   for n=1:M
    
   irow=irow+1;
  
    F(i)=F(i)+gamma(n+1)/gamma(n+niu-1)*(i-k).^(niu+n-2)*B(irow);
 
       end
      end
end
%figure(1)
%plot(F)
%figure(2)
%hold on;
%plot(P)
%eps=0;
dat=0;
%eps=zeros(1,1)
eps=0;
for i=1:t-1
   F(i)=P(i);
end
for i=t:m
  %  dat=dat+P(i).^2;
 % eps=eps+(P(i)-F(i)).^2/dat*100;
eps=eps+abs((P(i)-F(i))/P(i));
end
eps=eps*100/m;
%eps1(i1)=eps;
%eps
%fprintf(fileID,'%12.8f %12.8f \n',niu,eps  );
if(eps<min)
    min=eps;
    niumin=niu;
    minerror=min;
    finalF=F;
    Bfinal=B;
end

%end
%niu
%minerror
end
minerror
niumin
%fclose(fileID);
%e%ps1

 %figure(1)
% plot(finalF)
% 
%hold on;
% plot(P)
 %grid on
 %set(gca,'xtick',[1:1:21]);
%format long
%B


%for id=1:NOD
  
 
  
   m=m+1;
  % num1=num;
  % xsm=linspace(1,num,num1);
   x=linspace(1,m,m);
    Psm=zeros(m,1);
   F11=zeros(m,1);
   for i=1:m
    
   
   Psm(i)=P0(i);
   end

   for i=1:t-1
   %F11(i)=Psm(i);
    F11(i)=P(i);
   end

   
   
   for i=t:m
     
   F11(i)=Bfinal(1);
   irow=1;
      for k=1:l
   for n=1:M
    
   irow=irow+1;
  
    F11(i)=F11(i)+gamma(n+1)/gamma(n+niumin)*(i-k).^(niumin+n-1)*Bfinal(irow);
 
       end
      end
   
      
        irow=l*M+1;
      for k=1:l
   for n=1:M
    
   irow=irow+1;
  
  
    F11(i)=F11(i)+gamma(n+1)/gamma(n+niumin-1)*(i-k).^(niumin+n-2)*Bfinal(irow);
 
       end
      end
   
   end

end
 




  plot(x,F11)
   hold on;
   plot(x,Psm)
 abs((Psm(m)-F11(m))/Psm(m))*100

