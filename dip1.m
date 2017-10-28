% function [ out ] = downscaling(in,rx,ry)
    %intialization step
    clear all;close all;clc
    
    in=imread('pepper.bmp');
     in=imresize(in,1/32);
    rx=2;ry=2;
    [wi, hi, di] = size(in);
    wo = floor(wi/rx);ho = floor(hi/ry);
    n = wo*ho;
    out = zeros(wo,ho);
    sigma = cell(1,n);
    mu = cell(1,n);
    cpmu=cell(1,n);
    vk = cell(1,n);
    xk = zeros(n);
    yk = zeros(n);
    sig = ones(n);
    for i = 1:n
        sigma{i} = [rx/3, 0; 0, ry/3 ];
        yk(i) = mod(i-1,wo)+1;
        xk(i) = floor(i/wo)+1;
        mu{i}  = [xk(i), yk(i)].';
        cpmu{i}=mu{i};
        vk{i}  = [ 1/2 1/2 1/2 ].';
        sig(i) = 0.1;
    end
 
vyas=1;
while vyas>=1
%expectation step
cin =zeros(wi,hi);
for i=1:wo
    for j=1:ho
        cin(1+ry*(i-1),1+rx*(j-1))=(i-1)*wo+j;
    end
end
R=cell(1,n);
for i=1:wi
    for j=1:hi
        if(cin(i,j)>0)
            k=cin(i,j);
            for x=i-2*rx+1:i+2*rx-1
                for y=j-2*ry+1:j+2*ry-1
                    if(x>=1 && y>=1 && x<=wi && y<=hi)
                        R{k}=[R{k}, x+(y-1)*wi];
                    end
                end
            end
        end
    end
end
w=zeros(wo*ho,wi*hi);
for k=1:n
    for i=1:size(R{k},2)
        pi = [floor((R{k}(i)-1)/wi)+1, mod(R{k}(i)-1,wi)+1].';        
        temp_ci = rgb2lab(in(floor((R{k}(i)-1)/wi) +1, mod(R{k}(i)-1,wi)+1,:));
        tmp=pi-mu{k};
        ci = [temp_ci(:,:,1)/100,(temp_ci(:,:,2)+87)/186,(temp_ci(:,:,3)+108)/203].';
        lhv= -1/2*(tmp.')*inv(sigma{k})*tmp;
        rhv= ((ci(1,1)-vk{k}(1)).^2+ (ci(2,1)-vk{k}(2)).^2 +(ci(3,1)-vk{k}(3)).^2)/(2*sig(k)*sig(k));
         w(k,i)=exp(lhv - rhv);
    end
    wsum=sum(w(k,:));
 
    for i=1:size(R{k},2)
        w(k,i)=w(k,i)/wsum;
    end
end

 gamma=zeros(n,wi*hi);
 for i=1:wi*hi
     for k=1:n
         if ismember(i,R{k}(:))
             if(sum(w(:,i)~=0))
                gamma(k,i)=w(k,i)/sum(w(:,i));
             end
         end
     end
 end
 display "ayush"
 %Maximization step
 for k=1:n
     wsum=sum(gamma(k,:));
       if(wsum==0)
       display "Hi"
        break;
        end
     num1=0;
     num2=0;
     num3=0;
%      for i=wi*hi
     for i=1:size(R{k},2)
         pi = [floor((R{k}(i)-1)/wi)+1, mod(R{k}(i)-1,wi)+1].';        
         temp_ci = rgb2lab(in(floor((R{k}(i)-1)/wi) +1, mod(R{k}(i)-1,wi)+1,:));
         tmp=pi-mu{k};
         ci = [temp_ci(:,:,1)/100,(temp_ci(:,:,2)+87)/186,(temp_ci(:,:,3)+108)/203].';
         num1=num1+gamma(k,i)*tmp*(tmp.');
         num2=num2+gamma(k,i)*pi;
         num3=num3+gamma(k,i)*ci;
     end
     sigma{k}=num1/wsum;
     mu{k}=num2/wsum;
     vk{k}=num3/wsum;
 end
 %Correction step
 
 %Spatial constraints
 for k=1:n
      x0=mu{k}(1);y0=mu{k}(2);
      n4(1,1)=x0+1;n4(1,2)=y0+1;
      n4(2,1)=x0+1;n4(2,2)=y0-1;
      n4(3,1)=x0-1;n4(3,2)=y0+1;
      n4(4,1)=x0-1;n4(4,2)=y0-1;
      uk=0;
      sumx=0;sumy=0;
      for i=1:4
         if (n4(i,1)<=wo && n4(i,2)<=ho && n4(i,1)>0 && n4(i,2)>0 )     
             sumx=sumx+n4(i,1);sumy=sumy+n4(i,2);
         end
      end
      modi=sqrt(sumx*sumx + sumy*sumy);
      mubar{k}=[sumx sumy]/modi;
 end
  
 for k=1:n
   if(mu{k}(1)/2 + mubar{k}(1)/2 < xk(k)+rx/4)
     mu{k}(1) = xk(k)+rx/4;
   end
   if(mu{k}(1)/2 + mubar{k}(1)/2 > xk(k)-rx/4)
     mu{k}(1) = xk(k)-rx/4;
   end
   if(mu{k}(2)/2 + mubar{k}(2)/2 < yk(k)+ry/4)
     mu{k}(2) = yk(k)+ry/4;
   end
   if(mu{k}(2)/2 + mubar{k}(2)/2 > yk(k)-ry/4)
     mu{k}(2) = yk(k)-ry/4;
   end
 end
 
      
         
%Constrain spatial variance
 for k=1:n
     xmax=0.1;xmin=0.05;
     [U,S,V]=svd(sigma{k});
     if(S(1,1)<0.05)
         S(1,1)=0.05;
     end
     if(S(1,1)>0.1)
         S(1,1)=0.1;
     end
     
     
      if(S(2,2)<0.05)
         S(2,2)=0.05;
     end
     if(S(2,2)>0.1)
         S(2,2)=0.1;
     end
 
     sigma{k}=U*S*V;
 end
 %Shape constraints
 
 for k=1:n
     x0=cpmu{k}(1);y0=cpmu{k}(2);
     n8(1,1)=x0-1;n8(1,2)=y0-1;
     n8(2,1)=x0-1;n8(2,2)=y0;
     n8(3,1)=x0-1;n8(3,2)=y0+1;
     n8(4,1)=x0;n8(4,2)=y0-1;
     n8(5,1)=x0;n8(5,2)=y0+1;
     n8(6,1)=x0+1;n8(6,2)=y0-1;
     n8(7,1)=x0+1;n8(7,2)=y0;
     n8(8,1)=x0+1;n8(8,2)=y0+1;
      s=0;
     for ngbr =1:8
         if(n8(ngbr,1)>0 && n8(ngbr,2)>0 && n8(ngbr,1)<wo+1 && n8(ngbr,2)<ho+1) 
             d=[n8(ngbr,1)-x0 , n8(ngbr,2)-y0].';
             for i=1:size(R{k},2)
                  pi = [floor((R{k}(i)-1)/wi)+1, mod(R{k}(i)-1,wi)+1].';
                  tmp=[pi-mu{k}];
                  s=s+gamma(k,i)*max(0,(tmp.')*d).^2;
             end
             n1=wo*(n8(ngbr,2)-1)+n8(ngbr,1);
             f=0;
             for i=1:size(R{k},2)
                 f=f+gamma(k,i)*gamma(n1,i);
             end 
             if ( s>0.2*rx || f<0.08)
                 sig(k)=1.1*sig(k);
                 sig(n)=1.1*sig(n);
             end
         end
     end
 end
 
 %Break condition
  if(vyas>1)
        
        diff = 0.035
         flag=0;
          %Difference in M
          for k=1:n
            for i=1:2
                if(abs(cpymu{k}(i,1)-mu{k}(i,1))>1)
                   display "Hello1"
                    flag=1;break;
                end
            end
            for i=1:2
              if(abs(cpysigma{k}(i,1)-sigma{k}(i,1))>diff || abs(cpysigma{k}(i,2)-sigma{k}(i,2))>diff )
                   display "Hello2"
                flag=1;break;
              end
            end
            
            for i=1:3
              if(abs(cpyvk{k}(i,1)-vk{k}(i,1))>0.04)
                 display "cpyvk:"
                 abs(cpyvk{k}(i,1)-vk{k}(i,1))
                flag=1;break;
              end
            end
          end
          if(flag==0)
              break;
          end
 end
  cpymu=mu;cpysigma=sigma;cpyvk=vk;
 vyas=vyas+1
end

%%%%Main Code for final image%%%%%%

% for i=1:wi
%     for j=1:hi
%         if(cin(i,j)>0)
%             k=cin(i,j);
%             for x=i-2*rx+1:i+2*rx-1
%                 for y=j-2*ry+1:j+2*ry-1
%                     if(x>=1 && y>=1 && x<=wi && y<=hi)
%                         R{k}=[R{k}, x+(y-1)*wi];
%                     end
%                 end
%             end
%         end
%     end
% end
for clr=1:3
for k=1:wo*ho
    vyasum=0;
    for j=1:size(R{k},2)
          vyasum=vyasum + in(floor((R{k}(j)-1)/wi)+1, mod(R{k}(j)-1,wi)+1,clr)*w(k,j);
    end
    (k-1)/wo+1
    mod(k,wo)
    finalImage(floor((k-1)/wo)+1,mod(k,wo)+1, clr)=vyasum;
end
end
imshow(uint8(in));figure
imshow(uint8(finalImage));

               
                        
   