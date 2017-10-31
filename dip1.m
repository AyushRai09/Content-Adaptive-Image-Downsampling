% function [ out ] = downscaling(in,rx,ry)
    %intialization step
    clear all;close all;clc
    
    in=imread('test1.png');
    in=imresize(in,1/4);
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
        mu{i}  = [xk(i)+1/2, yk(i)+1/2].';
        cpmu{i}=[xk(i) yk(i)];
        vk{i}  = [ 1/2 1/2 1/2 ].';
        sig(i) = 0.0001;
    end
 
vyas=1;
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
                        R{k}=[R{k}, (x-1)*wi+y];
                    end
                end
            end
        end
    end
end
w=zeros(wo*ho,wi*hi);
gamma=zeros(n,wi*hi);

while vyas>=1
    
%expectation step
for k=1:wo*ho
    for i=1:wi*hi
        if(ismember(i,R{k}(:))==1)
            h1=floor((i-1)/wi)+1;h2=mod(i-1,wi)+1;
            pi = [h1, h2].';        
            temp_ci = rgb2lab(in(h1,h2,:));
            tmp=pi-mu{k};
            ci = [temp_ci(:,:,1)/100,(temp_ci(:,:,2)+87)/186,(temp_ci(:,:,3)+108)/203].';
            lhv= -1/2*(tmp.')*inv(sigma{k})*tmp
            rhv= ((ci(1,1)-vk{k}(1)).^2+ (ci(2,1)-vk{k}(2)).^2 +(ci(3,1)-vk{k}(3)).^2)/(2*sig(k)*sig(k))
            if(isfinite(exp(lhv-rhv))==1)
                w(k,i)=exp(lhv - rhv);
            end
        end
    end
    tempsum=0;
    for i=1:wi*hi
        if(ismember(i,R{k})==1)
            tempsum=tempsum+w(k,i);
        end
    end
    wsum=tempsum;
    for i=1:wi*hi
        if(isfinite(w(k,i)/wsum) == 1 && ismember(i,R{k}(:))==1)
            w(k,i)=w(k,i)/wsum;
        end
    end
end

 for i=1:wi*hi
     for k=1:wo*ho
         if (ismember(i,R{k}(:)) && sum(w(:,i))~=0)
%                  w(k,i)
%                  sum(w(:,i))
%                  k
%                  i
                gamma(k,i)=w(k,i)/sum(w(:,i));
         end
     end
 end

 %Maximization step
 for k=1:wo*ho
     wsum=sum(gamma(k,:));
     num1=[0 0; 0 0];
     num2=[0 ; 0];
     num3=[0; 0; 0];
%      for i=wi*hi
     for i=1:wi*hi
         h1=floor((i-1)/wi)+1;h2=mod(i-1,wi)+1;
         pi = [h1, h2].';        
         temp_ci = rgb2lab(in(h1,h2,:));
         tmp=pi-mu{k};
         ci = [temp_ci(:,:,1)/100,(temp_ci(:,:,2)+87)/186,(temp_ci(:,:,3)+108)/203].';
         num1=num1+gamma(k,i)*tmp*(tmp.');
         num2=num2+gamma(k,i)*pi;
         num3=num3+gamma(k,i)*ci;
     end
     if(wsum~=0)
         sigma{k}=num1/wsum;
         mu{k}=num2/wsum;
         vk{k}=num3/wsum;
     end
 end
 %Correction step
 
 %Spatial constraints
 for k=1:wo*ho
      x0=cpmu{k}(1);y0=cpmu{k}(2);
      n4(1,1)=x0+1;n4(1,2)=y0+1;
      n4(2,1)=x0+1;n4(2,2)=y0-1;
      n4(3,1)=x0-1;n4(3,2)=y0+1;
      n4(4,1)=x0-1;n4(4,2)=y0-1;
      uk=0;
      sumx=0;sumy=0;
      mutmp=[0;0];
      for i=1:4
           if (n4(i,1)<=wo && n4(i,2)<=ho && n4(i,1)>0 && n4(i,2)>0 )
               vlk=(n4(i,1)-1)*wo+n4(i,2);
               mutmp=mutmp+mu{vlk};
               sumx=sumx+n4(i,1);sumy=sumy+n4(i,2);
           end
      end

      modi=sqrt(sumx*sumx + sumy*sumy);
      if(isfinite(mutmp/modi)==1)
        mubar{k}=mutmp/modi;
      end
 end
  
 for k=1:wo*ho
   if(mu{k}(1)/2 + mubar{k}(1)/2 < xk(k)+rx/4)
     mu{k}(1) = xk(k)+rx/4;
   elseif(mu{k}(1)/2 + mubar{k}(1)/2 > xk(k)-rx/4)
     mu{k}(1) = xk(k)-rx/4;
   end
   if(mu{k}(2)/2 + mubar{k}(2)/2 < yk(k)+ry/4)
     mu{k}(2) = yk(k)+ry/4;
   elseif(mu{k}(2)/2 + mubar{k}(2)/2 > yk(k)-ry/4)
     mu{k}(2) = yk(k)-ry/4;
   end
 end
 
      
         
%Constrain spatial variance
 for k=1:wo*ho
     xmax=0.1;xmin=0.05;
     [U,S,V]=svd(sigma{k});
     if(S(1,1)<0.05)
         S(1,1)=0.05;
     elseif(S(1,1)>0.1)
         S(1,1)=0.1;
     end    
     if(S(2,2)<0.05)
         S(2,2)=0.05;
     elseif(S(2,2)>0.1)
         S(2,2)=0.1;
     end
     sigma{k}=U*S*V;
 end
 %Shape constraints
 
 for k=1:wo*ho
     x0=cpmu{k}(1);y0=cpmu{k}(2);
     n8(1,1)=x0-1;n8(1,2)=y0-1;
     n8(2,1)=x0-1;n8(2,2)=y0;
     n8(3,1)=x0-1;n8(3,2)=y0+1;
     n8(4,1)=x0;n8(4,2)=y0-1;
     n8(5,1)=x0;n8(5,2)=y0+1;
     n8(6,1)=x0+1;n8(6,2)=y0-1;
     n8(7,1)=x0+1;n8(7,2)=y0;
     n8(8,1)=x0+1;n8(8,2)=y0+1;
     
     for ngbr =1:8
         s=0;
         if(n8(ngbr,1)>0 && n8(ngbr,2)>0 && n8(ngbr,1)<wo+1 && n8(ngbr,2)<ho+1) 
             d=[n8(ngbr,1)-x0 , n8(ngbr,2)-y0].';
             for i=1:wi*hi
                  if(ismember(i,R{k}(:))==1)
                      pi = [floor((i-1)/wi)+1, mod(i-1,wi)+1].';
                      tmp=pi-mu{k};
                      s=s+gamma(k,i)*max(0,(tmp.')*d).^2;
                  end
             end
             n1=wo*(n8(ngbr,1)-1)+n8(ngbr,2);
             f=0;
             for i=1:wi*hi
                 if(ismember(i,R{k})==1)
                    f=f+gamma(k,i)*gamma(n1,i);
                 end
             end 
             if ( s>0.2*rx || f<0.08)
                 sig(k)=1.1*sig(k);
                 sig(n1)=1.1*sig(n1);
             end
         end
     end
 end
 
 %Break condition
    if(vyas>30)
        convergeFlag=1;
         for k=1:wo*ho
             if(abs(vk{k}(1,1)-cpyvk{k}(1,1))>0.002 || abs(vk{k}(2,1)-cpyvk{k}(2,1))>0.002 || abs(vk{k}(3,1)-cpyvk{k}(3,1))>0.002)
                 convergeFlag=0;
                 break;
             end
         end
         if(convergeFlag)
                  break;
         end
   end
    cpyvk=vk;
    vyas=vyas+1
end

%%%%Main Code for final image%%%%%%
%  ci = [temp_ci(:,:,1)/100,(temp_ci(:,:,2)+87)/186,(temp_ci(:,:,3)+108)/203].';
for i=1:size(vk,2)
    l=(vk{i}(1,1)*100);
    a=(vk{i}(2,1)*186-87);
    b=(vk{i}(3,1)*203-108);
    rgb=255*lab2rgb([l a b]);
    m=floor((i-1)/wo)+1;
    n=mod(i-1,wo)+1;
    finalImage(m,n,1)=rgb(1,1);finalImage(m,n,2)=rgb(1,2);finalImage(m,n,3)=rgb(1,3);
end
imshow(uint8(in));figure
imshow(uint8(finalImage));

               
                        
   