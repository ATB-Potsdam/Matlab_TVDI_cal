clear 
close all
clc

%% DN to Ref
Mrho=0.00002;
Arho=-0.1;
Prompt={'location of RED band',...
                'location of NIR  band', 'SunElvation degree'};
        Title='Red, NIR bands, sun elevation degree';
          DefaultValues={'b419323.tif',...
              'b519323.tif','56.73892978'};
           PARAMS=inputdlg(Prompt,Title,[1 200; 1 200; 1 50],DefaultValues);
red=imread(PARAMS{1});
nir=imread(PARAMS{2});
tetha=deg2rad(str2num(PARAMS{3})); %#ok

[h,v]=size(red);
for i=1:h
    for j=1:v
        if red(i,j)<0
            red(i,j)=0;
            nir(i,j)=0;
        end
    end
end
red=double(red);
nir=double(nir);
refred=((red.*Mrho)+Arho)./sin(tetha);
refnir=(nir.*Mrho)+Arho ./sin(tetha);
[hh,vv]=size(refred);
for ii=1:hh
    for jj=1:vv
        if refred(ii,jj)== (-0.1/sin(tetha))
            refred(ii,jj)=0;
            refnir(ii,jj)=0;
        end
    end
end
pause(0.1);
%% DN to Radiance
Prompt={'location of thermal band#10', 'location of thermal band#11'};
        Title='Thermal bands';
          DefaultValues={'b1019323.tif',...
              'b1119323.tif'};
           PARAMS=inputdlg(Prompt,Title,[1 100; 1 100],DefaultValues);

Band10=double(imread(PARAMS{1}));
Band11=double(imread(PARAMS{2}));
ML= double(3.3420E-04);
AL=0.1;
rad10=(Band10.*ML)+AL;
rad11=(Band11.*ML)+AL;
%% NDVI
NDVI=(refnir-refred)./(refnir+refred);
imtiff=geotiffinfo(PARAMS{1});
[A, R] = geotiffread(PARAMS{1}); 
mmm=imtiff.GeoTIFFCodes;
RRR=mmm.PCS;
geotiffwrite('NDVI.tif',NDVI,R,'CoordRefSysCode',RRR)
[NDVI, RR] = geotiffread('NDVI.tif');
figure
mapshow(NDVI, RR, 'DisplayType', 'mesh');
title('NDVI')
colorbar
axis image off

%% SAVI
% Soil line
figure
scatter(nonzeros(refred),nonzeros(refnir),'k.');
xlabel('Red band Reflection (\mum) ')
ylabel('NIR band Reflection (\mum)')
dif=((0.95.*nonzeros(refnir))-(1.05 .* nonzeros(refred)));
m=[dif nonzeros(refred) nonzeros(refnir)];
[hhh,vvv]=size(m);
for iii=1:hhh
     for jjj=1:vvv
         if m(iii,2)> m(iii,3)
             m(iii,jjj)=0;
         end
     end
 end
 pause(0.1)
 for iiii=1:hhh
     for jjjj=1:vvv
         if m(iiii,1)>=0
             m(iiii,jjjj)=0;
         end
     end
 end
 s1=nonzeros(m(:,1));
 s2=nonzeros(m(:,2));
 s3=nonzeros(m(:,3));
 s=[s1 s2 s3];

hold on
scatter(s2,s3,'rx');
grid on
X = [ones(size(s2)) s2];
b = double(regress(s3,X));
YFIT = b(1) + b(2).* s2;
plot(s2,YFIT);
legend('Red vs NIR','Bare soil pixels','Soil line')
grid on
pause(0.1)
% SAVI calculation
WDVI= refnir-(refred .* b(2));
L= 1- 2*b(2).* NDVI .* WDVI;
SAVI=((refnir-refred)./(refnir+refred+L)).* (1+L);
geotiffwrite('SAVI.tif',SAVI,R,'CoordRefSysCode',RRR)
[SAVI, RR] = geotiffread('SAVI.tif');
figure
mapshow(SAVI, RR, 'DisplayType', 'mesh');
title('SAVI')
colorbar
axis image off
pause(0.1);
%% Pv
NDVI1=NDVI;
[aa,bb]=size(NDVI);
for i=1:aa
    for j=1:bb
        if NDVI(i,j)<0.2
            NDVI1(i,j)=0.2;
        elseif NDVI(i,j)>0.5
            NDVI1(i,j)=0.5;
        else
             NDVI1(i,j)=NDVI(i,j);
        end
    end
end

Pv=((NDVI1-0.2)./(0.3)).^2;

%% Surface Emissivity
F=0.55;
es=0.966; 
ev=0.973;
Ci=(1-es)*ev*F.*(1-Pv);
em=zeros(size(Pv));
YY=ev.*Pv;
YY2=es.*(1-Pv);


for qq=1:aa
    for pp=1:bb
        if NDVI(qq,pp)<0.2
            em(qq,pp)= 0.973+(0.047.*refred(qq,pp)); 
  
        elseif 0.2<=NDVI(qq,pp)<=0.5
            em(qq,pp)=YY(qq,pp)+YY2(qq,pp)+Ci(qq,pp);
            

      
        else
            em(qq,pp)= ev +Ci(qq,pp); 
           
        end
    end
end
pause(0.1)
for qq=1:aa
    for pp=1:bb
        if em(qq,pp)>1
            em(qq,pp)= 1;
                   
        end
    end
end
geotiffwrite('em.tif',em,R,'CoordRefSysCode',RRR)
[em, RR] = geotiffread('em.tif');
figure
mapshow(em, RR, 'DisplayType', 'mesh');
title('Surface emissivity')
colorbar
axis image off
pause(0.1)
%% Transmittance
Prompt={'relative humidity (%)',...
                'temparture near the surface (K)'};
        Title='RH, T0';
          DefaultValues={'53.1','298.15'};
           PARAMS=inputdlg(Prompt,Title,[1 20; 1 20],DefaultValues);
RH=(str2num(PARAMS{1}))/100; %#ok
if RH>1
     msgbox('RH should be between 0-100')
     return
end
if RH<0
     msgbox('RH should be between 0-100')
     return
end
  T0=str2num(PARAMS{2});%#ok
  w0=(17.27*(T0-273.15))/(273.3+(T0-273.15));
  w1=exp(w0)*RH*6.108;
  w2=w1*0.0981;
  w=w2+0.1697;

if (w>0.2) && (w<3)
    tow10=(-0.0164*(w^2))-(0.04203*w)+0.9715;
    tow11=(-0.01218*(w^2))-(0.0735*w)+0.9603;
elseif  (w>=3) && (w<6)
    tow10=(-0.00168*(w^2))-(0.1329*w)+1.127;
    tow11=(0.00918*(w^2))-(0.2137*w)+1.181;
else
    msgbox('w is out of range, recheck the inputs')
    return
end
%% linear regression for a10, a11, b10, b11 
K1b10=774.8853;
K2b10=1321.0789;
K1b11=480.8883;
K2b11=1201.1442;
T10=log(((K1b10)./rad10)+1);
BT10=(K2b10./T10);
[rf,rf2]=size(BT10);
for i=1:rf
    for j=1:rf2
        if BT10(i,j)> 320
            BT10(i,j)=0;
        end
    end
end
T11=log(((K1b11)./rad11)+1);
BT11=(K2b11./T11);
[rf,rf2]=size(BT11);
for i=1:rf
    for j=1:rf2
        if BT11(i,j)> 320
            BT11(i,j)=0;
        end
    end
end
c1=1.19104356E08;
c2=14387.7;
lambda10=10.9;
lambda11=12;
f10=((BT10).^2).*(lambda10/c2);
ff10=1-exp(-c2./((lambda10).*(BT10)));
fff10=f10.*ff10;
x=[ones(size(nonzeros(BT10))) nonzeros(BT10)];
bb10=regress(nonzeros(fff10),x);
b10=bb10(2);a10=bb10(1);
figure;
scatter(nonzeros(BT10),nonzeros(fff10));
hold on
Yfit10=(b10.* nonzeros(BT10))+a10;
plot(nonzeros(BT10),Yfit10)
grid on
R210=(corr(nonzeros(fff10),Yfit10)).^2;
f11=((BT11).^2).*(lambda11/c2);
ff11=1-exp(-c2./((lambda11).*(BT11)));
fff11=f11.*ff11;
x=[ones(size(nonzeros(BT11))) nonzeros(BT11)];
bb11=regress(nonzeros(fff11),x);
b11=bb11(2);a11=bb11(1);
figure;
scatter(nonzeros(BT11),nonzeros(fff11));
hold on
Yfit11=(b11.* nonzeros(BT11))+a11;
plot(nonzeros(BT11),Yfit11)
grid on
R211=(corr(nonzeros(fff11),Yfit11)).^2;
a10 %#ok
b10 %#ok
a11 %#ok
b11 %#ok
%% SW algorithm for LST
D10=(1-tow10)*(1+((1-em).*tow10));
D11=(1-tow11)*(1+((1-em).*tow11));
C10=em.*tow10;
C11=em.*tow11;
E0=(D11.*C10)-(D10.*C11);
E2=(D10.*(1-C11-D11))./E0;
E1=(D11.*(1-C10-D10))./E0;
A=D10./E0;
A2=A+(E2.*b11);
A1=1+A+(E1.*b10);
A0=(E1.*a10)-(E2.*a11);
LSTSW=A0+(A1.*BT10)-(A2.*BT11);
[rf,rf2]=size(LSTSW);
for i=1:rf
    for j=1:rf2
        if LSTSW(i,j)<270
            LSTSW(i,j)=nan;
        end
    end
end

geotiffwrite('LSTSW.tif',LSTSW,R,'CoordRefSysCode',RRR)
[LSTSW, RR] = geotiffread('LSTSW.tif');
figure
mapshow(LSTSW, RR, 'DisplayType', 'mesh');
title('Land surface temperature')
colorbar
axis image off
MeanLST=nanmean(nanmean(LSTSW));


%% TVDI
LST=LSTSW;
for i=1:size(NDVI,1)
    for j=1:size(NDVI,2)
    if NDVI(i,j)<=0
        NDVI(i,j)=nan;
        LST(i,j)=nan;
    elseif isnan(LST(i,j))>0
        NDVI(i,j)=nan;
    end
    end
end
NDVI2=NDVI;LST2=LST; % will be used later
LST=reshape(LST,[],1);
NDVI=reshape(NDVI,[],1);
NDVI(isnan(NDVI))=[];
LST(isnan(LST))=[];
figure
scatter(NDVI, LST,'k.')

L=linspace(min(NDVI),max(NDVI),11);

k1 = find(NDVI>L(1) & NDVI<= L(2)); LSTmax1=max(LST(k1));LSTmin1=min(LST(k1));
k2 = find(NDVI>L(2) & NDVI<= L(3)); LSTmax2=max(LST(k2));LSTmin2=min(LST(k2));
k3 = find(NDVI>L(3) & NDVI<= L(4)); LSTmax3=max(LST(k3));LSTmin3=min(LST(k3));
k4 = find(NDVI>L(4) & NDVI<= L(5)); LSTmax4=max(LST(k4));LSTmin4=min(LST(k4));
k5 = find(NDVI>L(5) & NDVI<= L(6)); LSTmax5=max(LST(k5));LSTmin5=min(LST(k5));
k6 = find(NDVI>L(6) & NDVI<= L(7)); LSTmax6=max(LST(k6));LSTmin6=min(LST(k6));
k7 = find(NDVI>L(7) & NDVI<= L(8)); LSTmax7=max(LST(k7));LSTmin7=min(LST(k7));
k8 = find(NDVI>L(8) & NDVI<= L(9)); LSTmax8=max(LST(k8));LSTmin8=min(LST(k8));
k9 = find(NDVI>L(9) & NDVI<= L(10)); LSTmax9=max(LST(k9));LSTmin9=min(LST(k9));
k10 = find(NDVI>L(10) & NDVI<= L(11)); LSTmax10=max(LST(k10));LSTmin10=min(LST(k10));

LSTmin=[LSTmin1 LSTmin2 LSTmin3 LSTmin4 LSTmin5 LSTmin6 LSTmin7 LSTmin8 LSTmin9 LSTmin10];
LSTmax=[LSTmax1 LSTmax2 LSTmax3 LSTmax4 LSTmax5 LSTmax6 LSTmax7 LSTmax8 LSTmax9 LSTmax10];

for i=1:10
    L(i)=(L(i)+ L(i+1))*0.5;
    
end
NDVIData=L(1:10)';
x=[ones(length(NDVIData),1) NDVIData];
y1=LSTmax';
y2=LSTmin';
bmax=regress(y1,x);
a1=bmax(1) %#ok 
b1=bmax(2) %#ok
bmin=regress(y2,x);
a2=bmin(1) %#ok 
b2=bmin(2) %#ok
LSTmaxLine=a1+(b1.*NDVI2);
LSTminLine=a2+(b2.*NDVI2);
hold on
plot(NDVI2,LSTmaxLine,'r')
plot(NDVI2,LSTminLine,'b')
grid on
TVDI=((LST2-LSTminLine)./(LSTmaxLine-LSTminLine));
[a,b]=size(TVDI);
for i=1:a
    for j=1:b
        if TVDI(i,j)<0
            TVDI(i,j)=0;

        elseif NDVI2(i,j)<=0.2
            TVDI(i,j)=1;
        end
    end
end

geotiffwrite('TVDI.tif',TVDI,R,'CoordRefSysCode',RRR)
[TVDI, RR] = geotiffread('TVDI.tif');
figure
mapshow(TVDI, RR, 'DisplayType', 'mesh');
title('TVDI')
colorbar
axis image off
% 

