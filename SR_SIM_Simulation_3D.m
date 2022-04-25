clc;
clear all;
close all;

%% Set up image parameters
disp('Set up parameters...');
w1 = 256;
w2 = 128;
wo1 = w1/2;
wo2 = w2/2;
x = linspace(0,w1-1,w1);
z = linspace(0,w2-1,w2);

%% Spatial frequency(0D,120D,240D)
fx = 0.42;
fy = 0;
fz = sqrt(0.67^2-fx^2-fy^2);

f=sqrt(fx^2+fy^2+fz^2);

fx1 = 0.42*cos(pi/3);
fy1 = 0.42*sin(pi/3);
fz1 = sqrt(0.67^2-fx^2-fy^2);

f1=sqrt(fx1^2+fy1^2+fz1^2);

fx2 = 0.42*cos(pi*2/3);
fy2 = 0.42*sin(pi*2/3);
fz2 = sqrt(0.67^2-fx2^2-fy2^2);

f2=sqrt(fx2^2+fy2^2+fz2^2);

[Y,X,Z]=meshgrid(x,x,z);
% 


% figure,imagesc(reshape(Ip1o(:,129,:),256,128));

% freq=fft2(reshape(Ip1o(:,129,:),256,128));
% figure,imagesc((abs(fftshift(freq))));

%% Five phase images
disp('Generating 15 simulation raw images...');
phase1 = 0;
phase2 = 2*pi/5;
phase3 = 4*pi/5;
phase4 = 6*pi/5;
phase5 = 8*pi/5;

Ip1=0.5+cos(2*fx*X+2*fy*Y+2*phase1);
Ip2=0.5+cos(2*fx*X+2*fy*Y+2*phase2);
Ip3=0.5+cos(2*fx*X+2*fy*Y+2*phase3);
Ip4=0.5+cos(2*fx*X+2*fy*Y+2*phase4);
Ip5=0.5+cos(2*fx*X+2*fy*Y+2*phase5);

Ip1o=2*cos(fx*X+fy*Y+phase1);
Ip2o=2*cos(fx*X+fy*Y+phase2);
Ip3o=2*cos(fx*X+fy*Y+phase3);
Ip4o=2*cos(fx*X+fy*Y+phase4);
Ip5o=2*cos(fx*X+fy*Y+phase5);

Ipa1=1+cos(2*fx1*X+2*fy1*Y+2*phase1);
Ipa2=1+cos(2*fx1*X+2*fy1*Y+2*phase2);
Ipa3=1+cos(2*fx1*X+2*fy1*Y+2*phase3);
Ipa4=1+cos(2*fx1*X+2*fy1*Y+2*phase4);
Ipa5=1+cos(2*fx1*X+2*fy1*Y+2*phase5);

Ipa1o=cos(fx1*X+fy1*Y+phase1);
Ipa2o=cos(fx1*X+fy1*Y+phase2);
Ipa3o=cos(fx1*X+fy1*Y+phase3);
Ipa4o=cos(fx1*X+fy1*Y+phase4);
Ipa5o=cos(fx1*X+fy1*Y+phase5);

Ipb1=1+cos(2*fx2*X+2*fy2*Y+2*phase1);
Ipb2=1+cos(2*fx2*X+2*fy2*Y+2*phase2);
Ipb3=1+cos(2*fx2*X+2*fy2*Y+2*phase3);
Ipb4=1+cos(2*fx2*X+2*fy2*Y+2*phase4);
Ipb5=1+cos(2*fx2*X+2*fy2*Y+2*phase5);

Ipb1o=cos(fx2*X+fy2*Y+phase1);
Ipb2o=cos(fx2*X+fy2*Y+phase2);
Ipb3o=cos(fx2*X+fy2*Y+phase3);
Ipb4o=cos(fx2*X+fy2*Y+phase4);
Ipb5o=cos(fx2*X+fy2*Y+phase5);

%% Generate PSF/OTF
disp('Generating PSF/OTF...');
for i=1:128
    psf(:,:,i)=double(imread('psf3d2.tif',i));
     obj(:,:,i)=double(imread('ground-truth-tube.tif',i));
end

otf=fftshift(fftn(psf));
otf=otf/max(max(max(abs(otf))));

dist=2*fx/pi/2/(x(2)-x(1))*(w1-1);
distz=(f-fz)/pi/2/(z(2)-z(1))*(w2-1);

otfk = fftn( (ifftn(otf)).*exp( +1i.*2*pi*( distz/w2.*(Z-wo2)) ))+fftn( (ifftn(otf)).*exp( -1i.*2*pi*( distz/w2.*(Z-wo2)) ));

% ima0=convn(obj,psf,'same');
ima1=fftshift(ifftn(fftn(obj.*Ip1).*fftshift(otf)+fftn(obj.*Ip1o).*fftshift(otfk)));
ima2=fftshift(ifftn(fftn(obj.*Ip2).*fftshift(otf)+fftn(obj.*Ip2o).*fftshift(otfk)));
ima3=fftshift(ifftn(fftn(obj.*Ip3).*fftshift(otf)+fftn(obj.*Ip3o).*fftshift(otfk)));
ima4=fftshift(ifftn(fftn(obj.*Ip4).*fftshift(otf)+fftn(obj.*Ip4o).*fftshift(otfk)));
ima5=fftshift(ifftn(fftn(obj.*Ip5).*fftshift(otf)+fftn(obj.*Ip5o).*fftshift(otfk)));

ima1a=fftshift(ifftn(fftn(obj.*Ipa1).*fftshift(otf)+fftn(obj.*Ipa1o).*fftshift(otfk)));
ima2a=fftshift(ifftn(fftn(obj.*Ipa2).*fftshift(otf)+fftn(obj.*Ipa2o).*fftshift(otfk)));
ima3a=fftshift(ifftn(fftn(obj.*Ipa3).*fftshift(otf)+fftn(obj.*Ipa3o).*fftshift(otfk)));
ima4a=fftshift(ifftn(fftn(obj.*Ipa4).*fftshift(otf)+fftn(obj.*Ipa4o).*fftshift(otfk)));
ima5a=fftshift(ifftn(fftn(obj.*Ipa5).*fftshift(otf)+fftn(obj.*Ipa5o).*fftshift(otfk)));

ima1b=fftshift(ifftn(fftn(obj.*Ipb1).*fftshift(otf)+fftn(obj.*Ipb1o).*fftshift(otfk)));
ima2b=fftshift(ifftn(fftn(obj.*Ipb2).*fftshift(otf)+fftn(obj.*Ipb2o).*fftshift(otfk)));
ima3b=fftshift(ifftn(fftn(obj.*Ipb3).*fftshift(otf)+fftn(obj.*Ipb3o).*fftshift(otfk)));
ima4b=fftshift(ifftn(fftn(obj.*Ipb4).*fftshift(otf)+fftn(obj.*Ipb4o).*fftshift(otfk)));
ima5b=fftshift(ifftn(fftn(obj.*Ipb5).*fftshift(otf)+fftn(obj.*Ipb5o).*fftshift(otfk)));

FS1aT=fftshift(fftn(ima1));
FS2aT=fftshift(fftn(ima2));
FS3aT=fftshift(fftn(ima3));
FS4aT=fftshift(fftn(ima4));
FS5aT=fftshift(fftn(ima5));

FS1aTa=fftshift(fftn(ima1a));
FS2aTa=fftshift(fftn(ima2a));
FS3aTa=fftshift(fftn(ima3a));
FS4aTa=fftshift(fftn(ima4a));
FS5aTa=fftshift(fftn(ima5a));

FS1aTb=fftshift(fftn(ima1b));
FS2aTb=fftshift(fftn(ima2b));
FS3aTb=fftshift(fftn(ima3b));
FS4aTb=fftshift(fftn(ima4b));
FS5aTb=fftshift(fftn(ima5b));

%% Start SIM Reconstruction
disp('Multiplying the SIM matrix...');

ModFac=1;
M = [exp(-1i*2*phase1) exp(-1i*phase1) 1 exp(1i*phase1) exp(1i*2*phase1);
     exp(-1i*2*phase2) exp(-1i*phase2) 1 exp(1i*phase2) exp(1i*2*phase2);
     exp(-1i*2*phase3) exp(-1i*phase3) 1 exp(1i*phase3) exp(1i*2*phase3);
     exp(-1i*2*phase4) exp(-1i*phase4) 1 exp(1i*phase4) exp(1i*2*phase4);
     exp(-1i*2*phase5) exp(-1i*phase5) 1 exp(1i*phase5) exp(1i*2*phase5);];
 
Minv = inv(M);
FS1 = Minv(1,1)*FS1aT + Minv(1,2)*FS2aT + Minv(1,3)*FS3aT+ Minv(1,4)*FS4aT + Minv(1,5)*FS5aT;
FS2 = Minv(2,1)*FS1aT + Minv(2,2)*FS2aT + Minv(2,3)*FS3aT+ Minv(2,4)*FS4aT + Minv(2,5)*FS5aT;
FS3 = Minv(3,1)*FS1aT + Minv(3,2)*FS2aT + Minv(3,3)*FS3aT+ Minv(3,4)*FS4aT + Minv(3,5)*FS5aT;
FS4 = Minv(4,1)*FS1aT + Minv(4,2)*FS2aT + Minv(4,3)*FS3aT+ Minv(4,4)*FS4aT + Minv(4,5)*FS5aT;
FS5 = Minv(5,1)*FS1aT + Minv(5,2)*FS2aT + Minv(5,3)*FS3aT+ Minv(5,4)*FS4aT + Minv(5,5)*FS5aT;

FS1a = Minv(1,1)*FS1aTa + Minv(1,2)*FS2aTa + Minv(1,3)*FS3aTa+ Minv(1,4)*FS4aTa + Minv(1,5)*FS5aTa;
FS2a = Minv(2,1)*FS1aTa + Minv(2,2)*FS2aTa + Minv(2,3)*FS3aTa+ Minv(2,4)*FS4aTa + Minv(2,5)*FS5aTa;
FS3a = Minv(3,1)*FS1aTa + Minv(3,2)*FS2aTa + Minv(3,3)*FS3aTa+ Minv(3,4)*FS4aTa + Minv(3,5)*FS5aTa;
FS4a = Minv(4,1)*FS1aTa + Minv(4,2)*FS2aTa + Minv(4,3)*FS3aTa+ Minv(4,4)*FS4aTa + Minv(4,5)*FS5aTa;
FS5a = Minv(5,1)*FS1aTa + Minv(5,2)*FS2aTa + Minv(5,3)*FS3aTa+ Minv(5,4)*FS4aTa + Minv(5,5)*FS5aTa;

FS1b = Minv(1,1)*FS1aTb + Minv(1,2)*FS2aTb + Minv(1,3)*FS3aTb+ Minv(1,4)*FS4aTb + Minv(1,5)*FS5aTb;
FS2b = Minv(2,1)*FS1aTb + Minv(2,2)*FS2aTb + Minv(2,3)*FS3aTb+ Minv(2,4)*FS4aTb + Minv(2,5)*FS5aTb;
FS3b = Minv(3,1)*FS1aTb + Minv(3,2)*FS2aTb + Minv(3,3)*FS3aTb+ Minv(3,4)*FS4aTb + Minv(3,5)*FS5aTb;
FS4b = Minv(4,1)*FS1aTb + Minv(4,2)*FS2aTb + Minv(4,3)*FS3aTb+ Minv(4,4)*FS4aTb + Minv(4,5)*FS5aTb;
FS5b = Minv(5,1)*FS1aTb + Minv(5,2)*FS2aTb + Minv(5,3)*FS3aTb+ Minv(5,4)*FS4aTb + Minv(5,5)*FS5aTb;

clear FS1aT FS2aT FS3aT FS4aT FS5aT FS1aTa FS2aTa FS3aTa FS4aTa FS5aTa FS1aTb FS2aTb FS3aTb FS4aTb FS5aTb;

%% Shift frequency components
disp('Calculating and shifting frequency components...');
%  dist=64.5;
kA=[dist,0,0];
kB=[dist/2,0,0];

kAa=[dist*cos(pi/3),dist*sin(pi/3),0];
kBa=[dist/2*cos(pi/3),dist/2*sin(pi/3),0];

kAb=[dist*cos(pi*2/3),dist*sin(pi*2/3),0];
kBb=[dist/2*cos(pi*2/3),dist/2*sin(pi*2/3),0];

Shift_F2 = fftn( (ifftn(FS2)).*exp( +1i.*2*pi*(kB(1)/w1.*(X-wo1) + kB(2)/w1.*(Y-wo1)+ kB(3)/w2.*(Z-wo2)) ));
Shift_F4 = fftn( (ifftn(FS4)).*exp( -1i.*2*pi*(kB(1)/w1.*(X-wo1) + kB(2)/w1.*(Y-wo1)+ kB(3)/w2.*(Z-wo2)) ));

Shift_F1 = fftn( (ifftn(FS1)).*exp( +1i.*2*pi*(kA(1)/w1.*(X-wo1) + kA(2)/w1.*(Y-wo1)+ kA(3)/w2.*(Z-wo2)) ));
Shift_F5 = fftn( (ifftn(FS5)).*exp( -1i.*2*pi*(kA(1)/w1.*(X-wo1) + kA(2)/w1.*(Y-wo1)+ kA(3)/w2.*(Z-wo2)) ));

Shift_F2a = fftn( (ifftn(FS2a)).*exp( +1i.*2*pi*(kBa(1)/w1.*(X-wo1) + kBa(2)/w1.*(Y-wo1)+ kBa(3)/w2.*(Z-wo2)) ));
Shift_F4a = fftn( (ifftn(FS4a)).*exp( -1i.*2*pi*(kBa(1)/w1.*(X-wo1) + kBa(2)/w1.*(Y-wo1)+ kBa(3)/w2.*(Z-wo2)) ));

Shift_F1a = fftn( (ifftn(FS1a)).*exp( +1i.*2*pi*(kAa(1)/w1.*(X-wo1) + kAa(2)/w1.*(Y-wo1)+ kAa(3)/w2.*(Z-wo2)) ));
Shift_F5a = fftn( (ifftn(FS5a)).*exp( -1i.*2*pi*(kAa(1)/w1.*(X-wo1) + kAa(2)/w1.*(Y-wo1)+ kAa(3)/w2.*(Z-wo2)) ));

Shift_F2b = fftn( (ifftn(FS2b)).*exp( +1i.*2*pi*(kBb(1)/w1.*(X-wo1) + kBb(2)/w1.*(Y-wo1)+ kBb(3)/w2.*(Z-wo2)) ));
Shift_F4b = fftn( (ifftn(FS4b)).*exp( -1i.*2*pi*(kBb(1)/w1.*(X-wo1) + kBb(2)/w1.*(Y-wo1)+ kBb(3)/w2.*(Z-wo2)) ));

Shift_F1b = fftn( (ifftn(FS1b)).*exp( +1i.*2*pi*(kAb(1)/w1.*(X-wo1) + kAb(2)/w1.*(Y-wo1)+ kAb(3)/w2.*(Z-wo2)) ));
Shift_F5b = fftn( (ifftn(FS5b)).*exp( -1i.*2*pi*(kAb(1)/w1.*(X-wo1) + kAb(2)/w1.*(Y-wo1)+ kAb(3)/w2.*(Z-wo2)) ));

clear FS2 FS4 FS1 FS5 FS2a FS4a FS1a FS5a FS2b FS4b FS1b FS5b; 

%% Wiener filtering for deconvolution
disp('Wiener filtering the shifted frequency components...');
% otf2 = fftn( (ifftn(otf)).*exp( +1i.*2*pi*(kB(1)/w1.*(X-wo1) + kB(2)/w1.*(Y-wo1)+ kB(3)/w2.*(Z-wo2)) ));
% otf4 = fftn( (ifftn(otf)).*exp( -1i.*2*pi*(kB(1)/w1.*(X-wo1) + kB(2)/w1.*(Y-wo1)+ kB(3)/w2.*(Z-wo2)) ));

otf1 = fftn( (ifftn(otf)).*exp( +1i.*2*pi*(kA(1)/w1.*(X-wo1) + kA(2)/w1.*(Y-wo1)+ kA(3)/w2.*(Z-wo2)) ));
otf5 = fftn( (ifftn(otf)).*exp( -1i.*2*pi*(kA(1)/w1.*(X-wo1) + kA(2)/w1.*(Y-wo1)+ kA(3)/w2.*(Z-wo2)) ));

otf2 = fftn( (ifftn(otfk)).*exp( +1i.*2*pi*(kB(1)/w1.*(X-wo1) + kB(2)/w1.*(Y-wo1)+ kB(3)/w2.*(Z-wo2)) ));
otf4 = fftn( (ifftn(otfk)).*exp( -1i.*2*pi*(kB(1)/w1.*(X-wo1) + kB(2)/w1.*(Y-wo1)+ kB(3)/w2.*(Z-wo2)) ));

otf1a = fftn( (ifftn(otf)).*exp( +1i.*2*pi*(kAa(1)/w1.*(X-wo1) + kAa(2)/w1.*(Y-wo1)+ kAa(3)/w2.*(Z-wo2)) ));
otf5a = fftn( (ifftn(otf)).*exp( -1i.*2*pi*(kAa(1)/w1.*(X-wo1) + kAa(2)/w1.*(Y-wo1)+ kAa(3)/w2.*(Z-wo2)) ));

otf2a = fftn( (ifftn(otfk)).*exp( +1i.*2*pi*(kBa(1)/w1.*(X-wo1) + kBa(2)/w1.*(Y-wo1)+ kBa(3)/w2.*(Z-wo2)) ));
otf4a = fftn( (ifftn(otfk)).*exp( -1i.*2*pi*(kBa(1)/w1.*(X-wo1) + kBa(2)/w1.*(Y-wo1)+ kBa(3)/w2.*(Z-wo2)) ));

otf1b = fftn( (ifftn(otf)).*exp( +1i.*2*pi*(kAb(1)/w1.*(X-wo1) + kAb(2)/w1.*(Y-wo1)+ kAb(3)/w2.*(Z-wo2)) ));
otf5b = fftn( (ifftn(otf)).*exp( -1i.*2*pi*(kAb(1)/w1.*(X-wo1) + kAb(2)/w1.*(Y-wo1)+ kAb(3)/w2.*(Z-wo2)) ));

otf2b = fftn( (ifftn(otfk)).*exp( +1i.*2*pi*(kBb(1)/w1.*(X-wo1) + kBb(2)/w1.*(Y-wo1)+ kBb(3)/w2.*(Z-wo2)) ));
otf4b = fftn( (ifftn(otfk)).*exp( -1i.*2*pi*(kBb(1)/w1.*(X-wo1) + kBb(2)/w1.*(Y-wo1)+ kBb(3)/w2.*(Z-wo2)) ));

sigma=0.2;
FS30=FS3./(abs(otf).^2+sigma);
Shift_F10=Shift_F1./(abs(otf1).^2+sigma);
Shift_F50=Shift_F5./(abs(otf5).^2+sigma);
Shift_F20=Shift_F2./(abs(otf2).^2+sigma);
Shift_F40=Shift_F4./(abs(otf4).^2+sigma);

FS30a=FS3a./(abs(otf).^2+sigma);
Shift_F10a=Shift_F1a./(abs(otf1a).^2+sigma);
Shift_F50a=Shift_F5a./(abs(otf5a).^2+sigma);
Shift_F20a=Shift_F2a./(abs(otf2a).^2+sigma);
Shift_F40a=Shift_F4a./(abs(otf4a).^2+sigma);

FS30b=FS3b./(abs(otf).^2+sigma);
Shift_F10b=Shift_F1b./(abs(otf1b).^2+sigma);
Shift_F50b=Shift_F5b./(abs(otf5b).^2+sigma);
Shift_F20b=Shift_F2b./(abs(otf2b).^2+sigma);
Shift_F40b=Shift_F4b./(abs(otf4b).^2+sigma);

%% Generate final SIM image
disp('Generating final image and save into files...');

FSum=FS30+Shift_F10+Shift_F50+Shift_F20+Shift_F40+FS30a+Shift_F10a+Shift_F50a+Shift_F20a+Shift_F40a+...
    FS30b+Shift_F10b+Shift_F50b+Shift_F20b+Shift_F40b;
% FSum=FS30+Shift_F10+Shift_F50+Shift_F20+Shift_F40;



DSim=ifftn(fftshift(FSum));
WF=ima1+ima2+ima3+ima4+ima5+ima1a+ima2a+ima3a+ima4a+ima5a+ima1b+ima2b+ima3b+ima4b+ima5b;

DSim=DSim/max(max(max(DSim)))*65536;
DSim=uint16(DSim);

WF=WF/max(max(max(WF)))*65536;
WF=uint16(WF);

for i=1:128
    if i==1
      imwrite(DSim(:,:,i),'SR-SIM-tube.tif');
      imwrite(WF(:,:,i),'Widefield-tube.tif');
    else
         imwrite(DSim(:,:,i),'SR-SIM-tube.tif','tif','WriteMode', 'append');
         imwrite(WF(:,:,i),'Widefield-tube.tif','tif','WriteMode', 'append');
    end
end

