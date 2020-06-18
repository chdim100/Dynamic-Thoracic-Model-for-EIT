%%%%load the measurements%%%%%
N=16;
clear rawimagedata RLimagedata RRimagedata Himagedata MRLimagedata MRRimagedata MHimagedata
clear imageCN
%%%%skip-m current pattern
skipcurr=2;
%%%%skip-n voltage pattern
skipvolt=0;
if exist('Measurements')~=1
    error('please import/load the Measurements!')
end
%%%%%set directory to the EIDORS library:
%%%%%%% path2eidors='C:\.....\eidors-v3.9-ng\eidors\';

if exist('path2eidors')~=1
    error('please set a directory to the EIDORS library') 
else
    run ([path2eidors 'startup.m'])
end

%%%% add some Gaussian noise to the measurements
noise = randn(size(Measurements))*norm(Measurements)/10^(90/20);
Measurementsn=Measurements+noise;
clc
%%%%calculate their SNR, manually redefine noise until desired SNR is
%%%%reached
fprintf('SNR:\n')
20*log10(norm(Measurements)/norm(noise))
if N==16
    %heuristically selected values for the 50-60dB SNR levels
    lambda=0.25;
elseif N==32
    %%%% the large difference in hyperparameters has to do with the Jacobian's norms, 
    %%%% not the  problems' ill-conditioning 
    lambda=0.01;
end
for i=1:length(Measurements)/N^2-1
    %%%%Differential iterative Gauss-Newton Reconstruction Algorithm
    %%%%for more choices, see function 
[imageCN(i)]=Gauss_Newton_Reconstruction(N,skipcurr,skipvolt,lambda,Measurementsn(i*N^2+1:(i+1)*N^2)',0,2,2,2,1,Measurementsn(1:N^2)');
end
figure
%%%%show the temporal conductivity behavior
show_slices(imageCN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Define and plot ROIs 
model=imageCN(1);
elements=model.fwd_model.elems;
L=length(elements);
elementcentre=zeros(L,2);
for element=1:L
    nodes=[model.fwd_model.elems(element,1) model.fwd_model.elems(element,2) model.fwd_model.elems(element,3)];
    elementcentre(element,1)=(model.fwd_model.nodes(nodes(1),1)+model.fwd_model.nodes(nodes(2),1)+model.fwd_model.nodes(nodes(3),1))/3;
    elementcentre(element,2)=(model.fwd_model.nodes(nodes(1),2)+model.fwd_model.nodes(nodes(2),2)+model.fwd_model.nodes(nodes(3),2))/3;
end

figure
clf
H1=show_fem(imageCN(9));
set(H1, 'edgecolor', 'none');
axis off
hold on
a=35; % horizontal radius
b=60; % vertical radius
x0=-70; % x0,y0 ellipse centre coordinates
y0=0;
t=-pi:0.01:pi;
x=x0+a*cos(t);
y=y0+b*sin(t);
plot(x,y,'r','LineWidth',3)
axis square
hold on
a1=35; % horizontal radius
b1=60; % vertical radius
x01=70; % x0,y0 ellipse centre coordinates
y01=0;
t=-pi:0.01:pi;
x1=x01+a1*cos(t);
y1=y01+b1*sin(t);
plot(x1,y1,'r','LineWidth',3)
axis square
hold on
a2=35; % horizontal radius
b2=35; % vertical radius
x02=0; % x0,y0 ellipse centre coordinates
y02=-35;
t=-pi:0.01:pi;
x2=x02+a2*cos(t);
y2=y02+b2*sin(t);
plot(x2,y2,'r','LineWidth',3)
axis square

XX=elementcentre(:,1); YY=elementcentre(:,2);
RLROI=find((XX-x0).^2/a^2+(YY-y0).^2/b^2<1);
RRROI=find((XX-x01).^2/a1^2+(YY-y01).^2/b1^2<1);
HROI=find((XX-x02).^2/a2^2+(YY-y02).^2/b2^2<1);

for frame=1:length(imageCN)
    rawimagedata(frame,:)=imageCN(frame).elem_data;
    RLimagedata(frame,:)=rawimagedata(frame,RLROI);
    RRimagedata(frame,:)=rawimagedata(frame,RRROI);
    Himagedata(frame,:)=rawimagedata(frame,HROI);
    MRLimagedata(frame)=mean(RLimagedata(frame,:));
    MRRimagedata(frame)=mean(RRimagedata(frame,:));
    MHimagedata(frame)=mean(Himagedata(frame,:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%plot mean ROI element values 
figure
plot(MRLimagedata,'--*b','LineWidth',2)
hold on
plot(MRRimagedata,'--*r','LineWidth',2)
hold on
plot(MHimagedata,'--*','color',[255,215,0]/256,'LineWidth',2)
hold on