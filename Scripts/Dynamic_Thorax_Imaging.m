function [M,Gr,Patient_labeled_data]=Dynamic_Thorax_Imaging(Measurement_params,Model_params,Zelectrodes,randomstart,randomcond,paths)
%%%Implementation of the Dynamic Thoracic Model for EIT
%%%Creates a .lua script, to be executed in FEMM
%%%Automatically pauses, until FEMM performs the EM simulation ordered for
%%%the .lua
%%%Collects the EIT potential measurements

%Inputs: 
%Measurement_params: parameters that define the measuring properties:
%-N: the number of electrodes (single layer)
%-currentskip: skip-m current protocol
%-voltageskip: skip-n voltage measurement protocol
%-frequency: injected ac current frequency (in Hz)
%-total_frames: total image frames to be simulated
%-fps: measurement frames per second
%-Current: Injected ac current's Amplitude

%Model_params: parameters that define the measuring properties:
%1. end_Admittances: includes the inspiration/ expiration-end
%conductivities and permittivities of each tissue
    %-deflated lung conductivity (S/m)  end_Admittances(1)
    %-inflated lung conductivity (S/m)  end_Admittances(2)
    %-blood cycle-related conductivity variance in lungs (S/m) end_Admittances(3)
    %-deflated lung permittivity (F/m)  end_Admittances(4)
    %-inflated lung permittivity (F/m)  end_Admittances(5)
    %-blood cycle-related permittivity variance in lungs (F/m)  end_Admittances(6)
    %-base heart (champers, aorta, no-myocardium) conductivity (S/m) end_Admittances(7)
    %-base heart (champers, aorta, no-myocardium) conductivity variation
    %(S/m) end_Admittances(8)
    %-base heart (champers, aorta, no-myocardium) conductivity (S/m) end_Admittances(9)
    %-base heart (champers, aorta, no-myocardium) permittivity variation 
    %(F/m) end_Admittances(10)
%2. breath_time (seconds), Initial Breathing Time
%3. delec (cm), the electrodes width
%4. bpm, Initial Heart Rate
%5. Collapse_LL, Left Lung Collapsion state (1-4)
%6. Collapse_RL, Right Lung Collapsion state (1-4)

%Zelectrodes: A 2XN matrix. The first raw contains each electrode's contact
%relative conductivity. The second raw contains each electrode's contact
%relative permittivity. 

%random_start: 0 if simulation starts from pulmonary and cardiac states 1
%(full exhalation and start of heartbeat). 1 if simulation starts from
%random states at both cycles. 

%randomcond: if 1: breathing times and HR are randomly changed every cycle and 3
%cycles respectively. Else if 0: initial breathing time is 3 seconds,
%silent space is 0.3seconds and next breath lasts +20%. Initial heart beat
%rate is 75bpm and increases by 3% every 3 beats. 

%paths: A cell array. paths{1} denotes the directory to the Bio-lungs
%folder ('C:\.....\Bio-lungs\')
%paths{2} denotes the directory to the Eidors library 
%('C:\.....\eidors-v3.9-ng\eidors\')

%%%%%%%%%%%%

%Outputs: 
%M: The final, noise-free potential measurements. Automatically saved in 
%C:\....\Bio-lungs\Dynamic_Thorax_Model\Cardio-Pulmonary

%Gr: Simulated admittances graph. A matrix called Input_set_date() is automatically
%saved in C:\....\Bio-lungs\Dynamic_Thorax_Model\Cardio-Pulmonary



%%%% unpack Measurement_params

N=Measurement_params(1); currentskip=Measurement_params(2);
voltageskip=Measurement_params(3);
frequency=Measurement_params(4);
total_frames=Measurement_params(5);
fps=Measurement_params(6);
Current=Measurement_params(7);

%%%% unpack Model_params
end_Admittances=Model_params(1:10);
breath_time=Model_params(11);
delec=Model_params(12);
bpm=Model_params(13);
Collapse_LL=Model_params(14);
Collapse_RL=Model_params(15);

%%%% unpack Paths
path=paths{1};
path2EIDORS=paths{2};
path2save=paths{1};

%start EIDORS
run ([path2EIDORS 'startup.m'])

[initial_times,initial_total_positions,Collapsion_areas,electrodes,Tissues,Ventilation,Circulation,Graph,Silences,r,k]=...
    set_initial_parameters(N,fps,total_frames,frequency,breath_time,delec,bpm,Collapse_LL,Collapse_RL,end_Admittances,randomstart,randomcond);
clc
%% lua
%%%
[fileID,d,folds]=open_script(frequency,path);
set_materials_in_lua(fileID,Current,Tissues,Zelectrodes);
[Graph,timesteps,total_time_est,History_of_Tissues]=...
    frame_loop(N,currentskip,initial_times,initial_total_positions,...
    Collapsion_areas,electrodes,Silences,fileID,Graph, Circulation,...
    Ventilation, Tissues,d,r,k,frequency,end_Admittances,folds,randomcond);
[letterplus,Gr]=configure_graph(Graph,total_frames,fps,total_time_est,path2save);
M=configure_measurements(N,letterplus,currentskip,voltageskip,timesteps,d,path2save,path);
%imageC=thorax_GN_v3(32,4,1,0.02,Measurements',0,2,2,2,1,[]);
Patient_labeled_data=complete_patient_labeled_data(Graph,Collapsion_areas,History_of_Tissues);
end

function [initial_times,initial_total_positions,Collapsion_areas,electrodes,Tissues,Ventilation,Circulation,Graph,Silences,r,k]=...
    set_initial_parameters(N,fps,total_frames,f,breath_time,delec,bpm,Collapse_LL,Collapse_RL,end_Admittances,randomstart,randomcond)
%defines the initial parameters (time, step, collapsions,electrodes,tissue_configuration,all possible positions, initial cycle states)
initial_times=define_times(N,fps,total_frames,breath_time,bpm,randomcond);
initial_total_positions=define_total_positions(N,1,initial_times.tim,initial_times.tbr); %Dy=1
[Collapsion_areas.collapsions,Collapsion_areas.cross_points]=define_collapse_areas(Collapse_LL,Collapse_RL,randomcond);
initial_times.pulse_dur=(2*(initial_total_positions.total_positions-1)+1)*(initial_times.tbr/initial_times.Tpulse)^(-1);
if randomcond
    k=0.95 + (1.05-0.95)*randn(1,1);
else
    k=1;
end
[electrodes.electrode_average,electrodes.electrode_possible,...
    electrodes.delec]=electrodes_preset(N,1,initial_total_positions.dy,delec,k); %Dy=1
if randomcond
    r=0.95 + (1.05-0.95)*randn(1,1);
else
    r=1;
end
[Tissues.Left_lung,Tissues.Right_lung,Tissues.Heart,Tissues.Muscles,Tissues.Bones,r,k]=...
    define_tissues(r,k,initial_times,initial_total_positions,0,f,end_Admittances,randomcond);
%returns the initial states of breathing and cardiac cycles
%pulmocycle and cardio_cycle indicate the positions of lungs and heart
%pulmcounter and cardcounter indicate whether they are increasing or
%decreasing
%plot(Tissues.Bones.Bone3.shape.full(:,1),Tissues.Bones.Bone3.shape.full(:,2))
%hold on
%plot(electrodes.electrode_average.x(1,:),electrodes.electrode_average.y(1,:),'o')
[Ventilation.pulmocycle,Circulation.cardio_cycle,Ventilation.pulmcounter,Circulation.cardcounter,Circulation.pulse_dur]=...
    randomize(initial_times,initial_total_positions,randomstart); %randomstart==1
Graph=initialize_graph();
Silences=initialize_silences();
end

function initial_times=define_times(N,fps,total_frames,breath_time,bpm,randomcond)
%sets simulation times, breathing and blood cycle periods (all time
%parameters)
initial_times.tim=define_frame_time(N,fps);
initial_times.tbr=define_breath_time(breath_time,randomcond);
initial_times.total_time_est=estimate_total_simulation_time(initial_times.tim,total_frames);
initial_times.Tpulse=define_cardiac_cycle_time(bpm,randomcond);
end

function tim=define_frame_time(N,fps)
%estimates required time for a frame (in seconds)
if N==16||N==32
    tim=1/fps; %total time for a measurement cycle
elseif N==64
    tim=2; %when using 64 electrodes system is too slow to perform tdEIT
end
end

function tbr=define_breath_time(breath_time,randomcond)
%defined unitial breath period (randomly between 2 and 6 seconds)
%uniform distribution selection
if randomcond
    tbr=rand(1)*6;
    while tbr<2
        tbr=rand(1)*6; %re-defined breath period until passes 2 seconds
    end
else
    tbr=breath_time;
end
end

function total_time_est=estimate_total_simulation_time(tim,total_frames)
%estimates total simulation time, that will NOT change!
total_time_est=tim*total_frames;
end

function Tpulse=define_cardiac_cycle_time(bpm,randomcond)
%defines the initial cardiac frequency
if randomcond
    bpm=55+rand(1)*55;  %initial Bpm
else
    %bpm=75;
end
fpulse=bpm/60;
Tpulse=1/fpulse;
end

function initial_total_positions=define_total_positions(N,Dy,tim,tbr)
%initializes electrode movement step and total possible positions,
%according to the initialized breathing time
initial_total_positions.dy=2*Dy*tim/(N*tbr); %initializes electrode movement step
%total possible positions for tissues
initial_total_positions.total_positions=floor(Dy/initial_total_positions.dy+1);
end

function [collapsions,cross_points]=define_collapse_areas(Collapse_LL,Collapse_RL,randomcond)
if randomcond
    %randomize collapsing area for left lung
    Collapse_cases.left_lung=randi([1 4]);
    %randomize collapsing area for right lung
    Collapse_cases.right_lung=randi([1 4]);
    %degree of collaspion on left_lung
    collapsions.collapsion_left_lung=randi(5)-1;
    %degree of collaspion on right_lung
    collapsions.collapsion_right_lung=randi(4)-1;
else
    collapsions.collapsion_left_lung=Collapse_LL;
    collapsions.collapsion_right_lung=Collapse_RL;
%     collapsions.collapsion_left_lung=1;
%     collapsions.collapsion_right_lung=3;
end
switch collapsions.collapsion_left_lung
    case 1
        cross_points.left_lung=[4 14];
    case 2
        cross_points.left_lung=[3 14];
    case 3
        cross_points.left_lung=[3 15];
    case 4
        cross_points.left_lung=[2 16];
    otherwise
        cross_points.left_lung=[1 16];
end
switch collapsions.collapsion_right_lung
    case 1
        cross_points.right_lung=[2 12];
    case 2
        cross_points.right_lung=[2 14];
    case 3
        cross_points.right_lung=[1 14];
    otherwise
        cross_points.right_lung=[1 15];
end
end

function [electrode_average,electrode_possible,delec]=electrodes_preset(N,Dy,dy,delec,k)
if N==16
    F5=mk_common_model('d2T3',16);
    AAAx=zeros(1,14);
    BBBx=AAAx;
    AAAy=zeros(1,14);
    BBBy=AAAx;
    for iii=1:14
        AAAx(1,iii)=F5.fwd_model.nodes(486+4*iii,1);
        BBBx(1,iii)=F5.fwd_model.nodes(487+4*iii,1);
    end
    for iii=15:16
        AAAx(1,iii)=F5.fwd_model.nodes(478+4*(iii-14),1);
        BBBx(1,iii)=F5.fwd_model.nodes(479+4*(iii-14),1);
    end
    
    electrode_average.x=[AAAx;BBBx]/5.8;
    for iii=1:14
        AAAy(1,iii)=F5.fwd_model.nodes(486+4*iii,2);
        BBBy(1,iii)=F5.fwd_model.nodes(487+4*iii,2);
    end
    for iii=15:16
        AAAy(1,iii)=F5.fwd_model.nodes(478+4*(iii-14),2);
        BBBy(1,iii)=F5.fwd_model.nodes(479+4*(iii-14),2);
    end
    electrode_average.y=[AAAy;BBBy]/5.8;
    dimension1=floor(Dy/dy+1);
    electrode_possible.y1=zeros(dimension1,7); electrode_possible.y2=zeros(dimension1,7);
    for el=6:12
        %         electrode_possible.y1(:,el-5)=electrode_average.y(1,el)+Dy/2:-dy:electrode_average.y(1,el)-Dy/2;
        %         electrode_possible.y2(:,el-5)=electrode_average.y(2,el)+Dy/2:-dy:electrode_average.y(2,el)-Dy/2;
        electrode_possible.y1(:,el-5)=flip(linspace(electrode_average.y(1,el)-Dy/2,electrode_average.y(1,el)+Dy/2,dimension1));
        electrode_possible.y2(:,el-5)=flip(linspace(electrode_average.y(2,el)-Dy/2,electrode_average.y(2,el)+Dy/2,dimension1));
    end
%     plot(electrode_average.x(1,9)*ones(dimension1,1),electrode_possible.y1(:,4),'o')
%     hold on
%     plot(electrode_average.x(2,9)*ones(dimension1,1),electrode_possible.y2(:,4),'o')
    o=[];
    
    
elseif N==32
    F5=mk_common_model('d2T3',32);
    AAAx=zeros(1,14);
    BBBx=AAAx;
    AAAy=zeros(1,14);
    BBBy=AAAx;
    for iii=1:28
        AAAx(1,iii)=F5.fwd_model.nodes(488+2*iii,1);
        BBBx(1,iii)=F5.fwd_model.nodes(489+2*iii,1);
    end
    for iii=29:32
        AAAx(1,iii)=F5.fwd_model.nodes(480+2*(iii-28),1);
        BBBx(1,iii)=F5.fwd_model.nodes(481+2*(iii-28),1);
    end
    
    electrode_average.x=[AAAx;BBBx]/5.8;
    for iii=1:28
        AAAy(1,iii)=F5.fwd_model.nodes(488+2*iii,2);
        BBBy(1,iii)=F5.fwd_model.nodes(489+2*iii,2);
    end
    for iii=29:32
        AAAy(1,iii)=F5.fwd_model.nodes(480+2*(iii-28),2);
        BBBy(1,iii)=F5.fwd_model.nodes(481+2*(iii-28),2);
    end
    electrode_average.y=[AAAy;BBBy]/5.8;
    dimension1=floor(Dy/dy+1);
    electrode_possible.y1=zeros(dimension1,15); electrode_possible.y2=zeros(dimension1,15);
    for el=10:24
        %electrode_possible.y1(:,el-9)=electrode_average.y(1,el)+Dy/2:-dy:electrode_average.y(1,el)-Dy/2;
        %electrode_possible.y2(:,el-9)=electrode_average.y(2,el)+Dy/2:-dy:electrode_average.y(2,el)-Dy/2;
        electrode_possible.y1(:,el-9)=flip(linspace(electrode_average.y(1,el)-Dy/2,electrode_average.y(1,el)+Dy/2,dimension1));
        electrode_possible.y2(:,el-9)=flip(linspace(electrode_average.y(2,el)-Dy/2,electrode_average.y(2,el)+Dy/2,dimension1));
    end
end
electrode_average.x=k*electrode_average.x;
electrode_average.y=k*electrode_average.y;
electrode_possible.y1=k*electrode_possible.y1;
electrode_possible.y2=k*electrode_possible.y2;
o=[];
end

function [Left_lung,Right_lung,Heart,Muscles,Bones,r,k]=...
    define_tissues(r,k,times,discrete_points,firsttime,f,end_Admittances,randomcond,Tissues)
%defines tissues properties
if exist('Tissues')~=0 %tissues exist but need to be updated (end of cycle)
    %in order to be updated we have to temporarily store their positions
    Left_lung=Tissues.Left_lung;
    pos1=Left_lung.position;
    Right_lung=Tissues.Right_lung;
    pos2=Right_lung.position;
    Exterior=Tissues.Heart.exterior;
    pos3a=Exterior.position;
    Interior=Tissues.Heart.interior;
    pos3b=Interior.position;
    Muscles=Tissues.Muscles;
    Bones=Tissues.Bones;
    pos4=Bones.Bone1.position;
    pos5=Bones.Bone2.position;
    pos6=Bones.Bone3.position;
    pos7=Bones.Bone4.position;
    pos8=Bones.Bone5.position;
    pos9=Bones.Bone6.position;
    pos10=Bones.Bone7.position;
    pos11=Bones.Bone8.position;
end
%%%% unpack end_Admittances
slairmax=end_Admittances(1);
slairmin=end_Admittances(2);
dslblood=end_Admittances(3);
elairmax=end_Admittances(4);
elairmin=end_Admittances(5);
delblood=end_Admittances(6);
Sheart=end_Admittances(7);
dSheart=end_Admittances(8);
Eheart=end_Admittances(9);
dEheart=end_Admittances(10);

Left_lung=define_lung_properties(times,discrete_points,slairmax,slairmin,dslblood,elairmax,elairmin,delblood);
Right_lung=define_lung_properties(times,discrete_points,slairmax,slairmin,dslblood,elairmax,elairmin,delblood);
Heart=define_heart_properties(f,times,Sheart,dSheart,Eheart,dEheart);
Muscles=define_muscles_properties(times);
Bones=define_bones_properties();

if exist('Tissues')~=0 %tissues exist but need to be updated (end of cycle)
    Left_lung.position=pos1;
    Right_lung.position=pos2;
    Heart.exterior.position=pos3a;
    Heart.interior.position=pos3b;
    Bones.Bone1.position=pos4;
    Bones.Bone2.position=pos5;
    Bones.Bone3.position=pos6;
    Bones.Bone4.position=pos7;
    Bones.Bone5.position=pos8;
    Bones.Bone6.position=pos9;
    Bones.Bone7.position=pos10;
    Bones.Bone8.position=pos11;
end

[Left_lung,Right_lung,Heart,Bones,r,k]=...
    define_tissues_movement(r,k,Left_lung,Right_lung,Heart,Bones,times,discrete_points,firsttime,randomcond);
end


function Lung=define_lung_properties(times,discrete_points,slairmax,slairmin,dslblood,elairmax,elairmin,delblood)
Lung.props.Volume=set_Volume(discrete_points.total_positions);
Lung.props.Conductivity=set_lung_Conductivity(times,Lung.props.Volume,slairmax,slairmin,dslblood);
Lung.props.permittivity=set_lung_permittivity(times,Lung.props.Volume,elairmax,elairmin,delblood);
end


function Heart=define_heart_properties(f,times,Sheart,dSheart,Eheart,dEheart)
Heart.interior.props.Conductivity=set_heart_Conductivity(times,Sheart,dSheart);
Heart.interior.props.permittivity=set_heart_permittivity(times,Eheart,dEheart);
Heart.exterior.props.Conductivity=set_myocardium_Conductivity(f);
Heart.exterior.props.permittivity=set_myocardium_permittivity(f);
end

function Volume=set_Volume(total_positions)
%during the breath, each tissue and the electrodes, pass from each
%possible position twice
i2=total_positions:1:2*total_positions+1; %time samples of exhalation
i1=1:1:(total_positions-1); %time samples of inhalation
F2=-3/(total_positions-1)*i1+3/(total_positions-1)+4; %lungs volume during inhalation
F1=3/(total_positions+1)*i2+4-3*(2*total_positions+1)/(total_positions+1); %lungs volume during exhalation
i=[i1 i2]; %combine the two continuous time windows
Volume=[F1 F2]; %combine the two continuous volume functions
end

function S=set_lung_Conductivity(times,Volume,slairmax,slairmin,dslblood)
%lungs parameters
w=1.5;
sb=0.5;
si=2;
% conductivity changes during breath cycle due to air flow
s10lungs=18*(32*Volume+4.5)./(32*Volume+9).^2*((0.85/w*sb+0.03*si));
%scaling (K1, K2)
s1lungs=(slairmax+slairmin)/2+(slairmax-(slairmax+slairmin)/2)*(s10lungs-mean(s10lungs))/(max(s10lungs-mean(s10lungs)));
pulse_dur=times.pulse_dur;
%dslblood=0.008; %ds in lungs due to blood cycle
j1=1:1:pulse_dur/2;
j2=pulse_dur/2:1:pulse_dur;
%conductivity changes during blood cycle due to blood flow
sklungs=dslblood*2/pulse_dur*j1;
sllungs=-dslblood*2/pulse_dur*j2+2*dslblood;
j=[j1 j2];
s2lungs=[sklungs sllungs];
S.air=s1lungs;
S.blood=s2lungs;
end

function E=set_lung_permittivity(times,Volume,elairmax,elairmin,delblood)
erb=10000;
erm=10;
b=0.325;
w=1.5;
% permittivity changes during breath cycle due to air flow
er10lungs=18*(32*Volume+4.5)./(32*Volume+9).^2.*((0.85/w*erb+2400*b*Volume.^(1/3)*erm));
%scaling (L1, L2)
er1lungs=(elairmax+elairmin)/2+(elairmax-(elairmax+elairmin)/2)*(er10lungs-mean(er10lungs))/(max(er10lungs-mean(er10lungs)));
pulse_dur=times.pulse_dur;
dElungs=delblood;
%dElungs=500; %ÄE in lungs due to blood cycle
j1=1:1:pulse_dur/2;
j2=pulse_dur/2:1:pulse_dur;
j=[j1 j2];
% permittivity changes during blood cycle due to blood flow
er2klungs=dElungs*2/pulse_dur*j1-2*dElungs/pulse_dur;
er2llungs=-dElungs*2/pulse_dur*j2+2*dElungs;
er2lungs=[er2klungs er2llungs];
E.air=er1lungs;
E.blood=er2lungs;
end

function S=set_heart_Conductivity(times,Sinheart,dSheart)
pulse_dur=times.pulse_dur;
j1=1:1:pulse_dur/2;
j2=pulse_dur/2:1:pulse_dur;
j=[j1 j2];
%heart parameters
%dSheart=0.025; %conductivity change of heart during cardio cycle
%Sinheart=0.5; %sigma heart~0.5, around is ~0.2
sheart1=-2*dSheart/pulse_dur*j1+dSheart+Sinheart;
sheart2=2*dSheart/pulse_dur*j2-dSheart+Sinheart;
S=[sheart1 sheart2];
o=[];
end

function E=set_heart_permittivity(times,Eheartin,dEheart)
pulse_dur=times.pulse_dur;
j1=1:1:pulse_dur/2;
j2=pulse_dur/2:1:pulse_dur;
j=[j1 j2];
% dEheart=300;
% Eheartin=6000;
er1heart=-dEheart*2/pulse_dur*j1+dEheart+Eheartin;
er2heart=dEheart*2/pulse_dur*j2-dEheart+Eheartin;
E=[er1heart er2heart];
end

function S=set_myocardium_Conductivity(f)
S=9*10^(-5)*f/1000+0.168;
end

function E=set_myocardium_permittivity(f)
E=2.29*10^5*(f/1000)^(-0.95);
end


function Muscles=define_muscles_properties(times)
Muscles.props.Conductivity=define_muscle_conductivity(times);
Muscles.props.permittivity=define_muscle_permittivity(times);

end

function S=define_muscle_conductivity(times)
%soft_tissues parameters
pulse_dur=times.pulse_dur;
j1=1:1:pulse_dur/2;
j2=pulse_dur/2:1:pulse_dur;
dSsoft=0.015; %muscle-plasma conductivity change due to blood flow
Sinsoft=0.28; %average muscle-plasma-gap conductivity
ssoft1=2*dSsoft/times.pulse_dur*j1+dSsoft+Sinsoft-2*dSsoft/times.pulse_dur;
ssoft2=-2*dSsoft/times.pulse_dur*j2+dSsoft+Sinsoft+2*dSsoft;
S=[ssoft1 ssoft2];
end

function E=define_muscle_permittivity(times)
%soft_tissues parameters
pulse_dur=times.pulse_dur;
j1=1:1:pulse_dur/2;
j2=pulse_dur/2:1:pulse_dur;
dEsoft=300;
Esoftin=6000;
er1soft=-dEsoft*2/pulse_dur*j1+dEsoft+Esoftin+2*dEsoft/pulse_dur;
er2soft=dEsoft*2/pulse_dur*j2+dEsoft+Esoftin-2*dEsoft;
E=[er1soft er2soft];
end

function Bones=define_bones_properties()
%assuming constant conductivity and permittivity
Bones.props.Conductivity=0.07;
Bones.props.permittivity=250;
end

function [Left_lung,Right_lung,Heart,Bones,r,ks]=...
    define_tissues_movement(r,ks,Left_lung,Right_lung,Heart,Bones,times,discrete_points,firsttime,randomcond)
%Defines:
%initially the minimum and maximum shapes of the tissues (assuming they
%have constant limits)
%at every new blood or breath cycle the steps of each tissue's boundary
%change
%r = H + (H-L)*rand(DIM,1)
%Left Lung
Left_lung.shape.full=r*[-5.1 -7.6 -8.3 -8.1 -8 -5.6 -4 -1.3 -10.1 -16.6 -19.8 -20.9 -20 -15.8 -11.8 -8; -12.8 -10.3 -7.7 -3.2 -0.1 3.2 7 9.7 15.5 10.9 6.3 1.3 -5.2 -10.2 -12.5 -13.1]';
Left_lung.shape.empty=r*[-6.4 -7.6 -8.8 -9.7 -10.4 -8.4 -6.9 -5 -11 -14.9 -17.8 -18.5 -17.4 -14.8 -10.9 -8.7;-10.9 -8.2 -5.6 -2.9 -0.3 3.7 5.7 7.5 13.6 9.3 5.3 1.5 -4.4 -8.1 -11.6 -11.7]';
Left_lung.shape.step=define_steps(Left_lung,discrete_points.total_positions);
%Right lung
Right_lung.shape.full=r*[10.5 18.8 21.5 18.3 12 6.2 1.9 3.3 5.2 7.6 9.2 8 7.5 7.3 6.9; -15.6 -9.1 0.4 8.8 15.7 16.1 10.1 6.2 3.4 3.5 -0.7 -4.1 -5.2 -6.9 -10.5]';
Right_lung.shape.empty=r*[10.8 17.7 19 17 11.2 6.6 3.7 4.4 5.3 7.6 10.3 9.6 8.4 9 8.7;-12.8 -7 1.6 8.7 14.6 15 10.3 7.4 4.8 4.8 0.9 -4.1 -5.1 -6.6 -9.5]';
Right_lung.shape.step=define_steps(Right_lung,discrete_points.total_positions);
%Heart exterior
Heart.exterior.shape.full=[0 4.1 4.5 4.5 5.9 7.3 7.1 4.3 3 0.8 1.3 1.6 -1 -2.8 -3.6 -6.7 -5.9 -1.1;-14.6 -10.5 -6.6 -3.6 -3.3 -1.7 0.5 2.6 -0.4 1 1.3 6.1 6 4 1.3 0 -7.8 -12.4]';
Heart.exterior.shape.empty=[0 3.6 3.1 3.3 5 6.4 5.9 4.2 3.1 1 0.6 0.6 -0.9 -2.9 -4.1 -6.3 -4.9 -0.8;-13 -8.9 -6.2 -3.3 -2.3 -1.1 0.1 1.6 -1.5 0 2.2 4.2 5.9 3.8 1.7 -0.2 -7.2 -11.8]';

Heart.exterior.shape.full(:,2)=Heart.exterior.shape.full(:,2)+1;
Heart.exterior.shape.empty(:,2)=Heart.exterior.shape.empty(:,2)+1;

Heart.exterior.shape.full=0.8*r*Heart.exterior.shape.full;
Heart.exterior.shape.empty=0.8*r*Heart.exterior.shape.empty;

%Heart interior
aa=[1 2 4 9 15 17];
Heart.interior.shape.full=0.65*Heart.exterior.shape.full(aa,:);
Heart.interior.shape.empty=0.65*Heart.exterior.shape.empty(aa,:);

%heart exterior shifting
Heart.exterior.shape.full(:,2)=Heart.exterior.shape.full(:,2)-2;
Heart.exterior.shape.empty(:,2)=Heart.exterior.shape.empty(:,2)-2;

Heart.exterior.shape.step=define_steps(Heart.exterior,discrete_points.total_positions*times.tbr/times.Tpulse);

%Heart interior shifting
Heart.interior.shape.full(:,2)=Heart.interior.shape.full(:,2)-2.5;
Heart.interior.shape.empty(:,2)=Heart.interior.shape.empty(:,2)-2.5;

Heart.interior.shape.step=define_steps(Heart.interior,discrete_points.total_positions*times.tbr/times.Tpulse);



%Bones
Bones.Bone1.shape.empty=ks*[21.8 22.9 21.7 20.3; -7.7 -3.1 -2 -7.6]';
Bones.Bone2.shape.empty=ks*[21.8 21.7 16.1 14.3; 6 8.1 14.5 14.4]';
%Spondilus
Bones.Bone3.shape.empty=ks*[11.5 3.4 1.2 0.4 -0.3 -0.8 -2.6 -5.2 -12.3 -3.5 -0.4 2.5; 17.1 18.2 17.8 17.9 19.2 18 17.5 17.8 16.5 14.3 10.6 14.5]';
Bones.Bone4.shape.empty=ks*[-14.2 -18.1 -22.2 -22.7;14.5 12.7 6.9 3.8]';
Bones.Bone5.shape.empty=ks*[-20.8 -21.7 -22.8 -22.1;-5.9 -1.2 -1.3 -6]';
Bones.Bone6.shape.empty=ks*[-15.8 -17 -11.4 -11.8; -13.5 -13.6 -17.2 -15.9]';
Bones.Bone7.shape.empty=ks*[-2.8 -2.7 3.7 3.7;-15.6 -17.4 -17.1 -15.6]';
Bones.Bone8.shape.empty=ks*[12.5 12.7 18 16.8;-16.1 -17.4 -13.4 -12.7]';
Bones.Bone1.shape.full=ks*[22.9 23.7 22.8 21.4; -8 -3.8 -2.6 -7.7]';
Bones.Bone2.shape.full=ks*[22.1 22.5 16.8 14.8; 6.5 8.1 15.1 15]';
%Spondilus
Bones.Bone3.shape.full=ks*[12.1 3.3 1.1 0.4 -0.3 -0.8 -2.6 -5.2 -12.9 -3.2 -0.4 2.7; 17.9 19.1 18.4 18.2 19 18 17.9 18.4 17 15 10.6 15.5]';
Bones.Bone4.shape.full=ks*[-14.5 -18.7 -22.7 -23; 15 12.6 6.9 4]';
Bones.Bone5.shape.full=ks*[-21.3 -22.4 -23.4 -22.6; -6 -1.4 -1.4 -6.2]';
Bones.Bone6.shape.full=ks*[-16.6 -17.8 -11.2 -11.6;-13.7 -13.9 -17.8 -16.8]';
Bones.Bone7.shape.full=ks*[-2.6 -2.4 3.7 3.9; -16.3 -18 -17.9 -16.5]';
Bones.Bone8.shape.full=ks*[12.7 12.8 18.7 17.6; -17 -18.2 -13.8 -13.1]';
for k=1:8
    Bone=eval(['Bones.Bone' num2str(k)]);
    eval(['Bones.Bone' num2str(k) '.shape.step=define_steps(Bone,discrete_points.total_positions)'])
end
if firsttime==0&&randomcond==1
    while check_intersections(Left_lung,Right_lung,Heart,Bones)==1
        %if we have just initialized the tissues and we have intersections
        %re-initialize the positions
        r=0.95 + (1.05-0.95)*randn(1,1);
        ks=0.95 + (1.05-0.95)*randn(1,1);
        [Left_lung,Right_lung,Heart,Bones,r,ks]=define_tissues_movement(r,ks,Left_lung,Right_lung,Heart,Bones,times,discrete_points,firsttime,0);
    end
end
end

function step=define_steps(Tissue,total_positions)
len=length(Tissue.shape.empty(:,1));
step(1:len,1)=(Tissue.shape.empty(:,1)-Tissue.shape.full(:,1))/total_positions;
step(1:len,2)=(Tissue.shape.empty(:,2)-Tissue.shape.full(:,2))/total_positions;
end

function [fileID,d,folds]=open_script(frequency,path)
d5=1;
folder=[path '\Dynamic_Thorax_Model\LUA_files\'];
prefix_data=['File',num2str(date()),'_',num2str(d5)];
dataformat='.lua';
nol=strcat(folder,prefix_data,dataformat);
while exist(nol)>0
    folder=[path '\Dynamic_Thorax_Model\LUA_files\'];
    prefix_data=['File',num2str(date()),'_',num2str(d5)];
    dataformat='.lua';
    nol=strcat(folder,prefix_data,dataformat);
    d5=d5+1;
end
%% begin writing on script
fileID=fopen(nol,'w');
%% set parameters of the problem in FEMM
fprintf(fileID,'newdocument(3)\n');
fprintf(fileID,'ci_probdef("centimeters","planar",%4.2f,1.E-8,1,15)\n',frequency);
%% set file for the measurements to be written
d=1;
ex=1;
while ex==1
    folder=[path 'Dynamic_Thorax_Model\LUA_files\'];
    prefix_data=[num2str(date()),'_test',num2str(d)];
    dataformat='.fee';
    nl=strcat(folder,prefix_data,dataformat);
    if exist(nl)>0
        d=d+1;
        %disp('exist')
    else
        ex=0;
    end
end
%% save it
folds=strrep(folder, '\', '\\\\');
folds=[folds 'Temporary\\\\'];
fprintf(fileID,'ci_saveas("%s%s_test%.0f.fee")\n',folds,date(),d);
end

function []=set_materials_in_lua(fileID,Current,Tissues,Zelectrodes)
inhomogeneous_materials_stable(fileID,Tissues,Zelectrodes)
%%set current values
set_currents(fileID,Current)
%%set ground value
set_ground(fileID)
end

function[]=inhomogeneous_materials_stable(fileID,Tissues,Zelectrodes)
%sets some constant values on tissues; only bones will be used
%trabecular bone resistance 2000Ùcm cortical 10000Ùcm at 10kHz
%20 Ùm===>1/20Sm^-1
%10000 Ùcm=0.01Sm^-1
fprintf(fileID,'ci_addmaterial("Lungs",0.1,0.1,5000,5000,"<None>","<None>")\n');
fprintf(fileID,'ci_addmaterial("Bone",0.07,0.07,250,250,"<None>","<None>")\n');
fprintf(fileID,'ci_addmaterial("Heart",0.3,0.3,5000,5000,"<None>","<None>")\n');
fprintf(fileID,'ci_addmaterial("Softer tissues",3,3,5000,5000,"<None>","<None>")\n');
fprintf(fileID,'ci_addmaterial("Near_skin_",0.3,0.3,2000,2000,"<None>","<None>")\n');
Sigmaext=Tissues.Heart.exterior.props.Conductivity;
Eext=Tissues.Heart.exterior.props.permittivity;
fprintf(fileID,'ci_addmaterial("Miocardium",%4.2f,%4.2f,%4.2f,%4.2f,"<None>","<None>")\n',Sigmaext,Sigmaext,Eext,Eext);
fprintf(fileID,'ci_addmaterial("Fat",%4.2f,%4.2f,%4.2f,%4.2f,"<None>","<None>")\n',0.12,0.12,1000,1000);
for electrode=1:length(Zelectrodes(:,1))
    fprintf(fileID,'ci_addmaterial("Zelectrode%2.0f",%4.2f,%4.2f,%4.2f,%4.2f,"<None>","<None>")\n',...
        electrode,Zelectrodes(electrode,1),Zelectrodes(electrode,1),Zelectrodes(electrode,2),Zelectrodes(electrode,2));
end
end

function[]=set_currents(fileID,Current)
%sets the value of current injection amplitude
fprintf(fileID,'ci_addconductorprop("Iin",0,%1.3f,0)\n',Current);
fprintf(fileID,'ci_addconductorprop("Iout",0,-%1.3f,0)\n',Current);
end

function[]=set_ground(fileID)
%Ground is 0 volts!
fprintf(fileID,'ci_addboundprop("Ground", 0, 0, 0, 0, 0)\n');
end

function [pulmo,card,pulmcounter,cardcounter,pulse_dur]=...
    randomize(initial_times,initial_total_positions,random)

%returns the initial states of breathing and cardiac cycles
%pulmo and card indicate the positions of lungs and heart
%pulmcounter and cardcounter indicate whether they are increasing or
%decreasing
tot=initial_total_positions.total_positions;
tbr=initial_times.tbr;
Tpulse=initial_times.Tpulse;
if random==1 %pulmo cycle and cardio cycle are at random points when simulation starts
    kappa=tot*rand(1); %random ventilation point
    cof1=rand(1); %determine if we have exhalation (-1) or inhalation (+1)
    if cof1<0.5
        cof=1;
    else
        cof=-1;
    end
    cof2=rand(1);
    if cof2<0.5
        cardio=1;
    else
        cardio=-1;
    end
    i=round(kappa);
    pulse_dur=tot*(tbr/Tpulse)^(-1);
    kappa1=rand(1)*pulse_dur;
    j=round(kappa1);
    pulmo=i;
    card=j;
    pulmcounter=cof;
    cardcounter=cardio;
else  %if random is non-one when simulation starts, pulmocycle and cardiocycle
    %are set to starting point
    pulmo=1;
    card=1;
    pulmcounter=1;
    cardcounter=1;
    pulse_dur=tot*(tbr/Tpulse)^(-1);
end
end

function Graph=initialize_graph()
%input conductivities zero defining
Graph.heart.conductivity=[];
Graph.right_lung.conductivity=[];
Graph.left_lung.conductivity=[];
Graph.muscles.conductivity=[];
Graph.heart.permittivity=[];
Graph.right_lung.permittivity=[];
Graph.left_lung.permittivity=[];
Graph.muscles.permittivity=[];
end

function Silences=initialize_silences()
Silences.silence=0; %when starting at a random point, not in silence space
Silences.time_sil_start=-1; %start time (seconds) of the silence space
Silences.relaxing_time=-1; %duration (seconds) of the silence space
Silences.cardiobeats=1;
end
%% loop
function [Graph,timesteps,total_time_est,History_of_Tissues]=...
    frame_loop(N,currentskip,initial_times,initial_total_positions,...
Collapsion_areas,electrodes,Silences,fileID,Graph, Circulation,...
Ventilation, Tissues,d,r,k,frequency,end_Admittances,folds,randomcond)
times=initial_times;
tim=initial_times.tim;
total_positions=initial_total_positions;
total_time_est=initial_times.total_time_est;
count=1;
timesteps=1;
History_of_Tissues=initialize_history();
%step is (time required for a frame)/(number of electrodes)
%therefore (time) step is the time we set a current input
fprintf(fileID,'ci_zoom(%d,%d,%d,%d) \n',-30,-30,30,30);
for time=(tim/N):(tim/N):total_time_est
    %each iteration has a tim/N duration which corresponds
    %to the time the measurements are taken for a single
    %current source electrodes position
    if mod(count-1,N)==0
        count=1;
    end
    % find  SIGMA and E in (ith) state
    Graph=inhomogeneous_materials_var_setting(Graph, Circulation, Ventilation, Tissues, fileID);
    %detect positions of each electrode and tissue
    %and define the geometry which changes in every time step
    %find the electrodes and tissues POSITIONS in (ith) state
    electrodes.position=electrode_positions(N,Ventilation,electrodes);
    Tissues=tissue_positions(Ventilation,Circulation,Tissues,Collapsion_areas);
    History_of_Tissues=update_history(History_of_Tissues,Tissues);
    %set assist temporary points to measure later
    [assist_points,electrode_central_points,electrode_exterior_sides]=dynamic_geometry(N,electrodes,Tissues,Collapsion_areas,fileID);
    %add labels to each region according to current (ith) state
    [LLCC,RLCC,CALCC,CARCC,STCC,Xbone,Ybone,Xfat,Yfat,Xint,Yint,Xext,Yext]=...
        addmaterial_labels(Tissues,Ventilation,Circulation,Collapsion_areas,fileID,r,k,electrode_central_points);
    %set current source position for the (ith) state
    set_current_source(N,currentskip,assist_points,count,fileID,electrode_exterior_sides);
    %set reference (arbitrary ground) electrode
    set_reference(N,currentskip,assist_points,count,fileID,electrode_exterior_sides);
    % save the file
    save_process(d,fileID,folds);
    % measure the voltages for this state (this current source position)
    measure_voltages(N,electrodes.position,d,timesteps,fileID,electrode_exterior_sides,folds);
    %update cycle states; then check if any of them has come to an end
    %if so, change durations, total positions, properties, e.t.c.
    delec=electrodes.delec;
    [Ventilation,Circulation,total_positions,times,Silences,Tissues,electrodes,r,k]=...
        update_cycle_states(N,Ventilation,Circulation,total_positions,times,Silences,...
        Tissues,time,electrodes,r,k,frequency,end_Admittances,delec,randomcond);
    display=0;
%      plot(Tissues.Bones.Bone7.position(:,1),Tissues.Bones.Bone7.position(:,2))
%         hold on 
%      plot(electrodes.position.x(:,9),electrodes.position.y(:,9),'o')
%      hold on
    if isfield(Tissues.Left_lung,'position')==0
        o=[];
    end
    if display==1&&mod(timesteps,15)==0
        RTG=display_realtime_graph(N,Graph,Ventilation,Circulation,Tissues,electrodes,time,timesteps,Collapsion_areas,times);
    end
    %% at the iteration end unset FEMM geometry
    %because it's gonna change in the next step/iteration
    unset_current_source(N,fileID,currentskip,assist_points,count,electrode_exterior_sides)
    unset_reference(N,fileID,currentskip,assist_points,count,electrode_exterior_sides)
    undo_dynamic_geometry(Tissues,fileID,electrodes.position,electrode_exterior_sides)
    clear_labels(fileID,LLCC,RLCC,CALCC,CARCC,STCC,Xbone,Ybone,Xfat,Yfat,Xint,Yint,Xext,Yext,electrode_central_points)
    count=count+1;
    timesteps=timesteps+1;
    if mod(timesteps,100)==0
        fprintf('Loading script %2.0f %%\n',round(100*time/total_time_est))
    end
end
end
%%
function History_of_Tissues=initialize_history()
History_of_Tissues.Left_lung=[];
History_of_Tissues.Right_lung=[];
History_of_Tissues.Heart.exterior=[];
History_of_Tissues.Heart.interior=[];
for Bone=1:8
    B_k=[];
    switch Bone
        case 1
            History_of_Tissues.Bones.Bone1=B_k;
        case 2
            History_of_Tissues.Bones.Bone2=B_k;
        case 3
            History_of_Tissues.Bones.Bone3=B_k;
        case 4
            History_of_Tissues.Bones.Bone4=B_k;
        case 5
            History_of_Tissues.Bones.Bone5=B_k;
        case 6
            History_of_Tissues.Bones.Bone6=B_k;
        case 7
            History_of_Tissues.Bones.Bone7=B_k;
        case 8
            History_of_Tissues.Bones.Bone8=B_k;
    end
end
end

function History_of_Tissues=update_history(History_of_Tissues,Tissues)
History_of_Tissues.Left_lung=[History_of_Tissues.Left_lung Tissues.Left_lung.position];
History_of_Tissues.Right_lung=[History_of_Tissues.Right_lung Tissues.Right_lung.position];
History_of_Tissues.Heart.exterior=[History_of_Tissues.Heart.exterior Tissues.Heart.exterior.position];
History_of_Tissues.Heart.interior=[History_of_Tissues.Heart.interior Tissues.Heart.interior.position];
for Bone=1:8
    B_k=eval(['Tissues.Bones.Bone',num2str(Bone),'.position']);
    switch Bone
        case 1
            History_of_Tissues.Bones.Bone1=[History_of_Tissues.Bones.Bone1 B_k];
        case 2
            History_of_Tissues.Bones.Bone2=[History_of_Tissues.Bones.Bone2 B_k];
        case 3
            History_of_Tissues.Bones.Bone3=[History_of_Tissues.Bones.Bone3 B_k];
        case 4
            History_of_Tissues.Bones.Bone4=[History_of_Tissues.Bones.Bone4 B_k];
        case 5
            History_of_Tissues.Bones.Bone5=[History_of_Tissues.Bones.Bone5 B_k];
        case 6
            History_of_Tissues.Bones.Bone6=[History_of_Tissues.Bones.Bone6 B_k];
        case 7
            History_of_Tissues.Bones.Bone7=[History_of_Tissues.Bones.Bone7 B_k];
        case 8
            History_of_Tissues.Bones.Bone8=[History_of_Tissues.Bones.Bone8 B_k];
    end
end
end
%%
function Graph=inhomogeneous_materials_var_setting(Graph, Circulation, Ventilation, Tissues, fileID)
%Writes the tissues properties of the current state(s) as FEMM materials
%Returns an updated graph
%this is updated at each time step

s1lungs=Tissues.Left_lung.props.Conductivity.air;
s2lungs=Tissues.Left_lung.props.Conductivity.blood;
%fprintf('%2.0f\n',Ventilation.pulmocycle)
if Ventilation.pulmocycle<=0
    o=[];
end
left_lung_cond_value=s1lungs(Ventilation.pulmocycle)+s2lungs(Circulation.cardio_cycle);
Graph.left_lung.conductivity=[Graph.left_lung.conductivity left_lung_cond_value];

er1lungs=Tissues.Left_lung.props.permittivity.air;
er2lungs=Tissues.Left_lung.props.permittivity.blood;

left_lung_perm_value=er1lungs(Ventilation.pulmocycle)+er2lungs(Circulation.cardio_cycle);
Graph.left_lung.permittivity=[Graph.left_lung.permittivity left_lung_perm_value];

cheart=Tissues.Heart.interior.props.Conductivity(Circulation.cardio_cycle); errheart=Tissues.Heart.interior.props.permittivity(Circulation.cardio_cycle);
Graph.heart.conductivity=[Graph.heart.conductivity cheart];
Graph.heart.permittivity=[Graph.heart.permittivity errheart];

csoft=Tissues.Muscles.props.Conductivity(Circulation.cardio_cycle); errsoft=Tissues.Muscles.props.permittivity(Circulation.cardio_cycle);
Graph.muscles.conductivity=[Graph.muscles.conductivity csoft];
Graph.muscles.permittivity=[Graph.muscles.permittivity errsoft];


fprintf(fileID,'ci_addmaterial("Lungs_state_%d_%d",%4.2f,%4.2f,%4.2f,%4.2f,"<None>","<None>")\n',Ventilation.pulmocycle,Circulation.cardio_cycle,...
    left_lung_cond_value,left_lung_cond_value,left_lung_perm_value,left_lung_perm_value);
fprintf(fileID,'ci_addmaterial("Heart_state_%d",%4.2f,%4.2f,%4.2f,%4.2f,"<None>","<None>")\n',Circulation.cardio_cycle,cheart,cheart,errheart,errheart);
fprintf(fileID,'ci_addmaterial("Soft_tissue_state_%d",%4.2f,%4.2f,%4.2f,%4.2f,"<None>","<None>")\n',Circulation.cardio_cycle,csoft,csoft,errsoft,errsoft);
end

function electrode_position=electrode_positions(N,Ventilation,electrodes)
%Returns this state's electrode positions
electrode_average=electrodes.electrode_average;
electrode_possible=electrodes.electrode_possible;
i=Ventilation.pulmocycle;

if N==16
    electrode_position.y(:,1:5)=electrode_average.y(:,1:5);
    electrode_position.y(:,12:16)=electrode_average.y(:,12:16);
    electrode_position.x=electrode_average.x;
    electrode_position.y(1,6:12)=electrode_possible.y1(i,1:7);
    electrode_position.y(2,6:12)=electrode_possible.y2(i,1:7);
elseif N==32
    electrode_position.y(:,1:9)=electrode_average.y(:,1:9);
    electrode_position.y(:,24:32)=electrode_average.y(:,24:32);
    electrode_position.x=electrode_average.x;
    electrode_position.y(1,10:24)=electrode_possible.y1(i,1:15);
    electrode_position.y(2,10:24)=electrode_possible.y2(i,1:15);
end
end

%%
function Tissues=tissue_positions(Ventilation,Circulation,Tissues,Collapsion_areas)

%Returns each tissue position at current state (time) i
LL=Tissues.Left_lung;
LL.position=[LL.shape.empty(:,1)-LL.shape.step(:,1)*Ventilation.pulmocycle,LL.shape.empty(:,2)-LL.shape.step(:,2)*Ventilation.pulmocycle];
switch Collapsion_areas.collapsions.collapsion_left_lung
    case 1
        %fix collapsed area positions at the middle of empty and full
        %positions
        LL.position(14:16,1)=1/2*(LL.shape.empty(14:16,1)+LL.shape.full(14:16,1));
        LL.position(14:16,2)=1/2*(LL.shape.empty(14:16,2)+LL.shape.full(14:16,2));
        LL.position(1:4,1)=1/2*(LL.shape.empty(1:4,1)+LL.shape.full(1:4,1));
        LL.position(1:4,2)=1/2*(LL.shape.empty(1:4,2)+LL.shape.full(1:4,2));
    case 2
        LL.position(14:16,1)=1/2*(LL.shape.empty(14:16,1)+LL.shape.full(14:16,1));
        LL.position(14:16,2)=1/2*(LL.shape.empty(14:16,2)+LL.shape.full(14:16,2));
        LL.position(1:3,1)=1/2*(LL.shape.empty(1:3,1)+LL.shape.full(1:3,1));
        LL.position(1:3,2)=1/2*(LL.shape.empty(1:3,2)+LL.shape.full(1:3,2));
    case 3
        LL.position(15:16,1)=1/2*(LL.shape.empty(15:16,1)+LL.shape.full(15:16,1));
        LL.position(15:16,2)=1/2*(LL.shape.empty(15:16,2)+LL.shape.full(15:16,2));
        LL.position(1:3,1)=1/2*(LL.shape.empty(1:3,1)+LL.shape.full(1:3,1));
        LL.position(1:3,2)=1/2*(LL.shape.empty(1:3,2)+LL.shape.full(1:3,2));
    case 4
        LL.position(16,1)=1/2*(LL.shape.empty(16,1)+LL.shape.full(16,1));
        LL.position(16,2)=1/2*(LL.shape.empty(16,2)+LL.shape.full(16,2));
        LL.position(1:2,1)=1/2*(LL.shape.empty(1:2,1)+LL.shape.full(1:2,1));
        LL.position(1:2,2)=1/2*(LL.shape.empty(1:2,2)+LL.shape.full(1:2,2));
    otherwise
        %do nothing
end
%Substitution
Tissues.Left_lung=LL;

RL=Tissues.Right_lung;
RL.position=[RL.shape.empty(:,1)-RL.shape.step(:,1)*Ventilation.pulmocycle,RL.shape.empty(:,2)-RL.shape.step(:,2)*Ventilation.pulmocycle];
switch Collapsion_areas.collapsions.collapsion_right_lung
    case 1
        RL.position(12:15,1)=1/2*(RL.shape.empty(12:15,1)+RL.shape.full(12:15,1));
        RL.position(12:15,2)=1/2*(RL.shape.empty(12:15,2)+RL.shape.full(12:15,2));
        RL.position(1:2,1)=1/2*(RL.shape.empty(1:2,1)+RL.shape.full(1:2,1));
        RL.position(1:2,2)=1/2*(RL.shape.empty(1:2,2)+RL.shape.full(1:2,2));
    case 2
        RL.position(14:15,1)=1/2*(RL.shape.empty(14:15,1)+RL.shape.full(14:15,1));
        RL.position(14:15,2)=1/2*(RL.shape.empty(14:15,2)+RL.shape.full(14:15,2));
        RL.position(1:2,1)=1/2*(RL.shape.empty(1:2,1)+RL.shape.full(1:2,1));
        RL.position(1:2,2)=1/2*(RL.shape.empty(1:2,2)+RL.shape.full(1:2,2));
    case 3
        RL.position(14:15,1)=1/2*(RL.shape.empty(14:15,1)+RL.shape.full(14:15,1));
        RL.position(14:15,2)=1/2*(RL.shape.empty(14:15,2)+RL.shape.full(14:15,2));
        RL.position(1,1:2)=1/2*(RL.shape.empty(1,1:2)+RL.shape.full(1,1:2));
    otherwise
        %do nothing
end
%Substitution
Tissues.Right_lung=RL;

H=Tissues.Heart.exterior;
H.position=[H.shape.full(:,1)+H.shape.step(:,1)*Circulation.cardio_cycle,H.shape.full(:,2)+H.shape.step(:,2)*Circulation.cardio_cycle];
Tissues.Heart.exterior=H;
Hint=Tissues.Heart.interior;
Hint.position=[Hint.shape.full(:,1)+Hint.shape.step(:,1)*Circulation.cardio_cycle,Hint.shape.full(:,2)+Hint.shape.step(:,2)*Circulation.cardio_cycle];
Tissues.Heart.interior=Hint;

B=Tissues.Bones;
for bone=1:8
    B_k=eval(['B.Bone',num2str(bone)]);
    B_k.position=[B_k.shape.empty(:,1)-B_k.shape.step(:,1)*Ventilation.pulmocycle,B_k.shape.empty(:,2)-B_k.shape.step(:,2)*Ventilation.pulmocycle];
    switch bone
        case 1
            B.Bone1=B_k;
        case 2
            B.Bone2=B_k;
        case 3
            B.Bone3=B_k;
        case 4
            B.Bone4=B_k;
        case 5
            B.Bone5=B_k;
        case 6
            B.Bone6=B_k;
        case 7
            B.Bone7=B_k;
        case 8
            B.Bone8=B_k;
    end
    %eval(['B.Bone',num2str(k)])=B_k;
end
Tissues.Bones=B;
end


function [assist_points,electrode_central_points,electrode_exterior_sides]=dynamic_geometry(N,electrodes,Tissues,Collapsion_areas,fileID)
%Takes the current step electrode and Tissue positions
%Generates the LUA code needed to design them in FEMM
%Returns the electode assist points from where the measurements are going
%to be taken later
LL=Tissues.Left_lung;
RL=Tissues.Right_lung;
H=Tissues.Heart;

if N==16 %only if 16 electrodes are used, arcs are needed
    %otherwise they only increase the mesh nodes without important reason
    angles1=[1 40 10 20 10 20 10 20 10 26 5 26 5 22.5 1 22.5 1];
    len=length(angles1);
    angles=[angles1 fliplr(angles1(1:len-1))];
end
d=electrodes.delec;
electrode_position=electrodes.position;
len=length(electrode_position.x);
% figure
% plot(electrode_position.x,electrode_position.y,'o')
% hold on
%add the electrode points
%add the electrode points
electrode_central_points=zeros(len,2);
electrode_exterior_sides=zeros(len,4);
extarc=5;
for i=1:len
    fprintf(fileID,'ci_addnode(%4.2f,%4.2f)\n',electrode_position.x(1,i),electrode_position.y(1,i));
    fprintf(fileID,'ci_addnode(%4.2f,%4.2f)\n',electrode_position.x(2,i),electrode_position.y(2,i));
    %%%%
    %add the exterior electrode volume points
    %d=0.5;
    x1=electrode_position.x(1,i);
    y1=electrode_position.y(1,i);
    x2=electrode_position.x(2,i);
    y2=electrode_position.y(2,i);
    a=-d*abs(x2-x1)*(y2-y1)/((x2-x1)*sqrt((x2-x1)^2+(y2-y1)^2));
    b=d*abs(x2-x1)/sqrt((x2-x1)^2+(y2-y1)^2);
    if y2>-3
        x3=x1+a; y3=y1+b;
        x4=x2+a; y4=y2+b;
    else
        x3=x1-a; y3=y1-b;
        x4=x2-a; y4=y2-b;
    end
    electrode_central_points(i,:)=[(x1+x4)/2 (y1+y4)/2];
    electrode_exterior_sides(i,:)=[x3 y3 x4 y4];
    %     plot(x3,y3,'+')
    %     hold on
    %     plot(x4,y4,'+')
    %     hold on
    fprintf(fileID,'ci_addnode(%4.2f,%4.2f)\n',x3,y3);
    fprintf(fileID,'ci_addnode(%4.2f,%4.2f)\n',x4,y4);
    %connect (x1,y1) with (x3,y3), (x3,y3) with (x4,y4) and (x2,y2) with
    %(x4,y4) (give dimensions to the electrodes)
    fprintf(fileID,'ci_addsegment(%4.2f,%4.2f,%4.2f,%4.2f)\n',x1,y1,x3,y3);
    fprintf(fileID,'ci_addarc(%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,1)\n',x4,y4,x3,y3,extarc);
    fprintf(fileID,'ci_addsegment(%4.2f,%4.2f,%4.2f,%4.2f)\n',x2,y2,x4,y4);
    %%%%
end
%fprintf('%4.4f,%4.4f\n',electrode_central_points(8,1),electrode_central_points(8,2))

%connect the electrodes with arcs if N=16 or straight if N=32 or 64
for i=1:1:len
    inext=i+1;
    if i==len
        inext=1;
    end
    if N==16
        fprintf(fileID,'ci_addarc(%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,1)\n',electrode_position.x(1,inext),electrode_position.y(1,inext),electrode_position.x(2,i),electrode_position.y(2,i),angles(2*i));
    else
        fprintf(fileID,'ci_addsegment(%4.2f,%4.2f,%4.2f,%4.2f)\n',electrode_position.x(1,inext),electrode_position.y(1,inext),electrode_position.x(2,i),electrode_position.y(2,i));
    end
end
%connect the 2 points of each electrode with arcs
if N==16
    for i=1:1:len
        fprintf(fileID,'ci_addarc(%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,1)\n',electrode_position.x(2,i),electrode_position.y(2,i),electrode_position.x(1,i),electrode_position.y(1,i),angles(2*i-1));
    end
else
    for i=1:1:len
        fprintf(fileID,'ci_addarc(%4.2f,%4.2f,%4.2f,%4.2f,%4.2f,1)\n',electrode_position.x(2,i),electrode_position.y(2,i),electrode_position.x(1,i),electrode_position.y(1,i),10);
    end
end

%add the left_lung points and connect them
tissue_marker(fileID,LL.position);
%now add the segment that separates the main lung from the collapse area
%%add the LEFT lung's collapsed area
collapsed_boundary_maker(fileID,LL.position,Collapsion_areas.cross_points.left_lung);

%add the right_lung points and connect them
%add the left_lung points and connect them
tissue_marker(fileID,RL.position);
%%add the right lung's collapsed area
collapsed_boundary_maker(fileID,RL.position,Collapsion_areas.cross_points.right_lung);

%add the Heart's exterior points and connect them
tissue_marker(fileID,H.exterior.position);

%add the Heart's interior points and connect them
tissue_marker(fileID,H.interior.position);

%add each Bone's points and connect them
for k=1:8
    B_k=eval(['Tissues.Bones.Bone' num2str(k)]);
    tissue_marker(fileID,B_k.position);
end

%connect fat boundary points
B1p=Tissues.Bones.Bone1.position;
B2p=Tissues.Bones.Bone2.position;
%from 3rd point of 1st Bone to 1st point of 2nd Bone
fprintf(fileID,'ci_addsegment(%4.2f,%4.2f,%4.2f,%4.2f)\n',B1p(3,1),B1p(3,2),B2p(1,1),B2p(1,2));
B3p=Tissues.Bones.Bone3.position;
fprintf(fileID,'ci_addsegment(%4.2f,%4.2f,%4.2f,%4.2f)\n',B2p(4,1),B2p(4,2),B3p(1,1),B3p(1,2));
B4p=Tissues.Bones.Bone4.position;
fprintf(fileID,'ci_addsegment(%4.2f,%4.2f,%4.2f,%4.2f)\n',B3p(9,1),B3p(9,2),B4p(1,1),B4p(1,2));
B5p=Tissues.Bones.Bone5.position;
fprintf(fileID,'ci_addsegment(%4.2f,%4.2f,%4.2f,%4.2f)\n',B4p(4,1),B4p(4,2),B5p(2,1),B5p(2,2));
B6p=Tissues.Bones.Bone6.position;
fprintf(fileID,'ci_addsegment(%4.2f,%4.2f,%4.2f,%4.2f)\n',B5p(1,1),B5p(1,2),B6p(1,1),B6p(1,2));
B7p=Tissues.Bones.Bone7.position;
fprintf(fileID,'ci_addsegment(%4.2f,%4.2f,%4.2f,%4.2f)\n',B6p(4,1),B6p(4,2),B7p(1,1),B7p(1,2));
B8p=Tissues.Bones.Bone8.position;
fprintf(fileID,'ci_addsegment(%4.2f,%4.2f,%4.2f,%4.2f)\n',B7p(4,1),B7p(4,2),B8p(1,1),B8p(1,2));
fprintf(fileID,'ci_addsegment(%4.2f,%4.2f,%4.2f,%4.2f)\n',B8p(4,1),B8p(4,2),B1p(4,1),B1p(4,2));

%derive assist points
len=length(electrode_position.x);
assist_points=zeros(len,2);
for i=1:len
    mid1=(electrode_position.x(1,i)+electrode_position.x(2,i))/2;
    mid2=(electrode_position.y(1,i)+electrode_position.y(2,i))/2;
    assist_points(i,:)=[mid1 mid2];
end
end

function []=tissue_marker(fileID,Positions)
% takes the specific tissue position points at the current time step
% writes the LUA code for the specific tissue design in FEMM
len=length(Positions);
for i=1:len
    fprintf(fileID,'ci_addnode(%4.2f,%4.2f)\n',Positions(i,1),Positions(i,2));
    if i>1
        fprintf(fileID,'ci_addsegment(%4.2f,%4.2f,%4.2f,%4.2f)\n',Positions(i-1,1),Positions(i-1,2),Positions(i,1),Positions(i,2));
    end
    if i==len
        fprintf(fileID,'ci_addsegment(%4.2f,%4.2f,%4.2f,%4.2f)\n',Positions(i,1),Positions(i,2),Positions(1,1),Positions(1,2));
    end
end
end

function []=collapsed_boundary_maker(fileID,Position,Cross_points)
% takes the (R or L) lung's tissue position points at the current time step and
% the collapsion area's cross points
% writes the LUA code for the specific collapse area boundary design in FEMM
fprintf(fileID,'ci_addsegment(%4.2f,%4.2f,%4.2f,%4.2f)\n',Position(Cross_points(1),1),...
    Position(Cross_points(1),2),Position(Cross_points(2),1),Position(Cross_points(2),2));
end

function[LLCC,RLCC,CALCC,CARCC,STCC,Xbone,Ybone,Xfat,Yfat,Xint,Yint,Xext,Yext]=...
    addmaterial_labels(Tissues,Ventilation,Circulation,Collapsion_areas,fileID,r,k,electrode_central_points)
%adds material labels in a point that is in the corresponding tissue area
%(in FEMM)
%each tissue dynamic shape centroid is selected for the corresponding
%material label
%Left Lung Properties
[LLccx,LLccy]=find_centroid(Tissues.Left_lung.position);
fprintf(fileID,'ci_addblocklabel(%4.2f,%4.2f)\n',LLccx,LLccy);
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',LLccx, LLccy);
fprintf(fileID,'ci_setblockprop("Lungs_state_%d_%d",1,"<None>","<None>")\n',Ventilation.pulmocycle,Circulation.cardio_cycle);
fprintf(fileID,'ci_clearselected()\n');
%set left lung collapsion area
if Collapsion_areas.collapsions.collapsion_left_lung~=0
    CACCL=Collapsion_areas.cross_points.left_lung;
    CAL=[Tissues.Left_lung.position(CACCL(2):1:length(Tissues.Left_lung.position(:,1)),:); Tissues.Left_lung.position(1:CACCL(1),:)];
    [CALccx,CALccy]=find_centroid(CAL);
    CALCC=[CALccx CALccy];
    fprintf(fileID,'ci_addblocklabel(%4.2f,%4.2f)\n',CALccx,CALccy);
    fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',CALccx,CALccy);
    fprintf(fileID,'ci_setblockprop("Lungs",1,"<None>","<None>")\n');
    fprintf(fileID,'ci_clearselected()\n');
else
    CALCC=[];
end
%Right Lung Properties
[RLccx,RLccy]=find_centroid(Tissues.Right_lung.position);
fprintf(fileID,'ci_addblocklabel(%4.2f,%4.2f)\n',RLccx,RLccy);
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',RLccx,RLccy);
fprintf(fileID,'ci_setblockprop("Lungs_state_%d_%d",1,"<None>","<None>")\n',Ventilation.pulmocycle,Circulation.cardio_cycle);
fprintf(fileID,'ci_clearselected()\n');
if Collapsion_areas.collapsions.collapsion_right_lung~=0
    %set right lung collapsion area
    CACCR=Collapsion_areas.cross_points.right_lung;
    CAR=[Tissues.Right_lung.position(CACCR(2):1:length(Tissues.Right_lung.position(:,1)),:); Tissues.Right_lung.position(1:CACCR(1),:)];
    [CARccx,CARccy]=find_centroid(CAR);
    CARCC=[CARccx CARccy];
    fprintf(fileID,'ci_addblocklabel(%4.2f,%4.2f)\n',CARccx,CARccy);
    fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',CARccx,CARccy);
    fprintf(fileID,'ci_setblockprop("Lungs",1,"<None>","<None>")\n');
    fprintf(fileID,'ci_clearselected()\n');
else
    CARCC=[];
end

%Internal Heart Properties
IntPos=Tissues.Heart.interior.position;
Xint=(IntPos(1,1)+IntPos(5,1))/2;
Yint=(IntPos(1,2)+IntPos(5,2))/2;

fprintf(fileID,'ci_addblocklabel(%4.2f,%4.2f)\n',Xint,Yint);
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xint,Yint);
fprintf(fileID,'ci_setblockprop("Heart_state_%d",1,"<None>","<None>")\n',Circulation.cardio_cycle);
fprintf(fileID,'ci_clearselected()\n');

%External Heart Properties
ExPose=Tissues.Heart.exterior.position;
Xext=(ExPose(4,1)+ExPose(9,1))/2;
Yext=(ExPose(4,2)+ExPose(9,2))/2;

fprintf(fileID,'ci_addblocklabel(%4.2f,%4.2f)\n',Xext,Yext);
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xext,Yext);
fprintf(fileID,'ci_setblockprop("Miocardium",1,"<None>","<None>")\n');
fprintf(fileID,'ci_clearselected()\n');


%Bone properties
%Bone 7 (down)
Bonepos7=Tissues.Bones.Bone7.position;
[Xbone7,Ybone7]=find_centroid(Bonepos7);
fprintf(fileID,'ci_addblocklabel(%4.2f,%4.2f)\n',Xbone7,Ybone7);
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xbone7,Ybone7);
fprintf(fileID,'ci_setblockprop("Bone",1,"<None>","<None>")\n');
fprintf(fileID,'ci_clearselected()\n');

%Bone 6
Bonepos6=Tissues.Bones.Bone6.position;
[Xbone6,Ybone6]=find_centroid(Bonepos6);
fprintf(fileID,'ci_addblocklabel(%4.2f,%4.2f)\n',Xbone6,Ybone6);
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xbone6,Ybone6);
fprintf(fileID,'ci_setblockprop("Bone",1,"<None>","<None>")\n');
fprintf(fileID,'ci_clearselected()\n');

%Bone 5
Bonepos5=Tissues.Bones.Bone5.position;
[Xbone5,Ybone5]=find_centroid(Bonepos5);
fprintf(fileID,'ci_addblocklabel(%4.2f,%4.2f)\n',Xbone5,Ybone5);
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xbone5,Ybone5);
fprintf(fileID,'ci_setblockprop("Bone",1,"<None>","<None>")\n');
fprintf(fileID,'ci_clearselected()\n');

%Bone 4
Bonepos4=Tissues.Bones.Bone4.position;
[Xbone4,Ybone4]=find_centroid(Bonepos4);
fprintf(fileID,'ci_addblocklabel(%4.2f,%4.2f)\n',Xbone4,Ybone4);
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xbone4,Ybone4);
fprintf(fileID,'ci_setblockprop("Bone",1,"<None>","<None>")\n');
fprintf(fileID,'ci_clearselected()\n');

%Bone 3
Bonepos3=Tissues.Bones.Bone3.position;
[Xbone3,Ybone3]=find_centroid(Bonepos3);
fprintf(fileID,'ci_addblocklabel(%4.2f,%4.2f)\n',Xbone3,Ybone3);
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xbone3,Ybone3);
fprintf(fileID,'ci_setblockprop("Bone",1,"<None>","<None>")\n');
fprintf(fileID,'ci_clearselected()\n');

%Bone 2
Bonepos2=Tissues.Bones.Bone2.position;
[Xbone2,Ybone2]=find_centroid(Bonepos2);
fprintf(fileID,'ci_addblocklabel(%4.2f,%4.2f)\n',Xbone2,Ybone2);
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xbone2,Ybone2);
fprintf(fileID,'ci_setblockprop("Bone",1,"<None>","<None>")\n');
fprintf(fileID,'ci_clearselected()\n');

%Bone 1
Bonepos1=Tissues.Bones.Bone1.position;
[Xbone1,Ybone1]=find_centroid(Bonepos1);
fprintf(fileID,'ci_addblocklabel(%4.2f,%4.2f)\n',Xbone1,Ybone1);
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xbone1,Ybone1);
fprintf(fileID,'ci_setblockprop("Bone",1,"<None>","<None>")\n');
fprintf(fileID,'ci_clearselected()\n');

%Bone 8
Bonepos8=Tissues.Bones.Bone8.position;
[Xbone8,Ybone8]=find_centroid(Bonepos8);
fprintf(fileID,'ci_addblocklabel(%4.2f,%4.2f)\n',Xbone8,Ybone8);
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xbone8,Ybone8);
fprintf(fileID,'ci_setblockprop("Bone",1,"<None>","<None>")\n');
fprintf(fileID,'ci_clearselected()\n');

Xbone=[Xbone1 Xbone2 Xbone3 Xbone4 Xbone5 Xbone6 Xbone7 Xbone8];
Ybone=[Ybone1 Ybone2 Ybone3 Ybone4 Ybone5 Ybone6 Ybone7 Ybone8];

% STccx=-k*6.7;
% STccy=-k*15.5;
STccx=(Xbone(6)+Tissues.Left_lung.position(15,1))/2;
STccy=(Ybone(6)+Tissues.Left_lung.position(15,2))/2;
fprintf(fileID,'ci_addblocklabel(%4.2f,%4.2f)\n',STccx,STccy);
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',STccx,STccy);
fprintf(fileID,'ci_setblockprop("Soft_tissue_state_%d",1,"<None>","<None>")\n',Circulation.cardio_cycle);
fprintf(fileID,'ci_clearselected()\n');
STCC=[STccx STccy];

Xfat=(Bonepos6(3,1)+Bonepos7(2,1))/2;
Yfat=(Bonepos6(3,2)+Bonepos7(2,2))/2;
%fat properties
fprintf(fileID,'ci_addblocklabel(%4.2f,%4.2f)\n',Xfat,Yfat);
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xfat,Yfat);
fprintf(fileID,'ci_setblockprop("Fat",1,"<None>","<None>")\n');
fprintf(fileID,'ci_clearselected()\n');

%%%electrodes
for electrode=1:length(electrode_central_points(:,1))
    fprintf(fileID,'ci_addblocklabel(%4.2f,%4.2f)\n',...
        electrode_central_points(electrode,1),electrode_central_points(electrode,2));
    fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',...
        electrode_central_points(electrode,1),electrode_central_points(electrode,2));
    fprintf(fileID,'ci_setblockprop("Zelectrode%2.0f",1,"<None>","<None>")\n',electrode);
    fprintf(fileID,'ci_clearselected()\n');
end

LLCC=[LLccx LLccy];
RLCC=[RLccx RLccy];

end

function [ccx,ccy]=find_centroid(Shape)
polyShape=polyshape([Shape(:,1); Shape(1,1)],[Shape(:,2); Shape(1,2)]);
[ccx,ccy]=centroid(polyShape);
end

function set_current_source(N,skipcurr,assist_points,count,fileID,electrode_exterior_sides)
k=skipcurr+1;
%set input current electrode
% fprintf(fileID,'ci_selectarcsegment(%4.2f,%4.2f)\n',assist_points(count,1),assist_points(count,2));
% fprintf(fileID,'ci_setarcsegmentprop(4,"<None>",0,"<None>","Iin")\n');
% fprintf(fileID,'ci_clearselected()\n');
fprintf(fileID,'ci_selectarcsegment(%4.2f,%4.2f)\n',(electrode_exterior_sides(count,1)+electrode_exterior_sides(count,3))/2,...
    (electrode_exterior_sides(count,2)+electrode_exterior_sides(count,4))/2);
fprintf(fileID,'ci_setarcsegmentprop(4,"<None>",0,"<None>","Iin")\n');
fprintf(fileID,'ci_clearselected()\n');
% fprintf(fileID,'ci_selectsegment(%4.2f,%4.2f)\n',(electrode_exterior_sides(count,1)+electrode_exterior_sides(count,3))/2,...
%     (electrode_exterior_sides(count,2)+electrode_exterior_sides(count,4))/2);
% fprintf(fileID,'ci_setsegmentprop(4,"<None>",0,"<None>","Iin")\n');
% fprintf(fileID,'ci_clearselected()\n');
if count+k<=N
    next=count+k;
else
    next=count+k-N;
end
%set output current electrode
% fprintf(fileID,'ci_selectarcsegment(%4.2f,%4.2f)\n',assist_points(next,1),assist_points(next,2));
% fprintf(fileID,'ci_setarcsegmentprop(4,"<None>",0,"<None>","Iout")\n');
% fprintf(fileID,'ci_clearselected()\n');
fprintf(fileID,'ci_selectarcsegment(%4.2f,%4.2f)\n',(electrode_exterior_sides(next,1)+electrode_exterior_sides(next,3))/2,...
    (electrode_exterior_sides(next,2)+electrode_exterior_sides(next,4))/2);
fprintf(fileID,'ci_setarcsegmentprop(4,"<None>",0,"<None>","Iin")\n');
fprintf(fileID,'ci_clearselected()\n');
% fprintf(fileID,'ci_selectsegment(%4.2f,%4.2f)\n',(electrode_exterior_sides(next,1)+electrode_exterior_sides(next,3))/2,...
%     (electrode_exterior_sides(next,2)+electrode_exterior_sides(next,4))/2);
% fprintf(fileID,'ci_setsegmentprop(4,"<None>",0,"<None>","Iout")\n');
% fprintf(fileID,'ci_clearselected()\n');
end

function set_reference(N,skipcurr,assist_points,count,fileID,electrode_exterior_sides)
k=skipcurr+1;
%set a ground reference
if count+k+1<=N
    ref=count+k+1;
else
    ref=count+k+1-N;
end
% fprintf(fileID,'ci_selectarcsegment(%4.2f,%4.2f)\n',assist_points(ref,1),assist_points(ref,2));
% fprintf(fileID,'ci_setarcsegmentprop(4,"Ground",0,"<None>","<None>")\n');
% fprintf(fileID,'ci_clearselected()\n');
fprintf(fileID,'ci_selectarcsegment(%4.2f,%4.2f)\n',(electrode_exterior_sides(ref,1)+electrode_exterior_sides(ref,3))/2,...
    (electrode_exterior_sides(ref,2)+electrode_exterior_sides(ref,4))/2);
fprintf(fileID,'ci_setarcsegmentprop(4,"Ground",0,"<None>","<None>")\n');
fprintf(fileID,'ci_clearselected()\n');
% fprintf(fileID,'ci_selectsegment(%4.2f,%4.2f)\n',(electrode_exterior_sides(ref,1)+electrode_exterior_sides(ref,3))/2,...
%     (electrode_exterior_sides(ref,2)+electrode_exterior_sides(ref,4))/2);
% fprintf(fileID,'ci_setsegmentprop(4,"Ground",0,"<None>","<None>")\n');
% fprintf(fileID,'ci_clearselected()\n');
end

function []=save_process(d,fileID,folds)
fprintf(fileID,'ci_saveas("%s%s_test%.0f.fee")\n',folds,date(),d);
fprintf(fileID,'ci_createmesh()\n');
fprintf(fileID,'ci_saveas("%s%s_test%.0f.fee")\n',folds,date(),d);
fprintf(fileID,'ci_saveas("%s%s_test%.0f.fec")\n',folds,date(),d);
fprintf(fileID,'ci_analyze(0)\n');
fprintf(fileID,'ci_loadsolution()\n');
end

function measure_voltages(N,electrode_position,d,timestep,fileID,electrode_exterior_sides,folds)
%writes LUA commands in order to measure the electrode voltages
o=[];
for i=1:N
    fprintf(fileID,'co_addcontour(%4.2f,%4.2f)\n',electrode_exterior_sides(i,1),electrode_exterior_sides(i,2));
    fprintf(fileID,'co_addcontour(%4.2f,%4.2f)\n',electrode_exterior_sides(i,3),electrode_exterior_sides(i,4));
    fprintf(fileID,'co_makeplot(0,10,0)\n');
    fprintf(fileID,'co_makeplot(0,%.0f,"%s%s_test%.0f_Voltages%d_%d.txt",0)\n',N,folds,date(),d,timestep,i);
    fprintf(fileID,'co_clearcontour()\n');
end
fprintf(fileID,'co_close()\n');
fprintf(fileID,'ci_purgemesh()\n');
pause(0.1)
end

function [Ventilation,Circulation,total_positions,times,Silences,Tissues,electrodes,r,k]=...
    update_cycle_states(N,Ventilation,Circulation,...
    total_positions,times,Silences,Tissues,time,electrodes,r,k,frequency,end_Admittances,delec,randomcond)
%updates breath and cardiac cycle indicators
%updates breath and cardiac cycle durations and tissues discrete positions
%after the end of a correspoding cycle
%called in each for loop iteration

Ventilation.pulmocycle=Ventilation.pulmocycle+Ventilation.pulmcounter;
Circulation.cardio_cycle=Circulation.cardio_cycle+Circulation.cardcounter;

% plot(time,Ventilation.pulmocycle,'b*')
% hold on
% plot(time,Circulation.cardio_cycle,'r*')
% hold on
% if pulmonary cycle reaches left/right limit, change the indicator
if Ventilation.pulmocycle>total_positions.total_positions-1
    Ventilation.pulmcounter=-Ventilation.pulmcounter;
end
% if pulmonary cycle is at 1 and breath pausing has not started
% start the pause duration
if (Ventilation.pulmocycle<2)&&(Silences.silence==0)
    %relaxing time between 0 and 1 seconds
    if randomcond==1
        Silences.relaxing_time=rand(1);
    else
        Silences.relaxing_time=0.3;
    end
    %pause breathing
    Silences.silence=1;
    %store the time point pause started
    Silences.time_sil_start=time;
    % zero the indicator during relaxing time
    Ventilation.pulmcounter=0;
end
% if relaxing time has come to an end
% stop pausing and redefine breathig time
if (Silences.time_sil_start+Silences.relaxing_time<=time)&&(Silences.silence==1)
    %breath period will change up to 20% at upcoming cycle
    if randomcond==1
        times.tbr=times.tbr+0.2*2*(round(rand(1,1))-1)*rand(1)*times.tbr;
    else
        times.tbr=1.2*times.tbr;
    end
    %%positive pulmonary cycle indicator
    %to start the upcoming breath cycle
    Ventilation.pulmcounter=1;
    %stop pausing
    Silences.silence=0;
    %no more relax!
    Silences.time_sil_start=-1;
    Silences.relaxing_time=-1;
    %new position step for tissues
    total_positions.dy=2*1*times.tim/(N*times.tbr); %DY=1
    %which changes total possible positions
    %not their boundaries
    %but their possible positions!
    total_positions.total_positions=1/total_positions.dy+1; %Dy=1
    %this either changes the possible tissue conductivities states!!
    %and either the tissue possible positions!!
    %(lungs and front bones)
    [Tissues.Left_lung,Tissues.Right_lung,Tissues.Heart,Tissues.Muscles,Tissues.Bones,r,k]=...
        define_tissues(r,k,times,total_positions,1,frequency,end_Admittances,randomcond,Tissues);
    fprintf('Positions have been redefined!\n');
    %either the possible (but not the limits) of the electrode
    %positions!
    [electrodes.electrode_average,electrodes.electrode_possible]=electrodes_preset(N,1,total_positions.dy,delec,k); %Dy=1
    %and either the tissue possible positions!!
    %(lungs and front bones)
end
%if cardio cycle reaches left/right limit
%change the indicator
if Circulation.cardio_cycle>times.pulse_dur-1
    Circulation.cardcounter=-Circulation.cardcounter;
end
if Circulation.cardio_cycle<2
    %if left limit reached
    %change the indicator
    Circulation.cardcounter=-Circulation.cardcounter;
    %and add one cardio beat!
    Silences.cardiobeats=Silences.cardiobeats+1;
    fpulse=1/times.Tpulse;
    if mod(Silences.cardiobeats,3)==0 %heart rythm change up to 3% each 3 beats
        %define new heart rythm every start after 3 cardio beats
        if randomcond==1
            fpulse=fpulse+0.03*2*(round(rand(1,1))-1)*rand(1)*fpulse;
        else
            fpulse=fpulse*1.03;
        end
        times.Tpulse=1/fpulse;
        %redefining the cardio frequency
        %ALSO leads to redefine the possible tissue conductivity
        %states!
        %and positions! (of the heart)
        %but not any change at limits!
        %increasing frequency leads to less positions
        %they have less time to come across the limits
        %thus they have to "jump" more
        %this increases the step!
        %and decreases the kwantized positions
        fprintf('Cardio has been redefined!\n');
        times.pulse_dur=2*total_positions.total_positions*(times.tbr/times.Tpulse)^(-1);
        [Tissues.Left_lung,Tissues.Right_lung,Tissues.Heart,Tissues.Muscles,Tissues.Bones,r,k]=...
        define_tissues(r,k,times,total_positions,1,frequency,end_Admittances,randomcond,Tissues);
    end
end
end

function []=unset_current_source(N,fileID,k,assist_points,count,electrode_exterior_sides)
%unset input current electrode
% fprintf(fileID,'ci_selectarcsegment(%4.2f,%4.2f)\n',assist_points(count,1),assist_points(count,2));
% fprintf(fileID,'ci_setarcsegmentprop(4,"<None>",0,"<None>","<None>")\n');
% fprintf(fileID,'ci_clearselected()\n');
fprintf(fileID,'ci_selectarcsegment(%4.2f,%4.2f)\n',(electrode_exterior_sides(count,1)+electrode_exterior_sides(count,3))/2,...
    (electrode_exterior_sides(count,2)+electrode_exterior_sides(count,4))/2);
fprintf(fileID,'ci_setarcsegmentprop(4,"<None>",0,"<None>","None")\n');
fprintf(fileID,'ci_clearselected()\n');
if count+k<=N
    next=count+k;
else
    next=count+k-N;
end
%unset output current electrode
% fprintf(fileID,'ci_selectarcsegment(%4.2f,%4.2f)\n',assist_points(next,1),assist_points(next,2));
% fprintf(fileID,'ci_setarcsegmentprop(4,"<None>",0,"<None>","<None>")\n');
% fprintf(fileID,'ci_clearselected()\n');
fprintf(fileID,'ci_selectarcsegment(%4.2f,%4.2f)\n',(electrode_exterior_sides(next,1)+electrode_exterior_sides(next,3))/2,...
    (electrode_exterior_sides(next,2)+electrode_exterior_sides(next,4))/2);
fprintf(fileID,'ci_setarcsegmentprop(4,"<None>",0,"<None>","None")\n');
fprintf(fileID,'ci_clearselected()\n');
end

function unset_reference(N,fileID,k,assist_points,count,electrode_exterior_sides)
%unset  ground reference
if count+k+1<=N
    ref=count+k+1;
else
    ref=count+k+1-N;
end
% fprintf(fileID,'ci_selectarcsegment(%4.2f,%4.2f)\n',assist_points(ref,1),assist_points(ref,2));
% fprintf(fileID,'ci_setarcsegmentprop(4,"<None>",0,"<None>","<None>")\n');
% fprintf(fileID,'ci_clearselected()\n');
fprintf(fileID,'ci_selectarcsegment(%4.2f,%4.2f)\n',(electrode_exterior_sides(ref,1)+electrode_exterior_sides(ref,3))/2,...
    (electrode_exterior_sides(ref,2)+electrode_exterior_sides(ref,4))/2);
fprintf(fileID,'ci_setarcsegmentprop(4,"<None>",0,"<None>","None")\n');
fprintf(fileID,'ci_clearselected()\n');
end

function []=undo_dynamic_geometry(Tissues,fileID,electrode_position,electrode_exterior_sides)
len=length(electrode_position.x);
LL=Tissues.Left_lung.position;
RL=Tissues.Right_lung.position;
H=Tissues.Heart.exterior.position;
Hint=Tissues.Heart.interior.position;
B=Tissues.Bones;
%select the electrode points
for i=1:len
    fprintf(fileID,'ci_selectnode(%4.2f,%4.2f)\n',electrode_position.x(1,i),electrode_position.y(1,i));
    fprintf(fileID,'ci_selectnode(%4.2f,%4.2f)\n',electrode_position.x(2,i),electrode_position.y(2,i));
end
%select the left_lung points
len=length(LL);
for i=1:len
    fprintf(fileID,'ci_selectnode(%4.2f,%4.2f)\n',LL(i,1),LL(i,2));
end
%select the right_lung points
len=length(RL);
for i=1:len
    fprintf(fileID,'ci_selectnode(%4.2f,%4.2f)\n',RL(i,1),RL(i,2));
end
%select the Heart's points
len=length(H);
for i=1:len
    fprintf(fileID,'ci_selectnode(%4.2f,%4.2f)\n',H(i,1),H(i,2));
end
len=length(Hint);
for i=1:len
    fprintf(fileID,'ci_selectnode(%4.2f,%4.2f)\n',Hint(i,1),Hint(i,2));
end
%select each Bone's points
for k=1:8
    len=length(eval(['B.Bone' num2str(k) '.position(:,1)']));
    for i=1:len
        fprintf(fileID,'ci_selectnode(%4.2f,%4.2f)\n',eval(['B.Bone' num2str(k) '.position(i,1)']),eval(['B.Bone' num2str(k) '.position(i,2)']));
    end
end
%%%select electrode external points
for electrode=1:length(electrode_exterior_sides(:,1))
    fprintf(fileID,'ci_selectnode(%4.2f,%4.2f)\n',...
        electrode_exterior_sides(electrode,1),electrode_exterior_sides(electrode,2));
    fprintf(fileID,'ci_selectnode(%4.2f,%4.2f)\n',...
        electrode_exterior_sides(electrode,3),electrode_exterior_sides(electrode,4));
end
%delete all selected nodes (with their segments and arcs)
fprintf(fileID,'ci_deleteselectednodes()\n');
end

function []=clear_labels(fileID,LLCC,RLCC,CALCC,CARCC,STCC,Xbone,Ybone,...
    Xfat,Yfat,Xint,Yint,Xext,Yext,electrode_central_points)
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',LLCC(1), LLCC(2));
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',RLCC(1), RLCC(2));
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xint,Yint);
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xext,Yext);
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xbone(1),Ybone(1));
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xbone(2),Ybone(2));
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xbone(3),Ybone(3));
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xbone(4),Ybone(4));
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xbone(5),Ybone(5));
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xbone(6),Ybone(6));
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xbone(7),Ybone(7));
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xbone(8),Ybone(8));
if ~isempty(CALCC)
    fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',CALCC(1),CALCC(2));
end
if ~isempty(CARCC)
    fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',CARCC(1),CARCC(2));
end
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',STCC(1),STCC(2));
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',-6,-11.85);
fprintf(fileID,'ci_selectlabel(%4.2f,%4.2f)\n',Xfat,Yfat);
for electrode=1:length(electrode_central_points(:,1))
    fprintf(fileID,'ci_selectlabel(%4.4f,%4.4f)\n',...
        electrode_central_points(electrode,1),electrode_central_points(electrode,2));
end
%fprintf('%4.4f,%4.4f\n',electrode_central_points(8,1),electrode_central_points(8,2))
fprintf(fileID,'ci_deleteselectedlabels()\n');
end

function [letterplus,Gr]=configure_graph(Graph,frames,fps,total_time_est,path2save)
[Gr,time]=plot_graph(Graph,frames,fps,total_time_est);
letterplus=save_graph(Graph,time,path2save);
fprintf('Please open a FEMM project!\n')
fprintf('Then load the corresponding LUA file and press enter\n')
pause()
end

function [Gr,time]=plot_graph(Graph,frames,fps,total_time_est)
Gr=figure(1);
subplot(2,1,1)
time=linspace(0, total_time_est,length(Graph.heart.conductivity));
plot(time,Graph.heart.conductivity,'r','LineWidth',2)
hold on
Graph.right_lung.conductivity=Graph.left_lung.conductivity;
plot(time,Graph.right_lung.conductivity,'b','LineWidth',2)
hold on
plot(time,Graph.left_lung.conductivity,'g--','LineWidth',2)
hold on
plot(time,Graph.muscles.conductivity(1:end),'y','LineWidth',2)
hold on
TGr=title('Tissues conductivity (S/m)');
TGr.FontSize=16;
legend('\sigma_{heart}','\sigma_{right lung}','\sigma_{left lung}','\sigma_{muscles}')
ylabel('\sigma (Sm^-^1)','Fontsize',14)
xlabel('Time (seconds)','Fontsize',14)
subplot(2,1,2)
plot(time,Graph.heart.permittivity,'r','LineWidth',2)
hold on
Graph.right_lung.permittivity=Graph.left_lung.permittivity;
plot(time,Graph.right_lung.permittivity,'b','LineWidth',2)
hold on
plot(time,Graph.left_lung.permittivity,'b','LineWidth',2)
hold on
plot(time,Graph.muscles.permittivity(1:end),'y','LineWidth',2)
TGr=title('Tissues permittivity (F/m)');
TGr.FontSize=16;
legend('\epsilon_{heart}','\epsilon_{right lung}','\epsilon_{left lung}','\epsilon_{muscles}')
ylabel('\epsilon (Fm^-^1)','Fontsize',14)
xlabel('Time (seconds)','Fontsize',14)
end


function [letterplus]=save_graph(Graph,time,path2save)
letterplus='a';
for letter='a':'z'
    filetosearch=[path2save 'Cardio-Pulmonary\Input_set_',num2str(date()),'_',num2str(letter),'.mat'];
    if exist(filetosearch)
        letterplus=char(letter+1);
    end
end
Inputs_dat.Graph=Graph;
Inputs_dat.time=time;
save([path2save 'Dynamic_Thorax_Model\Cardio-Pulmonary\Input_set_',num2str(date()),'_',num2str(letterplus),'.mat'],'Inputs_dat');
end

function [Measurements]=configure_measurements(N,letterplus,currentskip,voltageskip,timesteps,d,path2save,path)
%% press everything to get the measurements!
% get the measurements from the files written during the timer loop
Measurements=get_measurements(N,voltageskip,d,timesteps,path);
l1=length(Measurements);
%disp(l1)
nk=mod(l1,N*N);
%hold integer cycles
Measurements=Measurements(1:end-nk);
%and save them
save([path2save 'Dynamic_Thorax_Model\Cardio-Pulmonary\Measurement_set_',num2str(date()),'_',num2str(letterplus),'.mat'],'Measurements');
end

function M=get_measurements(N,voltageskip,d,timestep,path)
M=[];
voltageskip=voltageskip+1;
for j=1:timestep-1
    M1=[];
    modulo=mod(j,N);
    if modulo==0
        modulo=N;
    end
    for i=1:N
        folder=[path '\Dynamic_Thorax_Model\LUA_files\Temporary\'];
        prefix_data=[date(),'_test',num2str(d),'_Voltages',num2str(j),'_',num2str(i)];
        dataformat='.txt';
        filetoread=strcat(folder,prefix_data,dataformat);
        while exist(filetoread)==0
            fprintf('ERROR! Measurements File missing!:')
            filetoread
            fprintf('\nPress anything to continue!')
            pause()
        end
        T1=dlmread(filetoread,'\t',4,0);
        delete(filetoread);
        N1=T1(:,2);
        meas=mean(N1(2:end-1));
        meas1=N1(length(N1)/2);
        meas=max(meas,meas1);
        %bad electrode!
        %if badelectrodes(i)==1||badelectrodes(modulo)==1
        %   meas=meas*2*rand(1,1);
        %end
        M1=[M1 meas];
    end
    if voltageskip==1
        M1=[-diff(M1) M1(N)-M1(1)];
    else
        M2=zeros(1,N);
        for i=1:N
            if i+voltageskip<=N
                M2(i)=M1(i)-M1(i+voltageskip);
            else
                M2(i)=M1(i)-M1(N-i+voltageskip);
            end
        end
        M1=M2;
    end
    M=[M M1];
end
oo=[];
end

function RTG=display_realtime_graph(N,Graph,Ventilation,Circulation,Tissues,electrodes,time,timesteps,Collapsion_areas,times)
subplot(1,2,1)
LL=Tissues.Left_lung.position;
CLL=Graph.left_lung.conductivity(end);
CALL=Collapsion_areas.cross_points.left_lung;
RL=Tissues.Right_lung.position;
if isempty(Graph.right_lung.conductivity)
    CRL=CLL;
else
    CRL=Graph.right_lung.conductivity(end);
end
CARL=Collapsion_areas.cross_points.right_lung;
H=Tissues.Heart.interior.position;
Hmioc=Tissues.Heart.exterior.position;
CH=Graph.heart.conductivity(end);
B=Tissues.Bones;

h1=plot(electrodes.position.x,electrodes.position.y,'bo');
hold on
h2a=plot([Hmioc(:,1); Hmioc(1,1)],[Hmioc(:,2) ;Hmioc(1,2)],'m');
h22a=fill([Hmioc(:,1); Hmioc(1,1)],[Hmioc(:,2) ;Hmioc(1,2)],'m');
h2=plot([H(:,1); H(1,1)],[H(:,2) ;H(1,2)],'r');
h22=fill([H(:,1); H(1,1)],[H(:,2) ;H(1,2)],[CH 0 0]);
h3=plot([LL(CALL(1):CALL(2),1); LL(CALL(1),1)],[LL(CALL(1):CALL(2),2) ;LL(CALL(1),2)],'b');
h33=fill([LL(CALL(1):CALL(2),1); LL(CALL(1),1)],[LL(CALL(1):CALL(2),2) ;LL(CALL(1),2)],[0 0 CLL]);
if ~isempty(CALL)
    h3a=plot([LL(CALL(1):-1:1,1); LL(16:-1:CALL(2),1); LL(CALL(1),1)],[LL(CALL(1):-1:1,2); LL(16:-1:CALL(2),2); LL(CALL(1),2)],'b');
    h33a=fill([LL(CALL(1):-1:1,1); LL(16:-1:CALL(2),1); LL(CALL(1),1)],[LL(CALL(1):-1:1,2); LL(16:-1:CALL(2),2); LL(CALL(1),2)],[0.7 0.1 0.4]);
end
h4=plot([RL(CARL(1):CARL(2),1); RL(CARL(1),1)],[RL(CARL(1):CARL(2),2) ;RL(CARL(1),2)],'b');
h44=fill([RL(CARL(1):CARL(2),1); RL(CARL(1),1)],[RL(CARL(1):CARL(2),2) ;RL(CARL(1),2)],[0 0 CRL]);
if ~isempty(CARL)
    h4a=plot([RL(CARL(1):-1:1,1); RL(15:-1:CARL(2),1); RL(CARL(1),1)],[RL(CARL(1):-1:1,2); RL(15:-1:CARL(2),2); RL(CARL(1),2)],'b');
    h44a=fill([RL(CARL(1):-1:1,1); RL(15:-1:CARL(2),1); RL(CARL(1),1)],[RL(CARL(1):-1:1,2); RL(15:-1:CARL(2),2); RL(CARL(1),2)],[0.7 0.1 0.4]);
end
h5=plot([B.Bone1.position(:,1); B.Bone1.position(1,1)],[B.Bone1.position(:,2) ;B.Bone1.position(1,2)],'k');
h6=plot([B.Bone2.position(:,1); B.Bone2.position(1,1)],[B.Bone2.position(:,2) ;B.Bone2.position(1,2)],'k');
h7=plot([B.Bone3.position(:,1); B.Bone3.position(1,1)],[B.Bone3.position(:,2) ;B.Bone3.position(1,2)],'k');
h8=plot([B.Bone4.position(:,1); B.Bone4.position(1,1)],[B.Bone4.position(:,2) ;B.Bone4.position(1,2)],'k');
h9=plot([B.Bone5.position(:,1); B.Bone5.position(1,1)],[B.Bone5.position(:,2) ;B.Bone5.position(1,2)],'k');
h10=plot([B.Bone6.position(:,1); B.Bone6.position(1,1)],[B.Bone6.position(:,2) ;B.Bone6.position(1,2)],'k');
h11=plot([B.Bone7.position(:,1); B.Bone7.position(1,1)],[B.Bone7.position(:,2) ;B.Bone7.position(1,2)],'k');
h12=plot([B.Bone8.position(:,1); B.Bone8.position(1,1)],[B.Bone8.position(:,2) ;B.Bone8.position(1,2)],'k');
hold off
axis off
RTG='done';
subplot(1,2,2)
t=linspace(0, time,timesteps);
plot(t,Graph.heart.conductivity,'r','LineWidth',2)
hold on
Graph.right_lung.conductivity=Graph.left_lung.conductivity;
plot(t,Graph.right_lung.conductivity,'b','LineWidth',2)
hold on
plot(t,Graph.left_lung.conductivity,'g--','LineWidth',2)
hold on
plot(t,Graph.muscles.conductivity(1:end),'y','LineWidth',2)
hold on
TGr=title('Tissues conductivity (S/m)');
TGr.FontSize=16;
legend('\sigma_{heart}','\sigma_{right lung}','\sigma_{left lung}','\sigma_{muscles}')
ylabel('\sigma (Sm^-^1)','Fontsize',14)
xlabel('Time (seconds)','Fontsize',14)
hold off
end

function intersect=check_intersections(Left_lung,Right_lung,Heart,Bones)
%first check
LLP=Left_lung.shape.full;
RLP=Right_lung.shape.full;
HP=Heart.exterior.shape.full;
Hint=Heart.interior.shape.full;
B1=Bones.Bone1.shape.full;
B2=Bones.Bone2.shape.full;
%B3 is spondilus
B3=Bones.Bone3.shape.full;
B4=Bones.Bone4.shape.full;
B5=Bones.Bone5.shape.full;
B6=Bones.Bone6.shape.full;
B7=Bones.Bone7.shape.full;
B8=Bones.Bone8.shape.full;
b = isintersect(LLP', B3');
if b==0
    b=isintersect(RLP',B3');
    if b==0
        b=isintersect(RLP',HP');
        if b==0
            b=isintersect(LLP',HP');
            if b==0
                b=isintersect(RLP',B1');
                if b==0
                    b=isintersect(RLP',B2');
                    if b==0
                        b=isintersect(LLP',B4');
                        if b==0
                            b=isintersect(LLP',B5');
                            if b==0
                                b=isintersect(LLP',B6');
                                if b==0
                                    b=isintersect(LLP',B7');
                                    if b==0
                                        b=isintersect(HP',B7');
                                        if b==0
                                            b=isintersect(RLP',B7');
                                            if b==0
                                                b=isintersect(RLP',B8');
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
intersect=b;
end

function Patient_labeled_data=complete_patient_labeled_data(Graph,Collapsion_areas,History_of_tissues)

Patient_labeled_data=1;
end