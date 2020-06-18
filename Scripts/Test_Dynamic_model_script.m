path='C:\Users\Chris\Desktop\Bio-lungs\';
path2eidors='C:\Users\Chris\Desktop\PHD_files\eidors-v3.9-ng\eidors\';
paths={path, path2eidors};

%%%%Measurement Parameters

%No of electrodes
N=16;
%current skip-m
skipcurr=2;
%voltage skip-n
skipvolt=0;
%ac current frequency (Hz), f=200kHz
f=200000;
%total simulated image frames 
total_frames=10;
%frames per second 
fps=10;
%ac current amplitude
Current_Amplitude=0.001; %Amperes
%pack measurement parameters
Measurement_params=[N skipcurr skipvolt f total_frames fps Current_Amplitude];

%%%%Model Parameters

breath_time=3; %seconds
Zelectrodes(:,1)=0.02*(randn(N,1)*0.1+1);
Zelectrodes(:,2)=300*(randn(N,1)*0.3+1);
delec=0.5; %electrode_width (cm)
random_start=0;
random_state=0;
%deflated lung conductivity at 200kHz ~0.25S/m []
%inflated lung conductivity at 200kHz ~0.1S/m []
slairmax=0.25; %(S/m)
end_Admittances(1)=slairmax;
slairmin=0.1; %(S/m)
end_Admittances(2)=slairmin;
%blood cycle-related conductivity variance in lungs ~0.008S/m (10-20 times
%less than air-related [])
dslblood=0.008; %(S/m)
end_Admittances(3)=dslblood;

%deflated lung permittivity at 200kHz ~4000F/m []
%inflated lung permittivity at 200kHz ~2000F/m []
elairmax=4000; %(S/m)
end_Admittances(4)=elairmax;
elairmin=2000; %(S/m)
end_Admittances(5)=elairmin;
%blood cycle-related permittivity variance in lungs ~100 F/m (10-20 times
%less than air-related [])
delblood=100;
end_Admittances(6)=delblood;

%base heart (champers, aorta, no-myocardium) conductivity ~0.55 S/m at
%200kHz
Sheart=0.55; %(S/m)
end_Admittances(7)=Sheart;
%base heart (champers, aorta, no-myocardium) conductivity variation ~0.025 S/m at
%200kHz (~15-20%)[]
dSheart=0.025;
end_Admittances(8)=dSheart;

%base heart (champers, aorta, no-myocardium) conductivity ~6000 F/m at
%200kHz
Eheart=6000;
end_Admittances(9)=Eheart;
%base heart (champers, aorta, no-myocardium) permittivity variation ~300 S/m at
%200kHz (~15-20%)[]
dEheart=300;
end_Admittances(10)=dEheart;

%HR
bpm=75;
%Collapsions states for Left and Right lung. From 1, which denotes an 20%
%lung collapsion to 4 which denotes a healthy lung. 
Collapse_LL=1;
Collapse_RL=3;

%pack model parameters
Model_params=[end_Admittances breath_time delec bpm Collapse_LL Collapse_RL];

%%%% call the model constructor
[M,Gr,Patient_labeled_data]=Dynamic_Thorax_Imaging(Measurement_params,...
    Model_params,Zelectrodes,random_start,random_state,paths);


%%%%%References

