# Dynamic-Thoracic-Model-for-EIT

Dynamic 2D Thoracic Electrical Impedance Tomography Simulation Model


//--This project has been published in 2020 IEEE 20th International Conference on Bioinformatics and Bioengineering (BIBE) Proceedings 
with the title "A highly tunable dynamic thoracic model for Electrical Impedance Tomography"
DOI: 10.1109/BIBE50027.2020.00163 and is available at https://ieeexplore.ieee.org/document/9288054  --//


Requirements:

MATLAB 2018b or later

EIDORS open source library tool (v.3.9 or later)  http://eidors3d.sourceforge.net/ [1]

FEMM open software for 2D F.E.M. meshing  http://www.femm.info/wiki/HomePage


Simulation times depend on the the total number of steps (electrode number X total frames) and the CPU's speed. 

It is strongly recommended to minimize the FEMM window during the simulation, in order to significantly increase simulation speeds. 


---------------INSTRUCTIONS------------------------

Step 1) Initialize the input parameters in Test_Dynamic_model_script.m  REMEMBER TO SET THE paths!!! (to Bio-lungs and EIDORS folders, see below at F)


Those include:

A. The measurement parameters:
	-No of electrodes (N)
	-current skip-m (skipcurr)
	-voltage skip-n (skipvolt)
	-ac current frequency (Hz)
	-total simulated image frames (total_frames)
	-frames per second (fps)
	-ac current amplitude (A)
	packed in Measurement_params input vector

B. The Model Parameters: parameters that define the measuring properties:
1. end_Admittances: includes the inspiration/ expiration-end conductivities and permittivities of each tissue

    -deflated lung conductivity (S/m)  end_Admittances(1)
    
    -inflated lung conductivity (S/m)  end_Admittances(2)
    
    -blood cycle-related conductivity variance in lungs (S/m) end_Admittances(3)
    
    -deflated lung permittivity (F/m)  end_Admittances(4)
    
    -inflated lung permittivity (F/m)  end_Admittances(5)
    
    -blood cycle-related permittivity variance in lungs (F/m)  end_Admittances(6)
    
    -base heart (champers, aorta, no-myocardium) conductivity (S/m) end_Admittances(7)
    
    -base heart (champers, aorta, no-myocardium) conductivity variation (S/m) end_Admittances(8)
    
    -base heart (champers, aorta, no-myocardium) conductivity (S/m) end_Admittances(9)
    
    -base heart (champers, aorta, no-myocardium) permittivity variation (F/m) end_Admittances(10)
    
2. breath_time (seconds), Initial Breathing Time
3. delec (cm), the electrodes width
4. bpm, Initial Heart Rate
5. Collapse_LL, Left Lung Collapsion state (1-4)
6. Collapse_RL, Right Lung Collapsion state (1-4)

C. Zelectrodes: A 2XN matrix. The first raw contains each electrode's contact relative conductivity. The second raw contains each electrode's contact relative permittivity. 

D. random_start: 0 if simulation starts from pulmonary and cardiac states 1 (full exhalation and start of heartbeat). 1 if simulation starts from random states at both cycles. 

E. randomcond: if 1: breathing times and HR are randomly changed every cycle and 3 cycles respectively. Else if 0: initial breathing time is 3 seconds, 
silent space is 0.3seconds and next breath lasts +20%. Initial heart beat rate is 75bpm and increases by 3% every 3 beats. 

F. paths: A cell array. paths{1} denotes the directory to the Bio-lungs folder ('C:\.....\Bio-lungs\'), paths{2} denotes the directory to the Eidors library ('C:\.....\eidors-v3.9-ng\eidors\')


Step 2) Run Test_Dynamic_model_script.m .  This calls the Dynamic_Thorax_Imaging function which takes the defined parameters and generates a .lua script. 

Step 3) Wait for the .lua script to LOAD

Step 4) Open a FEMM (open current-flow problem). The go to file->Open Lua script and open the created .lua file  

Step 5) Wait until FEMM performs the simulations for the dynamic 2D structures. For faster process, MINIMIZE the FEMM window. 
This Process may take from some minutes to a couple of hours, depended on the total frames and electrodes simulated, and the CPU

Step 6) When FEMM finishes, press enter in MATLAB console and wait until the final noise-free measurements load. 

Step 7) RUN reconstruction_and_functional_analysis.m WITH THE SAME MEASUREMENT PARAMETERS (N, SKIPCURR, SKIPVOLT).

Step 8) Wait till reconstructed images appear


[1] Adler, Andy, and William RB Lionheart. "Uses and abuses of EIDORS: an extensible software base for EIT." Physiological measurement 27.5 (2006): S25.
