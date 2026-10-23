Mathematical Modelling of Acute Myeloid Leukaemia (AML)

This repository contains MATLAB code for a mathematical model investigating bone marrow population dynamics in AML.

Model:
The model consists of four coupled ODEs representing:
H: Healthy myeloblasts
S: Cytarabine-sensitive AML blasts
R: Cytarabine-resistant AML blasts
B: Other healthy bone marrow cells

Executing the code:
The MATLAB script AML_Model_Final.m can be executed through GitHub using the Run AML MATLAB Model workflow. The system is solved using MATLAB's ode45 solver over a specified simulation period of 17520 hours. Chemotherapy initiation times are specified manually.
The following settings can be modified to investigate different scenarios:
y0: Initial populations of H, S, R and B
chemoON: Set to true or false to enable or disable cytarabine administration
secondChemoON: Set to true or false to enable or disable a second chemotherapy cycle
chemoStart1: Starting time of the first chemotherapy cycle (hours)
chemoStart2: Starting time of the second chemotherapy cycle (hours)
tEnd: Total simulation duration (hours)
threshold: Total blast percentage used to calculate upward threshold crossing times

Note: The threshold parameter is used to calculate the times at which the selected total blast percentage is crossed from below. Chemotherapy initiation times must then be entered manually.

Results:
The code generates 4 figures illustrating:
1. Bone marrow population dynamics.
2. Total normalised bone marrow population and carrying capacity.
3. Relative population composition
4. Cytarabine administration schedules.
These results are presented and discussed in the research report.
