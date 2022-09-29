## 1 About this program
usermat. f is a program written in Fortran to simulate the damage of brittle material in ANSYS.
## 2 Usage method
Step 1: Setting the operating environment. Click the [link](URL "https://www.bilibili.com/read/cv10433898?from=search&spm_id_from=333.337.0.0") for reference.   
Step 2: Put the usermat. f into ANSYS, and then relink ANSYS.    
Step 3: Set the physical parameters of the material, including strain damage threshold and strain damage limit, to simulate the damage of the material.
## 3 Ansys instruction
### 3.1 User defined material

Tb,user,1,1,6    
Tbdata,1,300e8,0.2,0.0003,0.0033,0.0008,0.002    
Tb,state,1,,1    
