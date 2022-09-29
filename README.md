## 1 About this program
usermat. f is a program written in Fortran to simulate the damage of brittle material in ANSYS.
## 2 Usage method
Step 1: Setting the operating environment. Click the [link](URL "https://www.bilibili.com/read/cv10433898?from=search&spm_id_from=333.337.0.0") for reference.
（https://www.bilibili.com/read/cv10433898?from=search&spm_id_from=333.337.0.0）
Step 2: Put the usermat. f into ANSYS, and then relink ANSYS.    
Step 3: Set the physical parameters of the material, including strain damage threshold and strain damage limit, to simulate the damage of the material.
## 3 User defined material
Define 6 parameters of the brittle material, corresponding to EMD, ENU, Etf, Etu, Ecf, and Ecu in the usermat. f. For example：    
Tb,user,1,1,6    
Tbdata,1,300e8,0.2,0.0003,0.0033,0.0008,0.0023    
Tb,state,1,,1    
