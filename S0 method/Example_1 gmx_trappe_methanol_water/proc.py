import shutil, subprocess
import os
 
x_met = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
file_path = os.getcwd()


for i in range(len(x_met)):
    met_mol = x_met[i]
    os.chdir(file_path)
    
    new_file = f"{met_mol}_met_water" 
    shutil.copy("trj_cal.py", new_file) 
    shutil.copy("convect.sh", new_file)

    os.chdir(new_file)
    command_1 =f"python trj_cal.py"
    command_2 =f"sh convect.sh"
    subprocess.run(command_1, shell=True)
    subprocess.run(command_2, shell=True)
    
    command_4 = "gmx energy -f eql.edr -o rho.xvg <<< \"20 0\" "
    subprocess.run(command_4, shell=True)
    print(f"{met_mol}, done! :D")
    
