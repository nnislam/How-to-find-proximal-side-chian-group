#!/anaconda3/bin/pythonw

######################################
## Radial Distribution Function ##
## Naeyma Islam ##
## UNO, 31st Mar 2018 ##
######################################
import pandas as pd
import numpy as np
import math as mth
#----- Importing the libraries -----#

input_file_name = "pmaa.lammpstrj"  # .lammpstrj file name along with location path
out_file_name = "num_proxy_deprox.out"   # output file name along with path to store rgyr values

print ("Start processing Trajectory file....\n")

#----- Reading .lammpstrj file (from lammps output) -----#

cols_read = ['TIMESTEP', 'ATOM_ID', 'ATOM_TYPE', 'x', 'y', 'z', 'xbox', 'ybox', 'zbox']
cols_o2 = ['R_DIST','x_o1', 'y_o1', 'z_o1','x_o2','y_o2','z_o2']
input_data = pd.DataFrame(columns=cols_read)    # Dataframe for the Timestep data
df_o2 = pd.DataFrame(columns=cols_read)
#df_na = pd.DataFrame(columns=cols_read)
df_o1_o2 = pd.DataFrame(columns=cols_o2)
df_o1_o2_gt = pd.DataFrame(columns=cols_o2)
#df_na_count = pd.DataFrame(columns=cols_o_na)
input_one_row = []

atom_type_o6 = 8  # Oxygen
#atom_type_o7 = 7  # Oxygen
#atom_type_na = 8  # Na; Water

dr = 0.1

xbox = 0
ybox = 0
zbox = 0
vol_box_sum = 0
time_step = -1
frame_count = 0



cols_output = ['Frame', 'Na_Count']
output_data = pd.DataFrame(columns=cols_output)  # Dataframe to store the calculated Rgyr values along with Time
new_output = []

try:
    input_file = open(input_file_name, "r")  #opens the .data file in the working directory in 'reading mode'
    out_file = open (out_file_name, "w")   # opening the out file to write rgyr values
    
    no_of_atoms = 0

    for line in input_file:          
        #print (line)                          # line points to the first line of the FRAME ('ITEM: TIMESTEP')
        next_line = input_file.readline()      # pointing to the 'ITEM: TIMESTEP' line 
        #print (next_line)
        line_as_list = next_line.split()
        time_step = int(line_as_list[0])       # value of the TIMESTEP/FRAME
        
        print ("->Processing Frame: %d" % time_step)
        frame_count = frame_count + 1          # Increasing frame counter
        
        text_line = input_file.readline()      # pointing to the 'ITEM: NUMBER OF ATOMS' line
        next_line = input_file.readline()      # pointing to the next line of 'ITEM: NUMBER OF ATOMS'
        #print (next_line)
        line_as_list = next_line.split()
        no_of_atoms = int(line_as_list[0])     # number of atoms per frame
        
        text_line = input_file.readline()      # pointing to the 'ITEM: BOX BOUNDS' line
       # print (text_line)
        
        next_line = input_file.readline()      # pointing to the line having X low and high for the box
        #print (next_line)
        line_as_list = next_line.split()
        xbox = float(line_as_list[1]) - float(line_as_list[0])  # X distance for Box
        
        next_line = input_file.readline()      # pointing to the line having Y low and high for the box
        line_as_list = next_line.split()
        ybox = float(line_as_list[1]) - float(line_as_list[0])  # Y distance for Box
        
        next_line = input_file.readline()      # pointing to the line having Z low and high for the box
        line_as_list = next_line.split()
        zbox = float(line_as_list[1]) - float(line_as_list[0])  # Z distance for Box
        
        vol_box_sum = vol_box_sum + (xbox * ybox * zbox)  # adding-up Volume of the Box
        
        text_line = input_file.readline()      # pointing to the 'ITEM: ATOMS' line
            # reading all atoms in the frame
        for i in range (0, no_of_atoms):
            next_line = input_file.readline()  # pointing to the line of each Atoms
            #print (next_line)
            line_as_list = next_line.split()
            if line_as_list and line_as_list[0].isdigit():   # check whether the line is not empty and the first word is a number
                atom_type = int(line_as_list[1])
                
                if (atom_type == atom_type_o6):
                    x = float(line_as_list[2])   
                    y = float(line_as_list[3])  
                    z = float(line_as_list[4])  
                    
                    input_one_row = [time_step, int(line_as_list[0]), int(line_as_list[1]), x, y, z, xbox, ybox, zbox]   # loading one row data
                    new_row_df = pd.DataFrame ([input_one_row], columns=cols_read, index=[int(line_as_list[0])])          # making dataframe for loaded row
                    #input_data = input_data.append(new_row_df)
                    df_o2 = df_o2.append(new_row_df)
                    input_one_row = []
                   
        
        for i in range(int(df_o2.shape[0] - 1)):          # for each Oxygen of Type - 6 & 7
            r_dist = 0
            dx = df_o2.iloc[i+1,3] - df_o2.iloc[i,3]
            dy = df_o2.iloc[i+1,4] - df_o2.iloc[i,4]
            dz = df_o2.iloc[i+1,5] - df_o2.iloc[i,5]
            dx = dx - (xbox * round (dx/xbox))
            dy = dy - (ybox * round (dy/ybox))
            dz = dz - (zbox * round (dz/zbox))
            r_dist = mth.sqrt (mth.pow(dx,2) + mth.pow(dy,2) + mth.pow(dz,2))
            r_dist = round(r_dist, 1)
            if r_dist < 3.7:
                one_row_data = [r_dist, df_o2.iloc[i,3], df_o2.iloc[i,4], df_o2.iloc[i,5],df_o2.iloc[i+1,3],df_o2.iloc[i+1,4],df_o2.iloc[i+1,5]]
                
                one_row_df = pd.DataFrame ([one_row_data], columns=cols_o2)
                df_o1_o2 = df_o1_o2.append(one_row_df)
             
            else:
                one_row_data = [r_dist, df_o2.iloc[i,3], df_o2.iloc[i,4], df_o2.iloc[i,5],df_o2.iloc[i+1,3],df_o2.iloc[i+1,4],df_o2.iloc[i+1,5]]
                
                one_row_df = pd.DataFrame ([one_row_data], columns=cols_o2)
                df_o1_o2_gt = df_o1_o2_gt.append(one_row_df)
                
        #print(df_o_na.shape[0]) 
                           
                 
        out_file.write(str(time_step) + "\t" +  str(df_o1_o2.shape[0]) + "\t" +  str(df_o1_o2_gt.shape[0])+ "\n")
        
        no_of_atoms = 0
        xbox = 0
        ybox = 0
        zbox = 0
        input_data = input_data.iloc[0:0]
        df_o2 = df_o2.iloc[0:0]
        #df_na = df_na.iloc[0:0]
        df_o1_o2 = df_o1_o2.iloc[0:0]
        df_o1_o2_gt = df_o1_o2_gt.iloc[0:0]
            
    

    out_file.close()                    # closing the out file
    input_file.close()                    # closing the opend trj file
    
except ValueError:
    print ("Cannot find the file or wrong file format. Please check and try again.")
    system.exit(os.EX_CONFIG)

#----- Graph plotting -----#

#print ("Ploting Radial Distribution Function ....\n")

#plt.plot(output_data['R'], output_data['RDF_O2'])
##plt.plot(output_data['R'], output_data['RDF_O7'])

#plt.ylabel('G(R)')
#plt.xlabel('R')
#plt.plot
#plt.show()



