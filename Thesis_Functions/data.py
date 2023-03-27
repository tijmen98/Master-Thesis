"""All Data specific functions"""

def Variable_Import(directory,variable):
    
    import pandas as pd
    import xarray as xr
    import numpy as np
    
    directory_information = pd.read_csv(directory+'Directory_information.csv',sep=';',names=['file','long_name','units','variable','plot'])
    directory_information.drop(index=np.nan,inplace=True)
    
    for num in range(len(directory_information)):
    
        if directory_information['variable'].values[num] == variable:
            dataset = xr.open_dataset(directory+directory_information['file'].values[num])
            dataset.rotated_pole
            
            return(dataset[directory_information['variable'].values[num]])
        
        
        
def Create_Directory_Information (directory,variable_split):
    
    print('The following excemptions occured when creating the Directory information file:')
    
    import xarray as xr
    import os
    import pandas as pd

    files = os.listdir(directory)
    
    information = [['file'],['long_name'],['units'],['variable'],['plot']]
    
    for num in range(len(files)):

        if files[num]!='Directory_information.csv':
            variable = files[num].split(variable_split)[0]
            try:
                ds = xr.open_dataset(directory+files[num])
            except:
                print('    File '+files[num]+' could not be opened, index is: '+str(num))
                ds=0
            else:
                ds = xr.open_dataset(directory+files[num])
            
        
                information[0].append(files[num])
            
                
                #print(ds[variable].long_name)
                
                try:
                    ds[variable].long_name
                except:
                    information[1].append(ds.title) 
                else:
                    information[1].append(ds[variable].long_name)
                
                
                try:
                    ds[variable].units
                except:
                    print('    No unit found, index is: ' + str(num))
                else:
                    information[2].append(ds[variable].units)
                
                information[3].append(variable)
            
            df = pd.DataFrame(information).transpose()
            df.rename({0: 'file', 1: 'long_name',2:'units',3:'variable',4:'plot'},axis='columns',inplace=True)
            df.drop(index=0,inplace=True)
        
            df.to_csv(directory+'Directory_information.csv',sep=';')
    print('Directory information saved at'+ directory)
    return()

"""returns sea mask for plotting"""
    
def Get_Seamask(directory,time=0):
    tile1 = Variable_Import(directory,'tilefrac1').isel(time=time)+1
    tile2 = Variable_Import(directory,'tilefrac3').isel(time=time)+1
    return(tile2.where(tile2==1)*tile1.where(tile1==1))