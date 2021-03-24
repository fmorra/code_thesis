import os
import pdb


filelist = []
path = 'C:\\Users\\fabri\\Desktop\\tu_delft\\Thesis\\ABAQUS\\test_run_python'
for root, dirs, files in os.walk(path, topdown=False):
    for name in files:
        if name == "Complete_file_EARTH.csv" or name == "Geographical_complete_file_EARTH.csv":
            filelist.append(os.path.join(root, name))
            print(os.path.join(root, name))
print(filelist)
