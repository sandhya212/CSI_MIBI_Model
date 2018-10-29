###3rd August 2018
## Code author SP

###########
### Python libraries
############

import os


#import matplotlib
#matplotlib.use('Agg')




############
#### define path variables for data, Python code, R code etc
############


path_data = '/Users/prabhaks/Documents/Pathology_images/SP_TEST_Image/Halo_archive_2017-04-21_15-04/MIBI_set2/TNBCShareDataWithPeerGroup'
path_t = os.getcwd() #'/Users/prabhaks/Desktop/Python_R_codebase'
path_python = path_t #os.path.join(path_t+"/MIBI_model"); just in case you have code and data in same file level, then use uncomment
path_R = os.path.join(path_t+ "/superpixel_MIBI_python")

execfile(os.path.join(path_python+"/init_file.py"))


