#17th Aug 2018

import rpy2.robjects as robjects

r_source = robjects.r['source']
r_source(os.path.join(path_R +"/start_file_R.R"))


# dummy R code to test bits
#ot=r_source(os.path.join(path_t + '/superpixel_MIBI_python/dummy.R'))
