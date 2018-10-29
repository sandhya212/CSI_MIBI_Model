# 17th Aug 2018
# Postprocessing MIBI model - spacial clustering based on clusters, segmentation mask update
# 24th Aug 2018
# colouring segmentation with Leeat's clusters

import matplotlib.pyplot as plt

from matplotlib import colors

import os

import numpy as np
np.set_printoptions(threshold=np.nan)





###colourmap
# choose more colours: https://matplotlib.org/examples/color/named_colors.html
colours=["black", "chocolate", "red", "orange", "limegreen", "blue", "purple", "seagreen","gold","lightpink","thistle","darkviolet","saddlebrown","slategrey",
            "palevioletred","mediumvioletred","yellowgreen","darkolivegreen","lemonchiffon","mistyrose","lightsalmon","lightcyan","lightblue","khaki","wheat","greenyellow","darkmagenta","maroon","dimgray","mediumseagreen"]
#cmap = colors.ListedColormap(colours)

###

##print out Leeat's clusters
path_temp = os.path.join(path_R + '/output_6_markers_Leeat_original_py/plots/extras/cluster_probabilities.csv')
cluster_original_t  = np.array(pd.read_csv(path_temp).as_matrix())[:,1]
col_ind = [int(i) for i in np.unique(cluster_original_t)]
col_ind.insert(0,0)
col_ind_all = [colours[i] for i in col_ind]
cmap = colors.ListedColormap(col_ind_all)

Lp = np.zeros([label.shape[0], label.shape[0]])

for i, regions_i in enumerate(regions_revised):
    cell_label = labelVec[i]
    coords_temp = regions_i['coords']
    Lp[coords_temp[:,0],coords_temp[:,1]] = cluster_original_t[i]

plt.figure()
plt.imshow(Lp,cmap)
plt.imshow(masked, cmap=plt.cm.cool, interpolation='none') #overlaying the segmentation mask on the clustered image

#plt.title("Leeat's CM clusters, K = "+str(int(cluster_original_t.max())))
plt.title("Leeat's CM clusters, K = "+ str(int(np.unique(cluster_original_t).shape[0])))


file_name = os.path.join(path_python + '/figures_python/Leeat_cluster.png')
plt.savefig(file_name)

print(Lp.max())

##print out clusters based on pixel reassignments
num_cells = df_labelvec.shape[0]
path_temp = os.path.join(path_R + '/count_data/meta_data.txt')
#file = open(path_temp, "r")

with open(path_temp) as file1:
    for _ in range(4):
        file1.readline()
    line = file1.readline()
iter = int(line)

for j in range(iter): #max_number_MCMC):
    j_c = j+1
    if(j_c ==1):
        path_temp = os.path.join(
            path_R+'/output_6_markers_updated_1_' + str(j_c) + '/plots/extras/cluster_probabilities.csv')
        cluster_original_t = np.array(pd.read_csv(path_temp).as_matrix())[:, 1]
        col_ind = [int(i) for i in np.unique(cluster_original_t)]
        col_ind.insert(0, 0)
        col_ind_all = [colours[i] for i in col_ind]
        cmap = colors.ListedColormap(col_ind_all)

        Lp = np.zeros([label.shape[0], label.shape[0]])
        for j_cc in range(num_cells):
            #area_cell = np.asscalar(area[j_cc])
            coords_temp = pixet_1_list[j_cc]
            #cumsum = cumsum + area_cell
            # cell_label = labelVec[i]
            Lp[coords_temp[:, 0], coords_temp[:, 1]] = cluster_original_t[j_cc]

    elif(j_c !=1):

        path_temp = os.path.join(path_R+'/output_6_markers_updated_1_'+str(j_c)+'/plots/extras/cluster_probabilities.csv')
        cluster_original_t = np.array(pd.read_csv(path_temp).as_matrix())[:, 1]
        col_ind = [int(i) for i in np.unique(cluster_original_t)]
        col_ind.insert(0, 0)
        col_ind_all = [colours[i] for i in col_ind]
        cmap = colors.ListedColormap(col_ind_all)


        Lp = np.zeros([label.shape[0], label.shape[0]])

        #make the regions_revised from R outputs
        path_area = os.path.join(path_R+ '/output_6_markers_updated_1_' + str(j_c) + '/plots/extras/regions_revised_area_updated.csv')
        path_coords = os.path.join(path_R + '/output_6_markers_updated_1_' + str(j_c) + '/plots/extras/regions_revised_pixelcoords_updated.csv')
        area = np.array(pd.read_csv(path_area, header=None).as_matrix(),dtype=int)
        p_coords = np.array(pd.read_csv(path_coords).as_matrix())
        cumsum = 0
        for j_cc in range(num_cells):
            area_cell = np.asscalar(area[j_cc])
            coords_temp = p_coords[cumsum:(cumsum+area_cell),:]
            cumsum = cumsum + area_cell
            #cell_label = labelVec[i]
            Lp[coords_temp[:, 0], coords_temp[:, 1]] = cluster_original_t[j_cc]

    print(Lp.max())
    plt.figure()
    plt.imshow(Lp, cmap)
    plt.imshow(masked, cmap=plt.cm.cool,
               interpolation='none')  # overlaying the segmentation mask on the clustered image

    #plt.title("Iter_"+str(j_c)+" clusters, K="+str(int(cluster_original_t.max())))
    plt.title("Iter_" + str(j_c) + " clusters, K=" + str(int(np.unique(cluster_original_t).shape[0])))

    file_name = os.path.join(path_python + '/figures_python/Iter_'+str(j_c)+'_clusters.png')
    plt.savefig(file_name)

'''
#### trying1
#Lp_mask = np.matrix(np.logical_and(Lp > 0, Lp < 0),dtype=int)
##effort to make the new seg mask
Lp_1 = np.ones([label.shape[0], label.shape[0]])*-1
Lp_1[Lp == 0] = 1
Lp_1[Lp != 0] = 0
mm = mark_boundaries(Lp, Lp_1, color=(1, 0, 1), mode='inner')
plt.figure()
plt.imshow(Lp, cmap, interpolation='none') #jet or gray
plt.imshow(mm,cmap=plt.cm.cool,interpolation='none')


masked_3 = np.ma.masked_where(Lp_1 <= 0, Lp_1)
plt.figure()
plt.imshow(Lp, cmap, interpolation='none') #jet or gray
plt.imshow(masked_3,cmap=plt.cm.cool,interpolation='none')


%%%try2

Lp_1 = np.zeros([label.shape[0], label.shape[0]])
Lp_1[Lp > 0] = 1
masked_new = np.ma.masked_where(Lp_1 == 0, Lp_1)
plt.imshow(masked_new)

label_L = measure.label(Lp, connectivity=neighbourhood_conn)
regions_L = measure.regionprops(label_L)
labelNum_L = label_L.max()

%%%try3

from scipy import ndimage
a = np.zeros((50,50))
a[10:30,10:30] = 1
a[35:45,35:45] = 2
distance = ndimage.distance_transform_edt(a)
distance[distance != 1] = 0
plt.imshow(distance)
plt.show()
np.where(distance == 1)



distance11 = ndimage.distance_transform_edt(Lp_1)
distance11 = ndimage.distance_transform_edt(Lp)
distance11[distance11 != 1] = 0
plt.imshow(distance11)
plt.show()
'''

#pandas.Series(cluster_original_t).value_counts()


