## MIBI segmentation correction code
## Code author SP
##
## Date - 30th July 2018
##
##
## Revision history
## 3rd aug - read the files in, built Leeat's count matrix, read the nuc markers in
## 7th Aug - built pixet1, CM based on pixet1
## Code with square tessellation

############
#### Main code
############
import numpy as np
np.set_printoptions(threshold=np.nan)

import time

from PIL import Image, ImageDraw
from skimage import measure #package is scikit-image

import pandas as pd
import os

import matplotlib.pyplot as plt


from skimage.segmentation import mark_boundaries
from skimage.color import rgb2gray
from skimage.filters import sobel


from scipy import ndimage as ndi
from skimage.feature import peak_local_max

from joblib import Parallel, delayed
import multiprocessing

#from skimage.segmentation import watershed


#import matplotlib
#matplotlib.use('Agg')

def build_pixet2(i):
    neigh_WSS_temp = []
    WSS_coords = regions_WSS_cleaned[i]
    l = setdiff2d(WSS_coords, sure_coords)
    if (len(l) != 0):
        print(i)
        pixet_2_list.append(l) #(l.tolist()) each (x,y) coord is one cell in the .csv, else l is the np.ndarray and is written as so
        data_WSS[i, :] = image_tensor[l[:, 0], l[:, 1], :].sum(axis=0)
        WSS_Sizes[i] = len(l)
        WSS_dataCentroids[i, :] = np.mean(l,axis=0)
        WSS_dataScaleSize[i, :] = data_WSS[i, :] / WSS_Sizes[i]
        WSS_dataScaleCentroids[i, :] = WSS_dataCentroids[i, :] / WSS_Sizes[i]
    else:
        pixet_2_list.append(fill_in)


    len_coords = len(WSS_coords)
    for j in range(len_coords):
        # print(j)
        neigh_WSS_temp.extend(cell_neighbors(label_sp2, WSS_coords[j, 0], WSS_coords[j, 1], dp))
        nt = np.unique(neigh_WSS_temp) - 1  # no idea why the cell labels are incremented by 1, check the cell_neighbourhood logic, for now I am just decrementing this
        nt = [y for y in nt if y != -1]  # remove -1. This means cell 0. a) real cell 0 will never be a neighbour b) why it is even there for every cell??
        nt = [y for y in nt if y != i]  # artefact of cell neighbourhppd function. Every cell is a neighbour of itself. So I remove this. But again the -1 from before needs to be tackled.
    neighbourhood_WSS_list.append(nt)

def assign_cell_id(i, r_sp2):
    print("i "+str(i))
    for j, r in enumerate(regions_revised):
        l = intersect2d(regions_WSS_cleaned[i], r['coords'])
        if (len(l)!=0):
            #print('j'+str(j))
            WSS_cell_label_assgn[0,i]=labelVec[j]
            #print(labelVec[j])
            break

t1 = time.time()

#read in the MIBI MassSpec data
massDS = np.array(pd.read_csv(os.path.join(path_data+'/TNBCpanelInfo'+"."+ 'csv')).as_matrix())
valid_channelNum = len(clusterChannels)
se = np.ones([window_neigh, window_neigh])

for p in range(0,points):

    print('processing patient # ', p+1)
    pointNumber = p+1

    #read the patient's .tiff file
    newLmod = Image.open(os.path.join(path_data + '/Point' + str(pointNumber)+ '/SegmentationInterior.tif'))
    newLmod_arr = np.array(newLmod)
    plt.figure()
    plt.imshow(newLmod_arr  * 200, 'gray', interpolation='none')  # gray
    plt.title('Segmentation mask full')
    file_name = os.path.join(path_python+'/figures_python/Seg_mask_full.png')
    plt.imsave(file_name, newLmod_arr * 200, cmap='gray')

    patch_image = int(newLmod_arr.shape[0] / divisor_patch)
    newLmod_arr = newLmod_arr[0:patch_image, 0: patch_image]

    plt.figure()
    plt.imshow(newLmod_arr * 200, 'gray', interpolation='none')  # gray
    plt.title('Segmentation mask patch')

    file_name = os.path.join(path_python + '/figures_python/Seg_mask_patch.png')
    plt.imsave(file_name, newLmod_arr * 200, cmap='gray')

    label = measure.label(newLmod_arr, connectivity=neighbourhood_conn)
    regions = measure.regionprops(label)
    labelNum = label.max() # number of segmented regions


    data = np.zeros([labelNum, valid_channelNum], dtype=np.float)
    dataScaleSize = np.zeros([labelNum, valid_channelNum], dtype=np.float)
    cellSizes = np.zeros([labelNum, 1], dtype=np.float)
    dataCentroids = np.zeros([labelNum, 2], dtype=np.float)
    dataScaleCentroids = np.zeros([labelNum, 2], dtype=np.float)

    #neighbourhood_list = dict()
    neighbourhood_list = []

    #build the image tensor from the marker channels

    image_tensor = np.zeros((newLmod_arr.shape[0], newLmod_arr.shape[0], valid_channelNum ))
    for t in range(0, valid_channelNum):
        image_temp = np.array(Image.open(os.path.join(path_data + '/Point' + str(pointNumber) + '/'+ clusterChannels[t] + '.tif')))
        image_temp = image_temp[0:patch_image, 0: patch_image]
        #image_tensor.append(image_temp)
        image_tensor[:,:,t] = image_temp

    #build the count matrix along with cell Area, centroids, scaled and unscaled as well as neighbourhoods for each cell

    start_neigh = time.time()

    neighbourhood_list = []
    print("Extracting neighbours")
    for i, regions_i in enumerate(regions):  # check for pixel entries too
        neigh_temp = []

        coords_temp = regions_i['coords']
        cellSizes[i] = regions_i['Area']
        dataCentroids[i, :] = regions_i['centroid']

        data[i, :] = image_tensor[coords_temp[:, 0], coords_temp[:, 1], :].sum(
            axis=0)  # sum the columns, per ith row, top down, length of vector
        dataScaleSize[i, :] = data[i, :] / cellSizes[i]
        dataScaleCentroids[i, :] = dataCentroids[i, :] / cellSizes[i]
        # find neighbours per region (aka cell)


        plt.text(dataCentroids[i, 0], dataCentroids[i, 1], str(i), fontsize=7, color='red', verticalalignment='center',
                 horizontalalignment='center')

        '''
        ####logic 1 - make a patch (se), superimpose on cell and get all the neighbours by logical &
        temp_arr = np.zeros([patch_image, patch_image])
        temp_arr[coords_temp] = 1
        im_dilate = ndi.binary_dilation(temp_arr, structure=se)  # .astype(temp_arr.dtype)
        neigh = np.logical_and(im_dilate, 1 - temp_arr)
        neigh_labels = np.unique(label[neigh])
        neigh_labels = neigh_labels[neigh_labels != 0]
        print('neighbouring cells for cell ' + str(i) + ':' + str(neigh_labels))
        neighbourhood_list.append(neigh_labels)
        # neighbourhood_list[str(i)] = neigh_labels
        # del neigh_labels
        '''
        ####logic 2 - make a patch (dp), slide it over cell and get all the neighbours
        len_coords = len(coords_temp)
        for j in range(len_coords):
            # print(j)
            neigh_temp.extend(cell_neighbors(label, coords_temp[j, 0], coords_temp[j, 1], dp))
            nt = np.unique(neigh_temp) - 1  # no idea why the cell labels are incremented by 1, check the cell_neighbourhood logic, for now I am just decrementing this
            nt = [y for y in nt if y != -1]  # remove -1. This means cell 0. a) real cell 0 will never be a neighbour b) why it is even there for every cell??
            nt = [y for y in nt if y != i]  # artefact of cell neighbourhppd function. Every cell is a neighbour of itself. So I remove this. But again the -1 from before needs to be tackled.
        neighbourhood_list.append(nt)
    end_neigh = time.time()
    print('Time to get neighbours for cells: ' + str(end_neigh - start_neigh))


    #clean up the count matrix
    #get the final information only for cell labels that are
    #1. positive nuclear identity (cells)
    #2. Cells that have enough information in the clustering channels to be clustered.


    labelIdentityNew2 = np.ones([1, labelNum])
    sumDataScaleSizeInClusterChannels = dataScaleSize.sum(axis=1)
    ind_1 = np.flatnonzero(sumDataScaleSizeInClusterChannels < 0.1)
    labelIdentityNew2[:,ind_1] = 2
    labelVec = np.flatnonzero(labelIdentityNew2 == 1) # labels of valid cells
    dataCells = data[labelVec,:] #Leeat's CM
    dataScaleSizeCells = dataScaleSize[labelVec,:]
    cellCentroids = dataCentroids[labelVec,:]
    cellScaleCentroids = dataScaleCentroids[labelVec,:]
    cellSizes_revised = cellSizes[labelVec] #area

    valid_cells = len(labelVec)


    ### writing to a .csv file

    df_data = pd.DataFrame(dataCells)
    df_data.to_csv(os.path.join(path_python + '/superpixel_MIBI_python/count_data/count_matrix_LM.csv'),header=None, index=None)
    df_datascaled = pd.DataFrame(dataScaleSizeCells)
    df_datascaled.to_csv(os.path.join(path_python + '/superpixel_MIBI_python/count_data/count_matrix_scaled_LM.csv'), header=None,
                   index=None)
    #labelVec - valid cells
    df_labelvec = pd.DataFrame(labelVec)
    df_labelvec.to_csv(os.path.join(path_python + '/superpixel_MIBI_python/count_data/valid_cells.csv'), header=None,
                   index=None)




    #create a new regions structure with just the valid cells
    regions_revised = [y for x,y in enumerate(regions) if x in labelVec] #can also use contains(labelvec,x)
    neighbourhood_list_revised = [y for x,y in enumerate(neighbourhood_list) if x in labelVec]
    #remove those invalid cells from neighbourhoods as well
    for i,neigh_i in enumerate(neighbourhood_list_revised):
        neighbourhood_list_revised[i] = np.intersect1d(neigh_i, labelVec)

    #area, centroid and pixel coordinates for regions_revised
    regions_revised_area = []
    regions_revised_centroid = []
    regions_revised_pixelcoords = []
    for i, r_i in enumerate(regions_revised):
        regions_revised_area.append(r_i['Area'])
        regions_revised_centroid.append(r_i['centroid'])
        regions_revised_pixelcoords.append(r_i['coords'])


    #regions area
    df_rra = pd.DataFrame(regions_revised_area)
    df_rra.to_csv(os.path.join(path_python + '/superpixel_MIBI_python/count_data/regions_revised_area.csv'), header=None,
                   index=None)

    #regions centroid
    df_rrcent = pd.DataFrame(np.array(regions_revised_centroid))
    df_rrcent.to_csv(os.path.join(path_python + '/superpixel_MIBI_python/count_data/regions_revised_centroid.csv'), header=None,
                   index=None)

    #regions pixel coords
    df_rrpcoord = pd.DataFrame(np.array(regions_revised_pixelcoords))
    df_rrpcoord.to_csv(os.path.join(path_python + '/superpixel_MIBI_python/count_data/regions_revised_pixelcoords.csv'), header=None,
                   index=None)



    ############################
    seg_mask_temp = np.array(Image.open(os.path.join(path_data + '/Point' + str(pointNumber)+ '/Segmentation.tif')))
    seg_mask = seg_mask_temp[0:patch_image, 0: patch_image]
    masked = np.ma.masked_where(seg_mask == 0, seg_mask)


    ############
    #### Building pixet1
    ############
    #read in the 3 nuclear markers - dsDNA, H3K9ac and H3K27me3 plus the segmentation mask


    len_cer_markers = len(certain_markers)
    cer_markers_tensor = np.zeros((newLmod_arr.shape[0], newLmod_arr.shape[1], len_cer_markers))
    certain_markers_all = np.zeros([newLmod_arr.shape[0], newLmod_arr.shape[1]],dtype=np.float)


    sure_coords = np.array([[0, 0]]) #pixet1

    for i in range(len_cer_markers):
        image_temp = np.array(Image.open(os.path.join(path_data + '/Point' + str(pointNumber) + '/' + certain_markers[i] + '.tif')))
        cer_markers_tensor[:, :, i] = image_temp[0:patch_image, 0: patch_image]
        sure_coords_temp = np.array(np.where(cer_markers_tensor[:, :, i] >= np.median(cer_markers_tensor[:, :, i].max(axis=0)) - threshold_pixels))
        sure_coords = np.concatenate((sure_coords, sure_coords_temp.T),axis=0)
        certain_markers_all = certain_markers_all + cer_markers_tensor[:, :, i]

    sure_coords = np.delete(sure_coords, (0), axis=0) #remove the first row of zeros

    # Efficiently Remove Duplicate Rows from a 2D Numpy Array
    x_dummy = np.random.rand(sure_coords.shape[1])
    y_dummy = sure_coords.dot(x_dummy)
    unique, index = np.unique(y_dummy, return_index=True)
    sure_coords = sure_coords[index]



    plt.figure(frameon=False)
    plt.imshow(certain_markers_all*20,'gray', interpolation='none') #jet or gray
    plt.imshow(masked,cmap=plt.cm.cool,interpolation='none')
    plt.title('Segmentation mask on all nuclear markers')
    file_name = os.path.join(path_python + '/figures_python/Seg_mask_nuc_marker.png')
    plt.savefig(file_name)

    pixet1_image = np.zeros([newLmod_arr.shape[0], newLmod_arr.shape[1]], dtype=np.float)
    pixet1_image[sure_coords[:,0],sure_coords[:,1]] = certain_markers_all[sure_coords[:,0],sure_coords[:,1]]
    plt.figure()
    plt.imshow(pixet1_image*200,'gray', interpolation='none') #gray
    plt.imshow(masked, cmap=plt.cm.cool, interpolation='none')
    plt.title('Segmentation mask on all pixet1')
    file_name = os.path.join(path_python + '/figures_python/Seg_mask_all_pixet1.png')
    plt.savefig(file_name)

    #total number of pixels per image
    pixels_num = newLmod_arr.shape[0] * newLmod_arr.shape[0]
    print('Total number of pixels = ',pixels_num)
    print('Total number of pixels in pixet1 = ', sure_coords.shape[0])
    print('% of pixet1 versus total pixels = ', sure_coords.shape[0]*100/pixels_num)


    ###########identifying pixels per cell that are sure_coords
    pixet_1_list = []
    fill_in = np.full((1,2),-10)
    dataCells_pixet1 = np.zeros([valid_cells, valid_channelNum], dtype=np.float)

    print('identifying pixels per cell that are sure_coords and building count matrix based on pixet1 (=superpixet1)')
    for i, regions_i in enumerate(regions_revised):

        l = intersect2d(regions_i['coords'], sure_coords)
        if (len(l) != 0):
            dataCells_pixet1[i, :] = image_tensor[l[:, 0], l[:, 1], :].sum(axis=0)
            pixet_1_list.append(l)
        else:
            pixet_1_list.append(fill_in)

    # superpixet1
    df_pixet1_data = pd.DataFrame(dataCells_pixet1)
    df_pixet1_data.to_csv(
        os.path.join(path_python + '/superpixel_MIBI_python/count_data/count_matrix_6_updated.csv'),
        header=None, index=None)


    #pixet1_coords_list
    df_pixet1_list = pd.DataFrame(np.array(pixet_1_list))
    df_pixet1_list.to_csv(os.path.join(path_python + '/superpixel_MIBI_python/count_data/pixet_1_struct.csv'), header=None,
                   index=None)

    #clusterChannels
    df_cc = pd.DataFrame(clusterChannels)
    df_cc.to_csv(os.path.join(path_python + '/superpixel_MIBI_python/count_data/clusterChannels.csv'), header=None,
                   index=None)


    #meta_data - paths etc
    file_f = open(os.path.join(path_python + '/superpixel_MIBI_python/count_data/meta_data.txt'),"w")
    file_f.write(path_data+"\n")
    file_f.write(str(pointNumber) + "\n")
    file_f.write(str(image_tensor.shape[0]) + "\n")
    file_f.write(str(max_number_MCMC) + "\n")
    file_f.close()


    # load the image and convert it to a floating point data type

    cer_markers_tensor = np.zeros((newLmod_arr.shape[0], newLmod_arr.shape[1], valid_channelNum))
    for i in range(valid_channelNum):
        image_temp = np.array(Image.open(os.path.join(path_data + '/Point' + str(pointNumber) + '/' + clusterChannels[i] + '.tif')))
        cer_markers_tensor[:, :, i] = image_temp[0:patch_image, 0: patch_image]
        #sure_coords_temp = np.array(np.where(cer_markers_tensor[:, :, i] >= np.median(cer_markers_tensor[:, :, i].max(axis=0)) - threshold_pixels))
        #sure_coords = np.concatenate((sure_coords, sure_coords_temp.T),axis=0)
        certain_markers_all = certain_markers_all + cer_markers_tensor[:, :, i]


    image_all_markers = certain_markers_all #img_as_float(io.imread(args["image"]))
    image_all_markers[sure_coords[:,0],sure_coords[:,1]] = 0




    ####build the grid image using square tessellation

    height = newLmod_arr.shape[0]
    width = newLmod_arr.shape[1]
    grid_image = Image.new(mode='L', size=(height, width), color=255)
    blank_image = np.array(grid_image)
    # Draw some lines
    draw = ImageDraw.Draw(grid_image)
    y_start = 0
    y_end = grid_image.height
    step_size = int(grid_image.width / step_count)

    for x in range(0, grid_image.width, step_size):
        line = ((x, y_start), (x, y_end))
        draw.line(line, fill=128)

    x_start = 0
    x_end = grid_image.width

    for y in range(0, grid_image.height, step_size):
        line = ((x_start, y), (x_end, y))
        draw.line(line, fill=128)

    del draw

    grid_image.show()

    masked1 = np.ma.masked_where(grid_image == 255, grid_image)# 255 is white
    plt.figure()
    plt.imshow(blank_image, 'gray', interpolation='none')  # gray
    plt.imshow(masked1, cmap=plt.cm.cool, interpolation='none')
    plt.title('Grid as mask')
    file_name = os.path.join(path_python + '/figures_python/Grid_as_mask.png')
    plt.savefig(file_name)
    image2 = np.array(grid_image)  # masked1
    #gradient = sobel(rgb2gray(image2))

    #distance1 = ndi.distance_transform_edt(image2)
    ###local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((3,3)), labels=image_all_markers) # introduced the line below to reduce the WSS by half
    #local_maxi1 = peak_local_max(distance1, min_distance=1, indices=False, footprint=np.ones((3, 3)), labels=image2)
    #markers = ndi.label(local_maxi1)[0]
    #labels2 = watershed(gradient, markers, mask=gradient)





    plt.figure()

    plt.imshow(mark_boundaries(image_all_markers, masked1, color=(1, 0, 1)))
    plt.title('Segmentation mask/grid on all WSS')
    file_name = os.path.join(path_python + '/figures_python/Seg_grid_all_WSS.png')
    plt.savefig(file_name)

    plt.figure()
    plt.imshow(mark_boundaries(image_all_markers, masked1, color=(138, 40, 200)))
    plt.figure()
    plt.imshow(mark_boundaries(image_all_markers, masked1, color=(5, 140, 120)))

    # show the plots

    plt.show()

    label_sp2 = measure.label(image2, connectivity=neighbourhood_conn)
    regions_sp2 = measure.regionprops(label_sp2)
    sp2_segments = label_sp2.max() # number

    print('Number of WSS segments is ', sp2_segments)



    #clean up the WSS segments to just have non-zero pixel intensities
    start_neigh = time.time()
    regions_WSS_cleaned = []
    fill_in_0 =  np.full((1, 2), 0,dtype=np.int8)


    for i, regions_wss in enumerate(regions_sp2):
        WSS_temp = regions_wss['coords']
        if (len(WSS_temp) != 0):
            ##append only those pixel coords that have a rowsum of non-zero
            check_rows = image_tensor[WSS_temp[:, 0], WSS_temp[:, 1], :].sum(axis=1)
            label_cr = np.ones([1, len(WSS_temp)])
            ind_l = np.flatnonzero(check_rows <= 1)  # are removed. Any WSS segment as a whole is checked for its contribution.
            # If < 0.1, that whole segment will fall off. But those segments with just one rowsum being 0.1 or more will thrive,
            # putting even those pixel coords whose rowsum = 0, in that segment. So when assignment happens, that entire chunk is getting reassigned, even those without signal!
            # Therefore clean the segments before they come to this point.
            label_cr[:, ind_l] = 2
            label_good = np.flatnonzero(label_cr == 1)
            if (len(label_good) != 0):
                WSS_temp = WSS_temp[label_good]
                regions_WSS_cleaned.append(WSS_temp)  # (l.tolist()) each (x,y) coord is one cell in the .csv, else l is the np.ndarray and is written as so
            else:
                regions_WSS_cleaned.append(fill_in_0)

        else:
            regions_WSS_cleaned.append(fill_in_0)

    end_neigh = time.time()
    print('Time to clean superpixels: ' + str(end_neigh - start_neigh))



    #####clean up the superpixels from pixet1 (now superpixel1) to create superpixel2

    ############
    #### Building superpixels from the image - superpixet2 that need reassignment
    ###build the count matrix along with cell Area, centroids, scaled and unscaled
    ############

    data_WSS = np.zeros([sp2_segments, valid_channelNum], dtype=np.float)
    WSS_dataScaleSize = np.zeros([sp2_segments, valid_channelNum], dtype=np.float)
    WSS_Sizes = np.zeros([sp2_segments, 1], dtype=np.float)
    WSS_dataCentroids = np.zeros([sp2_segments, 2], dtype=np.float)
    WSS_dataScaleCentroids = np.zeros([sp2_segments, 2], dtype=np.float)

    '''
    for i, regions_i in enumerate(regions_sp2):  # check for pixel entries too
    print(i)
    l = regions_i['coords'] - sure_coords
    if(len(l)!=0):
        data_WSS[i, :] = image_tensor[l[:, 0], l[:, 1], :].sum(axis=0)
    '''

    start_neigh = time.time()
    neighbourhood_WSS_list = []
    neighbourhood_WSS_list.append(fill_in)
    pixet_2_list = []
    pixet_2_list.append(fill_in)

    print('building pixet2')
    #for i, regions_i in enumerate(regions_sp2):

    inputs = range(1,sp2_segments)
    num_cores = multiprocessing.cpu_count()
    Parallel(n_jobs=num_cores)(delayed(build_pixet2)(i) for i in inputs)

    # for i in range(1,sp2_segments):
    #     neigh_WSS_temp = []
    #     WSS_coords = regions_WSS_cleaned[i]
    #     l = setdiff2d(WSS_coords, sure_coords)
    #     if (len(l) != 0):
    #         print(i)
    #         pixet_2_list.append(l) #(l.tolist()) each (x,y) coord is one cell in the .csv, else l is the np.ndarray and is written as so
    #         data_WSS[i, :] = image_tensor[l[:, 0], l[:, 1], :].sum(axis=0)
    #         WSS_Sizes[i] = len(l)
    #         WSS_dataCentroids[i, :] = np.mean(l,axis=0)
    #         WSS_dataScaleSize[i, :] = data_WSS[i, :] / WSS_Sizes[i]
    #         WSS_dataScaleCentroids[i, :] = WSS_dataCentroids[i, :] / WSS_Sizes[i]
    #     else:
    #         pixet_2_list.append(fill_in)


    #     len_coords = len(WSS_coords)
    #     for j in range(len_coords):
    #         # print(j)
    #         neigh_WSS_temp.extend(cell_neighbors(label_sp2, WSS_coords[j, 0], WSS_coords[j, 1], dp))
    #         nt = np.unique(neigh_WSS_temp) - 1  # no idea why the cell labels are incremented by 1, check the cell_neighbourhood logic, for now I am just decrementing this
    #         nt = [y for y in nt if y != -1]  # remove -1. This means cell 0. a) real cell 0 will never be a neighbour b) why it is even there for every cell??
    #         nt = [y for y in nt if y != i]  # artefact of cell neighbourhppd function. Every cell is a neighbour of itself. So I remove this. But again the -1 from before needs to be tackled.
    #     neighbourhood_WSS_list.append(nt)

    end_neigh = time.time()
    print('Time to get superpixet2 and its neighbours for cells: ' + str(end_neigh - start_neigh))

    #clean up the WSS count matrix
    #get the final information only for segments that have
    #1. positive nuclear identity
    #2. enough information in the clustering channels to be clustered.


    labelIdentityNew21 = np.ones([1, sp2_segments])
    sumDataScaleSizeInClusterChannels21 = WSS_dataScaleSize.sum(axis=1)
    ind_1 = np.flatnonzero(sumDataScaleSizeInClusterChannels21 < 0.1) #are removed. Any WSS segment as a whole is checked for its contribution.
    # If < 0.1, that whole segment will fall off. But those segments with just one rowsum being 0.1 or more will thrive,
    # putting even those pixel coords whose rowsum = 0, in that segment. So when assignment happens, that entire chunk is getting reassigned, even those without signal!
    # Therefore clean the segments before they come to this point.
    labelIdentityNew21[:,ind_1] = 2
    labelVec21 = np.flatnonzero(labelIdentityNew21 == 1) # labels of valid segments
    data_WSS_segments = data_WSS[labelVec21,:] # CM based on sp2
    dataScaleSizeWSS = WSS_dataScaleSize[labelVec21,:]
    cellCentroids_WSS = WSS_dataCentroids[labelVec21,:]
    cellScaleCentroids_WSS = WSS_dataScaleCentroids[labelVec21,:]
    WSS_Sizes_revised = WSS_Sizes[labelVec21] #area:q

    valid_segments = len(labelVec21) #from 41030 to 27428 when divisor patch = 2

    pixet_2_list_revised = [y for x,y in enumerate(pixet_2_list) if x in labelVec21]

    #create a new regions structure with just the valid cells
    regions_sp2_revised = [y for x,y in enumerate(regions_WSS_cleaned) if x in labelVec21] #can also use contains(labelvec,x) #27428 when divisor patch = 2
    neighbourhood_WSS_list_revised = [y for x, y in enumerate(neighbourhood_WSS_list) if x in labelVec21]
    # remove those invalid superpixels from neighbourhoods as well
    for i, neigh_i in enumerate(neighbourhood_WSS_list_revised):
        neighbourhood_WSS_list_revised[i] = np.intersect1d(neigh_i, labelVec21)




    print('Sum of original Leeat count matrix is ', dataCells.sum(axis=0))
    print('Sum of pixet1 count matrix is ', dataCells_pixet1.sum(axis=0))
    print('Sum of WSS count matrix is ', data_WSS_segments.sum(axis=0))

    #assign a celllabel to WSS segments
    print('assigning cell id, if at all, to superpixet2')
    WSS_cell_label_assgn = np.ones([1, valid_segments])*-1
    inputs = enumerate(regions_sp2_revised)
    num_cores = multiprocessing.cpu_count()
    Parallel(n_jobs=num_cores)(delayed(assign_cell_id)(i, r_sp2) for i, r_sp2 in inputs)
    
    # for i, r_sp2 in enumerate(regions_sp2_revised):
    #     print("i "+str(i))
    #     for j, r in enumerate(regions_revised):
    #         l = intersect2d(regions_WSS_cleaned[i], r['coords'])
    #         if (len(l)!=0):
    #             #print('j'+str(j))
    #             WSS_cell_label_assgn[0,i]=labelVec[j]
    #             #print(labelVec[j])
    #             break






    ### writing patient image structures to .csv files for R intake


    #superpixet2
    df_WSS_segments_data = pd.DataFrame(data_WSS_segments)
    df_WSS_segments_data.to_csv(os.path.join(path_python + '/superpixel_MIBI_python/count_data/count_matrix_super_pixel.csv'), header=None,
                   index=None)
    #feather.write_dataframe(df_WSS_segments_data, os.path.join(path_python +'/superpixel_MIBI_python/count_data/count_matrix_super_pixel.feather'))

    #neighbourhood superpixet2
    df_WSS_neigh_list = pd.DataFrame(np.array(neighbourhood_WSS_list_revised))
    df_WSS_neigh_list.to_csv(os.path.join(path_python + '/superpixel_MIBI_python/count_data/neighbours_super_pixel.csv'), header=None,
                   index=None)
    #feather.write_dataframe(df_WSS_neigh_list, os.path.join(path_python +'/superpixel_MIBI_python/count_data/neighbours_super_pixel.feather'))



    #labelVec21 - valid WSS segments
    df_labelvec21 = pd.DataFrame(labelVec21)
    df_labelvec21.to_csv(os.path.join(path_python + '/superpixel_MIBI_python/count_data/valid_WSS.csv'), header=None,
                   index=None)
    #feather.write_dataframe(df_labelvec21, os.path.join(path_python +'/superpixel_MIBI_python/count_data/valid_WSS.feather'))



    #WSS to cell assignments
    df_WSS_cell_assgn = pd.DataFrame(WSS_cell_label_assgn)
    df_WSS_cell_assgn.to_csv(os.path.join(path_python + '/superpixel_MIBI_python/count_data/WSS_celllab_assgn.csv'), header=None,
                   index=None)
    #feather.write_dataframe(df_WSS_cell_assgn, os.path.join(path_python +'/superpixel_MIBI_python/count_data/WSS_celllab_assgn.feather'))


    #cleaned and corrected WSS area
    df_rWSSraa = pd.DataFrame(WSS_Sizes_revised)
    df_rWSSraa.to_csv(os.path.join(path_python + '/superpixel_MIBI_python/count_data/WSS_revised_area.csv'), header=None,
                   index=None)
    #feather.write_dataframe(df_rWSSraa, os.path.join(path_python +'/superpixel_MIBI_python/count_data/WSS_revised_area.feather'))

    #cleaned and corrected WSS centroid
    df_rWSSrcentt = pd.DataFrame(np.array(cellCentroids_WSS))
    df_rWSSrcentt.to_csv(os.path.join(path_python + '/superpixel_MIBI_python/count_data/WSS_revised_centroid.csv'), header=None,
                   index=None)
    #feather.write_dataframe(df_rWSSrcentt, os.path.join(path_python +'/superpixel_MIBI_python/count_data/WSS_revised_centroid.feather'))

    #rcleaned and corrected pixel coords
    df_rWSSrpcoordd = pd.DataFrame(np.array(pixet_2_list_revised))
    df_rWSSrpcoordd.to_csv(os.path.join(path_python + '/superpixel_MIBI_python/count_data/WSS_revised_pixelcoords.csv'), header=None,
                   index=None)
    #feather.write_dataframe(df_rWSSrpcoordd, os.path.join(path_python +'/superpixel_MIBI_python/count_data/WSS_revised_pixelcoords.feather'))



print('Overall time for code execution ', format(time.time()-t1) )


###extra code
'''
# ###neighbourhood detection per cell
start_neigh = time.time()
d=2
neighbourhood_list=[]

for i, regions_i in enumerate(regions):  # check for pixel entries too
    neigh_temp = []
    # coords_temp = regions[i]['coords'] # props[0].centroid
    # cellSizes[i] = regions[i]['Area']
    print(i)
    coords_temp = regions_i['coords']
    cellSizes[i] = regions_i['Area']
    dataCentroids[i, :] = regions_i['centroid']
    # for c in range(valid_channelNum):
    # y0,x0=props.centroid
    # print(regions[0]['area'])
    data[i, :] = image_tensor[coords_temp[:, 0], coords_temp[:, 1], :].sum(
        axis=0)  # sum the columns, per ith row, top down, length of vector
    dataScaleSize[i, :] = data[i, :] / cellSizes[i]
    dataScaleCentroids[i, :] = dataCentroids[i, :] / cellSizes[i]
    # find neighbours per region (aka cell)

    # plt.hold()
    plt.text(dataCentroids[i, 0], dataCentroids[i, 1], str(i), fontsize=7, color='red', verticalalignment='center',
             horizontalalignment='center')

    ''''''
    temp_arr = np.zeros([patch_image, patch_image])
    temp_arr[coords_temp] = 1
    im_dilate = ndi.binary_dilation(temp_arr, structure=se)  # .astype(temp_arr.dtype)
    neigh = np.logical_and(im_dilate, 1 - temp_arr)
    neigh_labels = np.unique(label[neigh])
    neigh_labels = neigh_labels[neigh_labels != 0]
    print('neighbouring cells for cell ' + str(i) + ':' + str(neigh_labels))
    neighbourhood_list.append(neigh_labels)
    # neighbourhood_list[str(i)] = neigh_labels
    # del neigh_labels
    ''''''
    len_coords = len(coords_temp)
    for j in range(len_coords):
        #print(j)
        neigh_temp.extend(cell_neighbors(label, coords_temp[j, 0], coords_temp[j, 1], d=2))
        nt = np.unique(neigh_temp)-1 #no idea why the cell labels are incremented by 1, check the cell_neighbourhood logic, for now I am just decrementing this
        nt = [y for y in nt if y != -1] #remove -1. This means cell 0. a) real cell 0 will never be a neighbour b) why it is even there for every cell??
        nt = [y for y in nt if y != i] #artefact of cell neighbourhppd function. Every cell is a neighbour of itself. So I remove this. But again the -1 from before needs to be tackled.
    neighbourhood_list.append(nt)
end_neigh = time.time()
print('Time to get neighbours for cells: '+ str(end_neigh - start_neigh))

###just checking the centroid plots
plt.figure()
plt.imshow(newLmod_arr * 200, 'gray', interpolation='none')  # gray
# plt.imshow(masked, cmap=plt.cm.cool, interpolation='none')
plt.title('Segmentation mask patch new')
for i, regions_i in enumerate(regions):  # check for pixel entries too

    print(i)
    coords_temp = regions_i['coords']
    cellSizes[i] = regions_i['Area']
    dataCentroids[i, :] = regions_i['centroid'] #np.mean(coords_temp,axis=0) #

    # plt.hold()
    plt.text(int(dataCentroids[i, 0]), int(dataCentroids[i, 1]), str(i), fontsize=7, color='red', verticalalignment='center',horizontalalignment='center')



start_neigh = time.time()

neighbourhood_WSS_list = []
    #for i, regions_i in enumerate(regions_sp2):
for i in range(10):  # check for pixel entries too
    neigh_WSS_temp = []
    #l = regions_sp2[i]['coords'] - sure_coords
    WSS_coords = regions_sp2[i]['coords']
    l = setdiff2d(WSS_coords, sure_coords)
    if (len(l) != 0):
            print(i)
            data_WSS[i, :] = image_tensor[l[:, 0], l[:, 1], :].sum(axis=0)
            WSS_Sizes[i] = len(l)
            WSS_dataCentroids[i, :] = np.mean(l,axis=0)
            # for c in range(valid_channelNum):
            # y0,x0=props.centroid
            # print(regions[0]['area'])
            #data_WSS[i, :] = image_tensor[WSS_coords_temp[:, 0], WSS_coords_temp[:, 1], :].sum(
            #    axis=0)  # sum the columns, per ith row, top down, length of vector
            WSS_dataScaleSize[i, :] = data_WSS[i, :] / WSS_Sizes[i]
            WSS_dataScaleCentroids[i, :] = WSS_dataCentroids[i, :] / WSS_Sizes[i]
            len_coords = len(WSS_coords)
            for j in range(len_coords):
        # print(j)
                neigh_WSS_temp.extend(cell_neighbors(label_sp2, WSS_coords[j, 0], WSS_coords[j, 1], dp))
                nt = np.unique(neigh_WSS_temp) - 1  # no idea why the cell labels are incremented by 1, check the cell_neighbourhood logic, for now I am just decrementing this
                nt = [y for y in nt if y != -1]  # remove -1. This means cell 0. a) real cell 0 will never be a neighbour b) why it is even there for every cell??
                nt = [y for y in nt if y != i]  # artefact of cell neighbourhppd function. Every cell is a neighbour of itself. So I remove this. But again the -1 from before needs to be tackled.
            neighbourhood_WSS_list.append(nt)
end_neigh = time.time()
print('Time to get neighbours for cells: ' + str(end_neigh - start_neigh))

#list creation
pixet_1_list = []
for i in range(10):  # check for pixel entries too
    WSS_coords = regions_sp2[i]['coords']
    l = setdiff2d(WSS_coords, sure_coords)
    if (len(l) != 0):
            print(i)
            pixet_1_list.append(l.tolist())

#break usage
for i in range(1, 5):
    for j in range(1, 5):
        if (i / j != 0):
            print(i)
            break

'''
#plt.close('all')

'''
gradient = sobel(rgb2gray(image_all_markers))

distance = ndi.distance_transform_edt(image_all_markers)
local_maxi = peak_local_max(distance, min_distance= 1, indices=False, footprint=np.ones((1, 1)), labels=image_all_markers)
markers = ndi.label(local_maxi)[0]
# segments_watershed = watershed(gradient, markers=250, compactness=0.001)
# labels = watershed(-distance, markers, mask=image_all_markers)
# labels1 = watershed(-gradient, markers, mask=image_all_markers)
labels2 = watershed(gradient, markers, mask=image_all_markers)

plt.figure()
# plt.imshow(labels2, cmap=plt.cm.nipy_spectral, interpolation='nearest')
# plt.imshow(masked, cmap=plt.cm.cool, interpolation='none')
plt.imshow(mark_boundaries(image_all_markers, labels2, color=(1, 0, 1)))
plt.imshow(masked, cmap=plt.cm.cool, interpolation='none')
# show the plots
plt.show()

label_sp2 = measure.label(labels2, connectivity=neighbourhood_conn)
regions_sp2 = measure.regionprops(label_sp2)
sp2_segments = label_sp2.max()  # number
'''