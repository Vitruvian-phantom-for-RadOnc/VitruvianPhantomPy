#
#
# Djamal Boukerroui, 08-2022

import PIL as pil
import numpy as np
import skimage.io as skio
from scipy import ndimage
import os
import matplotlib.pyplot as plt


def create_random_noise(blob_shape, scale=4000):

    # blob_data = np.empty(blob_shape, np.uint16)
    blob_data = np.random.randint(0, scale, blob_shape, np.uint16)
    return blob_data


def create_image_vitruvians(blob_shape, ind=0, w_shape_pixels=200):

    vitruvius_info = [
        {
            "w": 268,
            "centerS": [0, 0],
            "imfile": "collage-2231082_640_corrected.png",
            "centerSquare": [318., 318.5],   # swapped  and corrected for matlab origin (-1)
        },
        {
            "w": 289,
            "centerS": [0, 0],
            # http://architectura.cesr.univ-tours.fr/Traite/Images/BPNME276Index.asp
            "imfile": "CollPriv_Pry_E276_0113_512m.png",
            "centerSquare": [313., 222.5],  # swapped  and corrected for matlab origin (-1)
        }
        ]

    pathfile = os.path.join(os.path.dirname(__file__), vitruvius_info[ind]['imfile'])
    img = skio.imread(pathfile, as_gray=True)

    if len(w_shape_pixels) == 1:
        w_shape_pixels = [w_shape_pixels, w_shape_pixels]

    scaleX = vitruvius_info[ind]['w'] / w_shape_pixels[1]  # Columns
    scaleY = vitruvius_info[ind]['w'] / w_shape_pixels[0]  # Rows

    dim = np.array(blob_shape[0:2])
    center_square = np.array(vitruvius_info[ind]['centerSquare'])
    offsetY = center_square[0] - scaleY * (dim[0] - 1) / 2
    offsetX = center_square[1] - scaleX * (dim[1] - 1) / 2
    Ts = np.array([[scaleY, 0], [0, scaleX]])

    mode = ('constant', 'mirror', 'wrap', 'reflect')
    imo = ndimage.affine_transform(img, Ts, offset=[offsetY, offsetX], output_shape=dim, prefilter=False,
                                       order=5, mode=mode[2])
    imo = imo*255 + 1024
    blob_data = imo
    for k in range(1, blob_shape[2]):
        blob_data = np.dstack((blob_data, imo + 10*(-1)**k))
    blob_data = blob_data.astype(np.uint16)

    return blob_data


def plot_shapes(point_dataset, name_list=None):

    for index in range(0, len(point_dataset)):
        current_contour = point_dataset[index]
        points_first_slice = current_contour['Contours'][0]
        np_tmp = np.array(points_first_slice)
        # reshape to get x and y coordinates
        np_tmp = np_tmp.reshape((3, -1), order='F')

        plt.figure(index)
        plt.plot(np_tmp[0], np_tmp[1], 'r+')
        plt.show()
    return


def transform_2d_shapes(shapes_list, transform):
    """
    A wrapper to  transform_2d_contour for all contours

    :param shapes_list: list of 3D shapes
    :param transform: dictionary with fields
        type in ('Rotation', 'Translation', 'Matrix') and
        data  corresponding data: angle value, translation vector or 2x2 matrix
    :return: transformed shapes
    """
    for index in range(0, len(shapes_list)):
        current_contour = shapes_list[index]
        coordinates = transform_2d_contour(current_contour['Contours'], transform)
        shapes_list[index]['Contours'] = coordinates
    return shapes_list


def transform_2d_contour(coordinates, transform):
    """
    Perform a 2D transformation of a list 3D points coordinates

    :param coordinates: a list of 3D points (contours)
    :param transform: dictionary with fields
        type in ('Rotation', 'Translation', 'Matrix') and
        data  corresponding data: angle value, translation vector or 2x2 matrix
    :return: transformed coordinates
    """
    for i in range(0, len(coordinates)):
        np_tmp = np.array(coordinates[i])
        # reshape to get x and y coordinates
        np_tmp = np_tmp.reshape((3, -1), order='F')
        point_xy = transform_2d_points(np_tmp[0:2, :], transform)
        np_tmp[0:2, :] = point_xy
        # put back the transformed points
        coordinates[i] = np_tmp.flatten(order='F')
    return coordinates


def transform_2d_points(input_xy, transform):
    """
    Perform a 2D transformation of a set of input points

    :param input_xy: 2xn Numpy array
    :param transform: dictionary with fileds
        type in ('Rotation', 'Translation', 'Matrix') and
        data  corresponding data: angle value, translation vector or 2x2 matrix
    :return: transformed points
    """
    value = transform['data']

    if transform['type'] == 'Rotation':
        rot = np.array([[np.math.cos(value), -np.math.sin(value)],[np.math.sin(value), np.math.cos(value)]])
        output_xy = rot @ input_xy
    elif transform['type'] == 'Translation':
        output_xy = input_xy + np.array(value).reshape(1,2)
    elif transform['type'] == 'Matrix':
        output_xy = value @ input_xy
    else:
        output_xy = None

    return output_xy


