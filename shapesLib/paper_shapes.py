#
#
# Djamal Boukerroui, 08-2022

from .shape_factory import make_cuboid_shape, make_cylinder_shape
from .shape_utils import transform_2d_contour
import numpy as np

def create_cube_cuboid_shapes(spacing, blob_shape, origin):
    point_data = []

    # Reference Cube - Sparse Sampled
    # ----------
    # 10 cm square in plane
    # 5 cm high
    # centred at 0,0,0
    # corners only
    contourset = {}
    contourset['Name'] = 'Reference Cuboid Sparse'
    contourset['Color'] = ['0', '255', '0']
    contourset['Contours'] = make_cuboid_shape([-50.0, -50.0, -25.0, 50.0, 50.0, 25.0],
                                               spacing, blob_shape, origin, 'corners')
    point_data.append(contourset)

    # Test Cuboid - Sparse Sampled
    # ----------
    # 10x11 cm in plane
    # 5 cm high
    # centred at 31.1,	21.2,	11.3
    # corners only
    contourset = {}
    contourset['Name'] = 'Test Cuboid Sparse'
    contourset['Color'] = ['255', '0', '0']
    contourset['Contours'] = make_cuboid_shape([-18.9, -33.8, -13.7, 81.1, 76.2, 36.3],
                                               spacing, blob_shape, origin, 'corners')
    point_data.append(contourset)

    # Test Cuboid - Dense Sampled
    # ----------
    # 10x11 cm in plane
    # 5 cm high
    # centred at 31.1,	21.2,	11.3
    # corners only
    contourset = {}
    contourset['Name'] = 'Test Cuboid Dense'
    contourset['Color'] = ['255', '0', '0']
    contourset['Contours'] = make_cuboid_shape([-18.9, -33.8, -13.7, 81.1, 76.2, 36.3],
                                               spacing, blob_shape, origin, 'pixels')
    point_data.append(contourset)
    return point_data


def create_cylinder_cuboid_shapes(spacing, blob_shape, origin, shape_list):
    """
    Create the cylinder and cuboid shapes are specified in the shape_list

    :param spacing: voxel spacing of the CT
    :param blob_shape: size of the CT
    :param origin: origin of the CT
    :param shape_list:  a list of dictionaries with parameters shape
    :return: list of shapes ready to be saved in an RTSS file
    """
    point_data = []

    for current in shape_list:
        # corners only cuboid of length 2w
        bounds = (-current['w'], -current['w'], -current['w'], current['w'], current['w'], current['w'])
        current_contour = make_cuboid_shape(bounds, spacing, blob_shape, origin, 'corners')
        # perform 2d  axial rotation if required
        if current['theta'] != 0.0:
            transform = {'type': 'Rotation', 'data': current['theta'] * np.math.pi / 180.0}
            current_contour = transform_2d_contour(current_contour, transform)
        # Store in shape list
        contourset = {}
        contourset['Name'] = current['nameA']
        contourset['Color'] = ['0', '255', '0']
        contourset['Contours'] = current_contour
        point_data.append(contourset)

        # Cylinder
        current_contour = make_cylinder_shape(current['r'], [0.0, 0.0], [-current['w'], current['w']], origin,
                                                       spacing, blob_shape)
        if current['theta'] != 0.0:
            transform = {'type': 'Rotation', 'data': current['theta']}
            current_contour = transform_2d_contour(current_contour, transform)
        # Store in shape list
        contourset = {}
        contourset['Name'] = current['nameB']
        contourset['Color'] = ['0', '255', '255']
        contourset['Contours'] = current_contour
        point_data.append(contourset)
    return point_data


def create_vitruvian_shapes(spacing, blob_shape, origin, shape_list, centerS=True, shift_fun=lambda x:0):
    """
    Create the Vitruvian man shape are specified in the shape_list

    :param shift_fun:
    :param centerS: origin of parametrisation of the shape (square if true)
    :param spacing: voxel spacing of the CT
    :param blob_shape: size of the CT
    :param origin: origin of the CT
    :param shape_list:  a list of dictionaries with parameters shape
    :return: list of shapes ready to be saved in an RTSS file
    """
    point_data = []

    for current in shape_list:
        # corners only cuboid of length 2w
        shift_y = current['w'] - current['r']
        if centerS:
            center_square = [0.0, shift_fun(current['w'])]
            circle_origin = [0.0, center_square[1] + shift_y]
        else:
            center_square = [0.0, shift_fun(current['w']) - shift_y]
            circle_origin = [0.0, shift_fun(current['w'])]

        bounds = (center_square[0] - current['w'], center_square[1] - current['w'], -current['w'],
                  center_square[0] + current['w'], center_square[1] + current['w'], current['w'])

        current_contour = make_cuboid_shape(bounds, spacing, blob_shape, origin, 'corners')
        # perform 2d  axial rotation if required
        if current['theta'] != 0.0:
            transform = {'type': 'Rotation', 'data': current['theta'] * np.math.pi / 180.0}
            current_contour = transform_2d_contour(current_contour, transform)
        # Store in shape list
        contourset = {'Name': current['nameA'], 'Color': ['0', '255', '0'], 'Contours': current_contour}
        point_data.append(contourset)

        # Cylinder
        current_contour = make_cylinder_shape(current['r'], circle_origin, [-current['w'], current['w']],
                                                       origin, spacing, blob_shape)
        if current['theta'] != 0.0:
            transform = {'type': 'Rotation', 'data': current['theta']}
            current_contour = transform_2d_contour(current_contour, transform)
        # Store in shape list
        contourset = {'Name': current['nameB'], 'Color': ['0', '255', '255'], 'Contours': current_contour}
        point_data.append(contourset)
    return point_data
