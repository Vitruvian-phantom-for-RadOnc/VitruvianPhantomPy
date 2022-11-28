#
#
# (c) Djamal Boukerroui, 11-2022

import argparse
import json
import textwrap

from shapesLib import write_a_ct as ctfile
from shapesLib import write_a_rtss as rtssfile
from shapesLib.shape_utils import create_image_vitruvians
import shapesLib.paper_shapes as sp
import numpy as np
import os



def save_shapes(ctfilepath, rtssfilepath, pointdataset, split=True):
    if split:
        ct_data = rtssfile.read_ct_dir(ctfilepath)
        for index, current in enumerate(pointdataset):
            name = current['Name']
            rtssfile.write_a_rtss(rtssfilepath, ct_data, [current], filename=name.replace(' ', '') + '.dcm',
                                  series_desc=name, series_num=index)
    else:
        rtssfile.write_a_rtss(rtssfilepath, ctfilepath, pointdataset, filename='allshapes.dcm',
                              series_desc='allshapes', series_num=1)


# Leonardo da Vinci Vitruvian
def create_ldv_vitruvian_cts(fullfilepath, spacing, blob_shape, origin,
                         study_desc=None, patient_name=None, patient_id=None,
                         series_desc=None, w_shape_pixels=160.):
    # this will create the data volume for the RTSS to live in

    # blob_data = np.empty(blob_shape, np.uint16)
    blob_data = create_image_vitruvians(blob_shape, 0, [w_shape_pixels/spacing[0], w_shape_pixels/spacing[1]])

    ctfile.write_a_ct_volume(fullfilepath, blob_data, spacing, origin, study_desc=study_desc,
                             patient_name=patient_name, patient_id=patient_id, series_desc=series_desc)

# Cesare Cesariano Vitruvian
def create_cc_vitruvian_cts(fullfilepath, spacing, blob_shape, origin,
                         study_desc=None, patient_name=None, patient_id=None,
                         series_desc=None, w_shape_pixels=160.):
    # this will create the data volume for the RTSS to live in

    # blob_data = np.empty(blob_shape, np.uint16)
    blob_data = create_image_vitruvians(blob_shape, 1, [w_shape_pixels/spacing[0], w_shape_pixels/spacing[1]])

    ctfile.write_a_ct_volume(fullfilepath, blob_data, spacing, origin, study_desc=study_desc,
                             patient_name=patient_name, patient_id=patient_id, series_desc=series_desc)

def create_cts(fullfilepath, spacing, blob_shape, origin, study_desc=None, patient_name=None, patient_id=None, series_desc=None):
    # this will create the data volume for the RTSS to live in

    # blob_data = np.empty(blob_shape, np.uint16)
    blob_data = np.random.randint(0, 4000, blob_shape, np.uint16)

    ctfile.write_a_ct_volume(fullfilepath, blob_data, spacing, origin, study_desc=study_desc,
                             patient_name=patient_name, patient_id=patient_id, series_desc=series_desc)

def usage():
    help_text = textwrap.dedent('''
        Generates data for :  
        Analytic calculations and synthetic shapes for validation of quantitative contour comparison software''')
    parser = argparse.ArgumentParser(description=help_text, formatter_class=argparse.RawTextHelpFormatter)

    help_text = 'Full path to output directory to save results'
    parser.add_argument('output_dir', type=str, help=help_text)

    help_text = 'Full path to json file specifying Cylinder Cuboid  shapes to create'
    parser.add_argument('-j', '--shapes_json', required=False,
                        help=help_text)  # type=str, this is not required as it's the default
    parser.add_argument('-c', '--cuboids', action='store_true', required=False, help='To create the Cuboid pair shapes')

    parser.add_argument('-n', '--noCT', action='store_true', required=False,
                        help='To not create a CT series if present')
    parser.add_argument('-nr', '--do_Typ', action='store_true', required=False, help='Process typical resolution')
    parser.add_argument('-hr', '--do_HR', action='store_true', required=False, help='Process high resolution')
    parser.add_argument('-lr', '--do_LR', action='store_true', required=False, help='Process low resolution')

    return parser.parse_args()


def load_dict_shapes(shapes_json):
    with open(shapes_json) as f:
        current = json.load(f)

    CC_VitruvianMan_shapes_list = LdV_VitruvianMan_shapes_list = None
    if 'CC_VitruvianMan' in current:
        CC_VitruvianMan_shapes_list = current.get('CC_VitruvianMan')
    if 'LdV_VitruvianMan' in current:
        LdV_VitruvianMan_shapes_list = current.get('LdV_VitruvianMan')

    return CC_VitruvianMan_shapes_list, LdV_VitruvianMan_shapes_list


def main():
    args = usage()

    shift_function = lambda x: 0.  # No shift, all shapes are centred at the origin
    if args.shapes_json != None:
        print('\nReading json file:', args.shapes_json)
        CC_VitruvianMan_shapes_list, LdV_VitruvianMan_shapes_list = load_dict_shapes(args.shapes_json)
        shift_function = lambda x: 80. - x  # The align the shapes at the inferior bound
    else:
        CC_VitruvianMan_shapes_list = LdV_VitruvianMan_shapes_list = None

    create_ct = not (args.noCT)
    create_rtss = True

    # # Common Tags to all resolution
    patient_id = patient_name = "VitruvianMan"
    if LdV_VitruvianMan_shapes_list != None:
        series_desc = "LdV_VitruvianMan"
        my_create_cts = create_ldv_vitruvian_cts
    else:
        series_desc = "CC_VitruvianMan"
        my_create_cts = create_cc_vitruvian_cts

    # # High Resolution
    if args.do_HR:
        print('\nProcessing High Resolution geometry: Fine ')
        try:
            spacing = [0.5, 0.5, 0.5]
            #blob_shape = [1024, 1024, 500]
            blob_shape = [512, 512, 350]
            origin = [-spacing[1] * (blob_shape[0] - 1) / 2, -spacing[0] * (blob_shape[1] - 1) / 2,
                      -spacing[2] * (blob_shape[2] - 1) / 2]  # origin at voxel centers

            if LdV_VitruvianMan_shapes_list != None:
                ctfilepath = os.path.join(args.output_dir, 'LdV_Hi-Res/CT')
                rtssfilepath = os.path.join(args.output_dir, 'LdV_Hi-Res/RTSS')
                study_desc = "LdV_VitruvianMan Hi-Res"
            else:
                ctfilepath = os.path.join(args.output_dir, 'CC_Hi-Res/CT')
                rtssfilepath = os.path.join(args.output_dir, 'CC_Hi-Res/RTSS')
                study_desc = "CC_VitruvianMan Hi-Res"

            if os.path.isdir(ctfilepath) == False:
                os.makedirs(ctfilepath)
            if os.path.isdir(rtssfilepath) == False:
                os.makedirs(rtssfilepath)

                # Create/Write a CT if needed
            if create_ct:
                print('\tGeneration of the CT data')
                my_create_cts(ctfilepath, spacing, blob_shape, origin, study_desc=study_desc,
                              patient_name=patient_name, patient_id=patient_id, series_desc=series_desc)
                # Create and save into RTSS files
            if create_rtss:
                if CC_VitruvianMan_shapes_list != None:
                    print('\tGeneration of Cesare Cesariano Vitruvian Man shapes')
                    shapes = sp.create_cylinder_cuboid_shapes(spacing, blob_shape, origin, CC_VitruvianMan_shapes_list)
                    save_shapes(ctfilepath, rtssfilepath, shapes, True)

                if LdV_VitruvianMan_shapes_list != None:
                    print('\tGeneration of Leonardo da Vinci Vitruvian Man shapes')
                    shapes = sp.create_vitruvian_shapes(spacing, blob_shape, origin, LdV_VitruvianMan_shapes_list,
                                                        shift_fun=shift_function)
                    save_shapes(ctfilepath, rtssfilepath, shapes, True)
                if args.cuboids:
                    print('\tGeneration of Cube vs shifted Cuboid shapes')
                    shapes = sp.create_cube_cuboid_shapes(spacing, blob_shape, origin)
                    save_shapes(ctfilepath, rtssfilepath, shapes, True)
        except Exception as ME:
            print(ME.message, flush=True)
            print("Something went wrong for the high resolution")

     #  sotropic
    if args.do_Typ:
        print('\nProcessing Typical resolution geometry: Normal')
        try:
            spacing = [0.96, 0.96, 1.0]
            blob_shape = [512, 512, 250]
            origin = [-spacing[1] * (blob_shape[0] - 1) / 2, -spacing[0] * (blob_shape[1] - 1) / 2,
                      -spacing[2] * (blob_shape[2] - 1) / 2]  # origin at voxel centers

            if LdV_VitruvianMan_shapes_list != None:
                ctfilepath = os.path.join(args.output_dir, 'LdV_Typical/CT')
                rtssfilepath = os.path.join(args.output_dir, 'LdV_Typical/RTSS')
                study_desc = "LdV_VitruvianMan Typical"
            else:
                ctfilepath = os.path.join(args.output_dir, 'CC_Typical/CT')
                rtssfilepath = os.path.join(args.output_dir, 'CC_Typical/RTSS')
                study_desc = "CC_VitruvianMan Typical"

            if os.path.isdir(ctfilepath) == False:
                os.makedirs(ctfilepath)
            if os.path.isdir(rtssfilepath) == False:
                os.makedirs(rtssfilepath)

            # Create/Write a CT if needed
            if create_ct:
                print('\tGeneration of the CT data')
                my_create_cts(ctfilepath, spacing, blob_shape, origin, study_desc=study_desc,
                              patient_name=patient_name, patient_id=patient_id, series_desc=series_desc)
            # Create and save into RTSS files
            if create_rtss:
                if CC_VitruvianMan_shapes_list != None:
                    print('\tGeneration of Cesare Cesariano Vitruvian Man shapes')
                    shapes = sp.create_cylinder_cuboid_shapes(spacing, blob_shape, origin, CC_VitruvianMan_shapes_list)
                    save_shapes(ctfilepath, rtssfilepath, shapes, True)

                if LdV_VitruvianMan_shapes_list != None:
                    print('\tGeneration of Leonardo da Vinci Vitruvian Man shapes')
                    shapes = sp.create_vitruvian_shapes(spacing, blob_shape, origin, LdV_VitruvianMan_shapes_list, shift_fun=shift_function)
                    save_shapes(ctfilepath, rtssfilepath, shapes, True)
                if args.cuboids:
                    print('\tGeneration of Cube vs shifted Cuboid shapes')
                    shapes = sp.create_cube_cuboid_shapes(spacing, blob_shape, origin)
                    save_shapes(ctfilepath, rtssfilepath, shapes, True)
        except Exception as ME:
            print(ME.message, flush=True)
            print("Something went wrong for the isotropic resolution")

    # Coarse
    if args.do_LR:
        print('\nProcessing lower resolution geometry: Coarse')
        try:
            spacing = [1.8, 2.2, 3.0]  # [Rows, columns, slice] [ Y X Z]
            blob_shape = [256, 256, 84]
            origin = [-spacing[1] * (blob_shape[0] - 1) / 2, -spacing[0] * (blob_shape[1] - 1) / 2,
                      -spacing[2] * (blob_shape[2] - 1) / 2]  # origin at voxel centers
            if LdV_VitruvianMan_shapes_list != None:
                ctfilepath = os.path.join(args.output_dir, 'LdV_Coarse/CT')
                rtssfilepath = os.path.join(args.output_dir, 'LdV_Coarse/RTSS')
                study_desc = "LdV_VitruvianMan Coarse"
            else:
                ctfilepath = os.path.join(args.output_dir, 'CC_Coarse/CT')
                rtssfilepath = os.path.join(args.output_dir, 'CC_Coarse/RTSS')
                study_desc = "CC_VitruvianMan Coarse"

            if os.path.isdir(ctfilepath) == False:
                os.makedirs(ctfilepath)
            if os.path.isdir(rtssfilepath) == False:
                os.makedirs(rtssfilepath)

                # Create/Write a CT if needed
            if create_ct:
                print('\tGeneration of the CT data')
                my_create_cts(ctfilepath, spacing, blob_shape, origin, study_desc=study_desc,
                              patient_name=patient_name, patient_id=patient_id, series_desc=series_desc)
                # Create and save into RTSS files
            if create_rtss:
                if CC_VitruvianMan_shapes_list != None:
                    print('\tGeneration of Cesare Cesariano Vitruvian Man shapes')
                    shapes = sp.create_cylinder_cuboid_shapes(spacing, blob_shape, origin, CC_VitruvianMan_shapes_list)
                    save_shapes(ctfilepath, rtssfilepath, shapes, True)

                if LdV_VitruvianMan_shapes_list != None:
                    print('\tGeneration of Leonardo da Vinci Vitruvian Man shapes')
                    shapes = sp.create_vitruvian_shapes(spacing, blob_shape, origin, LdV_VitruvianMan_shapes_list,
                                                        shift_fun=shift_function)
                    save_shapes(ctfilepath, rtssfilepath, shapes, True)
                if args.cuboids:
                    print('\tGeneration of Cube vs shifted Cuboid shapes')
                    shapes = sp.create_cube_cuboid_shapes(spacing, blob_shape, origin)
                    save_shapes(ctfilepath, rtssfilepath, shapes, True)
        except Exception as ME:
            print(ME.message, flush=True)
            print("Something went wrong for the coarse resolution")


if __name__ == '__main__':
    main()
