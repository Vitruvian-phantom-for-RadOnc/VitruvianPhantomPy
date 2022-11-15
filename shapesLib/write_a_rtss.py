#
#
# Mark Gooding, 08, 2022
# Revisions
#
# Djamal Boukerroui, 11,2022
#       Minor corrections after running DICOM Validation tool 3.0

# Coded version of DICOM file 'IM1.DCM'
# Produced by pydicom codify utility script
from __future__ import unicode_literals  # Only for python2.7 and save_as unicode filename
import pydicom
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from pydicom import uid
from datetime import datetime
import numpy as np
import os


def read_ct_dir(input_directory):
    non_dicom_ignored_files = 0
    ct_data_list = []
    for curr_dir, sub_dirs, files in os.walk(input_directory):
        if len(files) > 0:
            # Simple check:
            for testfile in files:
                full_testfile_name = os.path.join(curr_dir, testfile)
                try:
                    ct_data = pydicom.dcmread(full_testfile_name)
                except pydicom.InvalidDicomError:
                    # Not a dicom file. Ignore
                    non_dicom_ignored_files = non_dicom_ignored_files + 1
                    continue

                if ct_data.SOPClassUID == '1.2.840.10008.5.1.4.1.1.2':
                    # This is one of the CTs, so add to the pile
                    ct_data_list.append(ct_data)

    #print('{} ct files read'.format(len(ct_data_list)))

    return ct_data_list


def write_a_rtss(fullfilepath, ct_data_location, point_data, filename=None, series_desc=None, series_num=None):

    currently = datetime.now()

    if filename is None:
        fullfilepathandname = os.path.join(fullfilepath, 'IM1.dcm')
    else:
        fullfilepathandname = os.path.join(fullfilepath, filename)

    if type(ct_data_location) == list:
        ct_data = ct_data_location # assume input is the read CT
    elif os.path.isdir(ct_data_location):
        ct_data = read_ct_dir(ct_data_location)
    else:
        raise RuntimeError("No CT or not created yet")

    series_instance_uid = uid.generate_uid()
    sop_instance_uid = uid.generate_uid()

    # File meta info data elements
    file_meta = Dataset()
    file_meta.FileMetaInformationGroupLength = 204
    file_meta.FileMetaInformationVersion = b'\x00\x01'
    file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3'
    file_meta.MediaStorageSOPInstanceUID = sop_instance_uid
    file_meta.TransferSyntaxUID = '1.2.840.10008.1.2.1'
    file_meta.ImplementationClassUID = '1.2.826.0.1.3680043.8.691.0.50'
    file_meta.ImplementationVersionName = 'pydicomVitruvius'

    # Main data elements
    ds = Dataset()
    ds.InstanceCreationDate = currently.strftime("%Y%m%d")
    ds.InstanceCreationTime = currently.strftime("%H%M%S.%f")
    ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.481.3'
    ds.SOPInstanceUID = sop_instance_uid
    ds.StudyDate = ct_data[0].StudyDate
    ds.StudyTime = ct_data[0].StudyTime
    ds.AccessionNumber = ''
    ds.Modality = 'RTSTRUCT'
    ds.Manufacturer = 'Mirada Medical Ltd'
    ds.ReferringPhysicianName = 'Vitruvius'
    ds.StationName = os.environ['COMPUTERNAME']
    ds.StudyDescription = ct_data[0].StudyDescription
    if series_desc is None:
        ds.SeriesDescription = 'RTSS for Synthetic Shape Test Set'
    else:
        ds.SeriesDescription = series_desc
    ds.PatientName = ct_data[0].PatientName
    ds.PatientID = ct_data[0].PatientID
    ds.PatientBirthDate = ''
    ds.PatientSex = 'O'
    ds.PatientAge = ''
    ds.StudyInstanceUID = ct_data[0].StudyInstanceUID
    ds.SeriesInstanceUID = series_instance_uid
    ds.StudyID = '1'
    if series_num is None:
        ds.SeriesNumber = '1'
    else:
        ds.SeriesNumber = series_num
    ds.InstanceNumber = "1"
    ds.StructureSetLabel = 'TESTDATA'
    ds.StructureSetName = 'TESTDATA'
    if series_desc is None:
        ds.StructureSetDescription = 'TESTDATA'
    else:
        ds.StructureSetDescription = series_desc
    ds.StructureSetDate = currently.strftime("%Y%m%d")
    ds.StructureSetTime = currently.strftime("%H%M%S.%f")

    # Referenced Frame of Reference Sequence
    refd_frame_of_ref_sequence = Sequence()
    ds.ReferencedFrameOfReferenceSequence = refd_frame_of_ref_sequence

    # Referenced Frame of Reference Sequence: Referenced Frame of Reference
    refd_frame_of_ref = Dataset()
    refd_frame_of_ref.FrameOfReferenceUID = ct_data[0].FrameOfReferenceUID

    # RT Referenced Study Sequence
    rt_refd_study_sequence = Sequence()
    refd_frame_of_ref.RTReferencedStudySequence = rt_refd_study_sequence

    # RT Referenced Study Sequence: RT Referenced Study
    rt_refd_study = Dataset()
    rt_refd_study.ReferencedSOPClassUID = '1.2.840.10008.3.1.2.3.2'  # not sure what this relates too
    rt_refd_study.ReferencedSOPInstanceUID = ct_data[0].StudyInstanceUID

    # RT Referenced Series Sequence
    rt_refd_series_sequence = Sequence()
    rt_refd_study.RTReferencedSeriesSequence = rt_refd_series_sequence

    # RT Referenced Series Sequence: RT Referenced Series
    rt_refd_series = Dataset()
    rt_refd_series.SeriesInstanceUID = ct_data[0].SeriesInstanceUID

    # Contour Image Sequence
    contour_image_sequence = Sequence()
    rt_refd_series.ContourImageSequence = contour_image_sequence

    # Contour Image Sequence for each CT image
    for i in range(0, len(ct_data)):
        contour_image = Dataset()
        contour_image.ReferencedSOPClassUID = ct_data[i].SOPClassUID
        contour_image.ReferencedSOPInstanceUID = ct_data[i].SOPInstanceUID
        contour_image_sequence.append(contour_image)

    rt_refd_series_sequence.append(rt_refd_series)
    rt_refd_study_sequence.append(rt_refd_study)
    refd_frame_of_ref_sequence.append(refd_frame_of_ref)

    # Structure Set ROI Sequence
    structure_set_roi_sequence = Sequence()
    ds.StructureSetROISequence = structure_set_roi_sequence

    # Structure Set ROI Sequence: for each point set
    for i in range(0, len(point_data)):
        structure_set_roi = Dataset()
        structure_set_roi.ROINumber = '{:d}'.format(i+1)
        structure_set_roi.ReferencedFrameOfReferenceUID = ct_data[0].FrameOfReferenceUID
        structure_set_roi.ROIName = point_data[i]['Name']
        structure_set_roi.ROIGenerationAlgorithm = 'AUTOMATIC'
        structure_set_roi.ROIGenerationDescription = 'Synthetic analytical shapes'
        structure_set_roi_sequence.append(structure_set_roi)

    # ROI Contour Sequence
    roi_contour_sequence = Sequence()
    ds.ROIContourSequence = roi_contour_sequence

    # ROI Contour Sequence: for each ROI Contour
    for i in range(0, len(point_data)):

        roi_contour = Dataset()
        roi_contour.ROIDisplayColor = point_data[i]['Color']

        # Contour Sequence
        contour_sequence = Sequence()
        roi_contour.ContourSequence = contour_sequence

        # Contour Sequence: For each slice/Contour
        for j in range(0, len(point_data[i]['Contours'])):
            contour = Dataset()

            # Contour Image Sequence
            contour_image_sequence = Sequence()
            contour.ContourImageSequence = contour_image_sequence

            # Contour Image Sequence: Contour Image 1
            contour_image = Dataset()
            contour_image.ReferencedSOPClassUID = ct_data[0].SOPClassUID
            # need to find a ct slice corresponding to the z-location of this slice
            uid_to_use = None
            for k in range(0, len(ct_data)):
                if point_data[i]['Contours'][j][2] == ct_data[k].ImagePositionPatient[2]:
                    uid_to_use = ct_data[k].SOPInstanceUID
                    continue
            if uid_to_use is None:
                print('CT for contour not found')
                continue
            contour_image.ReferencedSOPInstanceUID = uid_to_use
            contour_image_sequence.append(contour_image)

            contour.ContourGeometricType = 'CLOSED_PLANAR'
            contour.NumberOfContourPoints = '{:d}'.format(int(len(point_data[i]['Contours'][j])/3))
            contour.ContourNumber = '{:d}'.format(j+1)
            contour_as_string = []
            for k in range(0, len(point_data[i]['Contours'][j])):
                contour_as_string.append('{:f}'.format(point_data[i]['Contours'][j][k]))
            contour.ContourData = contour_as_string
            contour_sequence.append(contour)
        roi_contour.ReferencedROINumber = i+1
        roi_contour_sequence.append(roi_contour)

    # RT ROI Observations Sequence
    rtroi_observations_sequence = Sequence()
    ds.RTROIObservationsSequence = rtroi_observations_sequence

    # RT ROI Observations Sequence: for each RT ROI Observations
    for i in range(0, len(point_data)):
        rtroi_observations = Dataset()
        rtroi_observations.ObservationNumber = '{:d}'.format(i+1)
        rtroi_observations.ReferencedROINumber = '{:d}'.format(i+1)
        rtroi_observations.RTROIInterpretedType = 'ORGAN'
        rtroi_observations.ROIInterpreter = ''
        rtroi_observations_sequence.append(rtroi_observations)

    ds.ApprovalStatus = 'UNAPPROVED'

    ds.file_meta = file_meta
    ds.is_implicit_VR = False
    ds.is_little_endian = True
    ds.save_as(fullfilepathandname, write_like_original=False)
