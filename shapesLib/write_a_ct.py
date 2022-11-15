#
#
# Mark Gooding, 08, 2022
# Revisions
#
# Djamal Boukerroui, 11,2022
#       Add SeriesNumber Tag: Required, Empty if Unknown (2)
#       Minor correctionsafter running DICOM Validation tool 3.0

import pydicom
from pydicom.dataset import Dataset
from pydicom import uid
from pydicom.sequence import Sequence
from datetime import datetime
import numpy as np
import os


def write_a_ct_volume(fullfilepath, blob, spacing, origin, patient_name=None, patient_id=None, study_desc=None, series_desc=None):

    dims = blob.shape
    study_instance_uid = uid.generate_uid()
    series_instance_uid = uid.generate_uid()
    frame_or_ref_uid = uid.generate_uid()  # TODO might need to pass back or receive in for RTSS

    for i in range(0, dims[2]):
        sop_instance_uid = uid.generate_uid()
        filename = os.path.join(fullfilepath, '{:03d}.dcm'.format(i))
        uids_to_pass = [study_instance_uid, series_instance_uid, sop_instance_uid, frame_or_ref_uid]
        slice_position = [origin[0], origin[1], origin[2]+i*spacing[2]]
        write_a_ct_slice(filename, uids_to_pass, blob[:, :, i].copy(order='C'), spacing, slice_position, patient_name=patient_name,
                         patient_id=patient_id, slice_number=i, study_desc=study_desc, series_desc=series_desc)

    # end of loop


def write_a_ct_slice(fullfilepathandname, uid_list, blob_slice, spacing, slice_position, patient_name=None,
                     patient_id=None, slice_number=None, study_desc=None, series_desc=None):

    currently = datetime.now()

    # File meta info data elements
    file_meta = Dataset()
    file_meta.FileMetaInformationGroupLength = 196
    file_meta.FileMetaInformationVersion = b'\x00\x01'
    file_meta.MediaStorageSOPClassUID = '1.2.840.10008.5.1.4.1.1.2'
    file_meta.MediaStorageSOPInstanceUID = uid_list[2]
    file_meta.TransferSyntaxUID = '1.2.840.10008.1.2.1'
    file_meta.ImplementationClassUID = '1.2.40.0.13.1.1.1'
    file_meta.ImplementationVersionName = 'pydicomVitruvius'

    # Main data elements
    ds = Dataset()
    ds.ImageType = 'DERIVED\\PRIMARY\\AXIAL'
    ds.SpecificCharacterSet = 'ISO_IR 100'
    ds.InstanceCreationDate = currently.strftime("%Y%m%d")
    ds.InstanceCreationTime = currently.strftime("%H%M%S.%f")
    ds.SOPClassUID = '1.2.840.10008.5.1.4.1.1.2'
    ds.SOPInstanceUID = uid_list[2]
    ds.StudyDate = currently.strftime("%Y%m%d")
    ds.StudyTime = currently.strftime("%H%M%S.%f")
    ds.AccessionNumber = ''
    ds.Modality = 'CT'
    ds.Manufacturer = 'Mirada Medical Ltd'
    ds.ReferringPhysicianName = ''
    if study_desc is not None:
        ds.StudyDescription = study_desc
    else:
        ds.StudyDescription = 'Synthetic Shape Test Set'
    if study_desc is not None:
        ds.SeriesDescription = series_desc
    else:
        ds.SeriesDescription = 'CT for Synthetic Shape Test Set'
    if patient_name is not None:
        ds.PatientName = patient_name
    else:
        ds.PatientName = 'Synthetic'
    if patient_id is not None:
        ds.PatientID = patient_id
    else:
        ds.PatientID = 'Synthetic'
    ds.PatientBirthDate = ''
    ds.PatientSex = 'O'
    ds.SliceThickness = '{:.6f}'.format(spacing[2])
    ds.KVP = ''
    ds.PatientPosition = 'HFS'
    ds.StudyInstanceUID = uid_list[0]
    ds.StudyID = ''
    ds.SeriesInstanceUID = uid_list[1]
    ds.SeriesNumber = '1'
    ds.AcquisitionNumber = ''
    if slice_number is None:
        ds.InstanceNumber = ''
    else:
        ds.InstanceNumber = '{:d}'.format(slice_number)
    ds.ImagePositionPatient = ['{:.6f}'.format(slice_position[0]), '{:.6f}'.format(slice_position[1]),
                               '{:.6f}'.format(slice_position[2])]
    ds.ImageOrientationPatient = ['1.000000', '0.000000', '0.000000', '0.000000', '1.000000', '0.000000']
    ds.FrameOfReferenceUID = uid_list[3]
    ds.PositionReferenceIndicator = ''
    ds.SliceLocation = '{:.6f}'.format(slice_position[2])
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = 'MONOCHROME2'
    ds.Rows = blob_slice.shape[0]
    ds.Columns = blob_slice.shape[1]
    ds.PixelSpacing = ['{:.6f}'.format(spacing[0]), '{:.6f}'.format(spacing[1])]
    ds.BitsAllocated = 16
    ds.BitsStored = 16
    ds.HighBit = 15
    ds.PixelRepresentation = 1
    ds.RescaleIntercept = "-1024.000000"
    ds.RescaleSlope = "1.000000"
    ds.RescaleType = 'HU'
    ds.PixelData = blob_slice.tobytes()

    ds.file_meta = file_meta
    ds.is_implicit_VR = False
    ds.is_little_endian = True
    ds.save_as(fullfilepathandname, write_like_original=False)
