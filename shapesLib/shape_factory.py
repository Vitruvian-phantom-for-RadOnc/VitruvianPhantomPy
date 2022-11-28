import numpy as np
import math
import os
#
#
# Djamal Boukerroui & Mark Gooding
# 08-2022



def make_cuboid_shape(bounds, spacing, blob_shape, origin, sampling='pixels'):
    '''
    Make a cuboid shape
    :param bounds: coordinates of  cuboid boundaries [x0,y0,z0, x1,y1,z1]
    :param spacing:
    :param blob_shape:
    :param origin:
    :param sampling: pixels or corners only sampling
    :return:
    '''
    contours_list = []

    for k in range(0, blob_shape[2]):
        z = origin[2] + spacing[2] * k
        contour = []

        if sampling == 'pixels':
            if (z > bounds[2]) & (z <= bounds[5]):
                for i in range(0, blob_shape[0] + 1, 1):
                    y = bounds[4]
                    x = (i - 0.5) * spacing[0] + origin[0]
                    if (x > bounds[0]) & (x <= bounds[3]):
                        contour.append(x)
                        contour.append(y)
                        contour.append(z)
                for j in range(blob_shape[1] + 1, 0, -1):
                    y = (j - 0.5) * spacing[1] + origin[1]
                    x = bounds[3]
                    if (y > bounds[1]) & (y <= bounds[4]):
                        contour.append(x)
                        contour.append(y)
                        contour.append(z)
                for i2 in range(blob_shape[0] + 1, 0, -1):
                    y = bounds[1]
                    x = (i2 - 0.5) * spacing[0] + origin[0]
                    if (x > bounds[0]) & (x <= bounds[3]):
                        contour.append(x)
                        contour.append(y)
                        contour.append(z)
                for j2 in range(0, blob_shape[1] + 1, 1):
                    y = (j2 - 0.5) * spacing[1] + origin[1]
                    x = bounds[0]
                    if (y > bounds[1]) & (y <= bounds[4]):
                        contour.append(x)
                        contour.append(y)
                        contour.append(z)
        else:
            if (z > bounds[2]) & (z <= bounds[5]):
                contour = [bounds[0], bounds[4], z,
                           bounds[3], bounds[4], z,
                           bounds[3], bounds[1], z,
                           bounds[0], bounds[1], z]
        if len(contour) > 0:
            contours_list.append(contour)

    return contours_list


def make_cylinder_shape(radius, centre, zbounds,
                        origin, spacing, blob_shape):
    """
    Make a cylinder shape as polygons to be saved in an RTSS

    :param radius: radius of the cylinder shape
    :param centre: [x, y] coordinates of the center
    :param zbounds:[z_start, z_end] slice bounds
    :param origin: [x, y, z] coordinates of the field of view (CT)
    :param spacing:[x, y, z] voxel spacing of the
    :param blob_shape: [dimx, dimy, dimz] of the blob.
    :return: contour list
    """
    contours_list = []
    vlines = np.arange(origin[0] - .5,
                       origin[0] + (blob_shape[0] - 0.5) * spacing[0],
                       spacing[0])
    hlines = np.arange(origin[1] - .5,
                       origin[1] + (blob_shape[1] - 0.5) * spacing[1],
                       spacing[1])

    for k in range(0, blob_shape[2]):
        z = origin[2] + spacing[2] * k
        if (z <= zbounds[0]) | (z > zbounds[1]):
            continue

        contour = []
        pts = make_circle_as_polygon(centre, radius, hlines, vlines)
        if len(pts) == 0:
            continue
        for value in pts:
            contour.append(value[1])
            contour.append(value[2])
            contour.append(z)

        contours_list.append(contour)
    return contours_list


# General version
def make_lnorm_shape(radius, centre, zbounds,
                     origin, spacing, blob_shape, order=1.0):
    """
    Make a L-norm shape as polygons to be saved in an RTSS

    :param radius: radius of the cylinder shape
    :param centre: [x, y] coordinates of the center
    :param zbounds:[z_start, z_end] slice bounds
    :param origin: [x, y, z] coordinates of the field of view (CT)
    :param spacing:[x, y, z] voxel spacing of the
    :param blob_shape: [dimx, dimy, dimz] of the blob.
    :param order: order of the L-norm shape
    :return: contour list
    """
    contours_list = []
    vlines = np.arange(origin[0] - .5 * spacing[0],
                       origin[0] + (blob_shape[0] - 0.5) * spacing[0],
                       spacing[0])
    hlines = np.arange(origin[1] - .5 * spacing[1],
                       origin[1] + (blob_shape[1] - 0.5) * spacing[1],
                       spacing[1])

    for k in range(0, blob_shape[2]):
        z = origin[2] + spacing[2] * k
        if (z <= zbounds[0]) | (z > zbounds[1]):
            continue

        contour = []
        pts = make_lnorm_as_polygon(centre, radius, hlines, vlines, order)
        if len(pts) == 0:
            continue
        for value in pts:
            contour.append(value[1])
            contour.append(value[2])
            contour.append(z)

        contours_list.append(contour)
    return contours_list


def make_circle_as_polygon(centre, radius, Hlines, VLines):
    """
    compute the intersection of a circle with a grid and sort the points as a polygon

    :param centre: [x, y] center coordinates
    :param radius: of the circle
    :param Hlines: a list of horizontal lines
    :param VLines: a list of vertical lines
    :return: a list of sorted intersection point (x,y)
    """
    pts = []
    for y in Hlines:
        tmp = intersect_circle_hline(centre[0], centre[1], radius, y)
        pts = pts + tmp

    for x in VLines:
        tmp = intersect_circle_vline(centre[0], centre[1], radius, x)
        pts = pts + tmp
    # Sort points using angle: This is ok only for convex shapes
    pts.sort(key=lambda elem: elem[0])
    return pts


def make_lnorm_as_polygon(centre, radius, Hlines, VLines, order=1.0):
    """
    compute the intersection of a l-norm shape with a grid and sort the points as a polygon

    :param centre: [x, y] center coordinates
    :param radius: of the circle
    :param Hlines: a list of horizontal lines
    :param VLines: a list of vertical lines
    :param order: oder of the L-norm shape
    :return: a list of sorted intersection point (x,y)
    """
    pts = []
    for y in Hlines:
        tmp = intersect_lnorm_hline(centre[0], centre[1], radius, y, order)
        pts = pts + tmp

    for x in VLines:
        tmp = intersect_lnorm_vline(centre[0], centre[1], radius, x, order)
        pts = pts + tmp
    # Sort points using angle: This is ok only for convex shapes
    pts.sort(key=lambda elem: elem[0])
    return pts


# slicewise functions

def intersect_circle_hline(x0, y0, radius, cst):
    """
    Compute the intersection of horizontal line y=cst with a circle

    :param x0: x coordinate of  the center
    :param y0: y coordinate of the center
    :param radius: radius of the of the L-norm shape
    :param cst: horizontal line coefficient
    :return: list of intersection points with angle [(theta, x, y)]
    """
    # line y = cst
    pts = []
    d = math.fabs(cst - y0)
    if d > radius:
        return pts
    if d == radius:
        x = x0
        y = cst
        theta = math.atan2((x - x0) / radius, (y - y0) / radius)
        pts.append((theta, x, y))
    else:
        delta = math.sqrt(radius * radius - d * d)
        x = x0 - delta
        y = cst
        theta = math.atan2((x - x0) / radius, (y - y0) / radius)
        pts.append((theta, x, y))

        x = x0 + delta
        y = cst
        theta = math.atan2((x - x0) / radius, (y - y0) / radius)
        pts.append((theta, x, y))
    return pts


def intersect_lnorm_hline(x0, y0, radius, cst, order=1.0):
    """
    Compute the intersection of horizontal line y=cst with the
    L-norm shape

    :param x0: x coordinate of  the center
    :param y0: y coordinate of the center
    :param radius: radius of the of the L-norm shape
    :param cst: horizontal line coefficient
    :param order: order of the of the L-norm shape
    :return: list of intersection points with angle [(theta, x, y)]
    """
    # line y = cst
    pts = []
    d = math.fabs(cst - y0)
    if d > radius:
        return pts
    if d == radius:
        x = x0
        y = cst
        theta = math.atan2((x - x0) / radius, (y - y0) / radius)
        pts.append((theta, x, y))
    else:
        delta = math.pow(radius ** order - d ** order, 1.0 / order)
        x = x0 - delta
        y = cst
        theta = math.atan2((x - x0) / radius, (y - y0) / radius)
        pts.append((theta, x, y))

        x = x0 + delta
        y = cst
        theta = math.atan2((x - x0) / radius, (y - y0) / radius)
        pts.append((theta, x, y))
    return pts


def intersect_circle_vline(x0, y0, radius, cst):
    """
    Compute the intersection of vertical line x=cst with a circle

    :param x0: x coordinate of  the center
    :param y0: y coordinate of the center
    :param radius: radius of the of the L-norm shape
    :param cst: vertical line coefficient
    :return: list of intersection points with angle [(theta, x, y)]
    """
    # line x = cst
    pts = []
    d = math.fabs(cst - x0)
    if d > radius:
        return pts
    if d == radius:
        x = cst
        y = y0
        theta = math.atan2((x - x0) / radius, (y - y0) / radius)
        pts.append((theta, x, y))
    else:
        delta = math.sqrt(radius * radius - d * d)
        x = cst
        y = y0 - delta
        theta = math.atan2((x - x0) / radius, (y - y0) / radius)
        pts.append((theta, x, y))

        x = cst
        y = y0 + delta
        theta = math.atan2((x - x0) / radius, (y - y0) / radius)
        pts.append((theta, x, y))
    return pts


def intersect_lnorm_vline(x0, y0, radius, cst, order=1.0):
    """
        Compute the intersection of vertical line x=cst with a L-norm shape

        :param x0: x coordinate of  the center
        :param y0: y coordinate of the center
        :param radius: radius of the of the L-norm shape
        :param cst: vertical line coefficient
        :param order: order of the of the L-norm shape
        :return: list of intersection points with angle [(theta, x, y)]
        """
    # line x = cst
    pts = []
    d = math.fabs(cst - x0)
    if d > radius:
        return pts
    if d == radius:
        x = cst
        y = y0
        theta = math.atan2((x - x0) / radius, (y - y0) / radius)
        pts.append((theta, x, y))
    else:
        delta = math.pow(radius ** order - d ** order, 1.0 / order)
        x = cst
        y = y0 - delta
        theta = math.atan2((x - x0) / radius, (y - y0) / radius)
        pts.append((theta, x, y))

        x = cst
        y = y0 + delta
        theta = math.atan2((x - x0) / radius, (y - y0) / radius)
        pts.append((theta, x, y))
    return pts
