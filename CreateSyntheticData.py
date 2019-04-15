import vtk
from vtk.util.misc import vtkGetDataRoot
VTK_DATA_ROOT = vtkGetDataRoot()
import sys
import os
import argparse
import numpy as np
import csv

def main(args):

    resolution            = args.resolution
    nb_groups             = args.nb_groups
    nb_timepoints         = args.nb_timepoints
    nb_subjects           = args.nb_subjects
    base_stratch          = args.base_stratch
    base_rotation         = args.base_rotation
    stretch_factor        = args.stretch_factor
    rotation_factor       = args.rotation_factor
    outputfile            = args.outputfile

    # use a sphere as a basis of the shape
    sphere = vtk.vtkSphereSource()
    sphere.SetPhiResolution(resolution)
    sphere.SetThetaResolution(resolution)
    sphere.Update()
    sphereData = sphere.GetOutput()
    writer = vtk.vtkPolyDataWriter()
    writer.SetFileName('Shape_Tamplate.vtk')
    writer.SetInputData(sphereData)
    writer.Write()

    # create a data array to hold the weighting coefficients
    tfarray = vtk.vtkFloatArray()
    npoints = sphereData.GetNumberOfPoints()
    tfarray.SetNumberOfComponents(2)
    tfarray.SetNumberOfTuples(npoints)

    # parameterize the sphere along the z axis, and fill the weights
    # with (1.0-a, a) to linearly interpolate across the shape
    i = 0
    while i < npoints:
        pt = sphereData.GetPoint(i)
        x = pt[0]
        y = pt[1]
        z = pt[2]
        zn = z + 0.5
        zn1 = 1.0 - zn
        if (zn > 1.0):
            zn = 1.0
        if (zn1 < 0.0):
            zn1 = 0.0
        tfarray.SetComponent(i, 0, zn1)
        tfarray.SetComponent(i, 1, zn)
        i += 1

    # create field data to hold the array, and bind it to the sphere
    fd = vtk.vtkFieldData()
    tfarray.SetName("weights")
    sphereData.GetPointData().AddArray(tfarray)

    Shapes = []
    Cov = np.zeros((nb_groups*nb_timepoints*nb_subjects,2))
    idx = 0

    for g in range(nb_groups):

        for t in range(nb_timepoints):

            for s in range(nb_subjects):

                stretch_noise  = np.random.normal(0, 0.05, 1)
                rotation_noise = np.random.uniform(-2,2,1)

                # use an ordinary transform to stretch the shape
                stretch = vtk.vtkTransform()
                stretch.Scale(1, 1, base_stratch + stretch_factor * (t+1) + stretch_noise)
                stretchFilter = vtk.vtkTransformFilter()
                stretchFilter.SetInputData(sphereData)
                stretchFilter.SetTransform(stretch)
                stretchFilter.Update()
                stretchedSphereData = stretchFilter.GetOutput()

                # now, for the weighted transform stuff
                weightedTrans = vtk.vtkWeightedTransformFilter()
                # create two transforms to interpolate between
                identity = vtk.vtkTransform()
                identity.Identity()
                rotated = vtk.vtkTransform()

                rotated.RotateX(base_rotation + rotation_factor*(g+1)*(t+1) + rotation_noise)
                weightedTrans.SetNumberOfTransforms(2)
                weightedTrans.SetTransform(identity, 0)
                weightedTrans.SetTransform(rotated, 1)

                weightedTrans.SetWeightArray("weights")
                weightedTrans.SetInputConnection(stretchFilter.GetOutputPort())
                weightedTrans.Update()
                TransStretchedSphereData = weightedTrans.GetOutput()

                outputfilename = './Shapes/' + outputfile + '_' + str(g) + '_' + str(t) + '_' + str(s) + '.vtk'
                # save stretched shape
                writer = vtk.vtkPolyDataWriter()
                writer.SetFileName(outputfilename)
                writer.SetInputData(TransStretchedSphereData)
                writer.Write()

                Shapes.append(outputfilename)
                Cov[idx,0] = int(t)
                Cov[idx,1] = int(g)
                idx+=1

    with open('Shapes.txt', 'w') as f:
        for item in Shapes:
            f.write("%s\n" % item)

    np.savetxt("Cov.csv", Cov, delimiter=",")


def parse_args(argv):   
    parser = argparse.ArgumentParser()

    parser.add_argument("--resolution", type=int, help="Sphere Resolution", default=40)
    parser.add_argument("--nb_groups", type=int, help="Number of Groups", default=3)
    parser.add_argument("--nb_timepoints", type=int, help="Number of Time Points", default=3)
    parser.add_argument("--nb_subjects", type=int, help="Number of Subjects", default=30)
    parser.add_argument("--base_stratch", type=int, help="Initial Stretch", default=3)
    parser.add_argument("--base_rotation", type=int, help="Initial Rotation", default=30)
    parser.add_argument("--stretch_factor", type=float, help="Progression Rotation Factor", default=0.5)
    parser.add_argument("--rotation_factor", type=float, help="Progression Rotation Factor", default=5)
    parser.add_argument("--outputfile", help="Output file base", default='Shape')

    return parser.parse_args(argv)

if __name__ == "__main__":
    main(parse_args(sys.argv[1:]))

