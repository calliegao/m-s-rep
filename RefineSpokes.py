#!/usr/bin/env python
from __future__ import print_function
import numpy as np
import math
import vtk
from vtk.util.numpy_support import vtk_to_numpy
import pyvista as pv
# noinspection PyUnresolvedReferences
import vtkmodules.vtkInteractionStyle
# noinspection PyUnresolvedReferences
import vtkmodules.vtkRenderingOpenGL2
from vtkmodules.vtkCommonCore import (
    mutable,
    vtkPoints
)
from vtkmodules.vtkCommonDataModel import vtkPolygon
import os
import subprocess
import sys
import argparse
import re

def IsInsideCheck(pX,pY,pZ,mesh):
    select = vtk.vtkSelectEnclosedPoints()
    select.SetSurfaceData(mesh)
    select.SetTolerance(1e-4)
    
    pts = vtk.vtkPoints()
    pts.InsertNextPoint((pX),(pY),(pZ))
    pts_pd = vtk.vtkPolyData()
    pts_pd.SetPoints(pts)
    select.SetInputData(pts_pd)
    select.Update()
    return select.IsInside(0)    

def CalculateNormalVectorofIntersection(pt,surf_vtk):
    #Calculate normal vectors of the surface
    normFilter = vtk.vtkPolyDataNormals()
    normFilter.SetInputData(surf_vtk)
    normFilter.SetComputePointNormals(1)
    normFilter.SetComputeCellNormals(0)
    normFilter.SetAutoOrientNormals(1)
    normFilter.SetSplitting(0)
    normFilter.Update()
    NormalVector = normFilter.GetOutput()#已经是VTK文件      
    all_Normals = np.array(NormalVector.GetPointData().GetNormals())
    #Calculate the closest point of tp on the boundary surface
    cell_locator = vtk.vtkCellLocator()
    cell_locator.SetDataSet(surf_vtk)
    cell_locator.BuildLocator()
    cellId = vtk.reference(0)
    c = [0.0, 0.0, 0.0]
    subId = vtk.reference(0)
    d = vtk.reference(0.0)
    cell_locator.FindClosestPoint(pt, c, cellId, subId, d)
    pt_ids = vtk.vtkIdList()
    surf_vtk.GetCellPoints(cellId, pt_ids)
    num_cell = pt_ids.GetNumberOfIds()
    #print(num_cell)
    vector_array = np.zeros((num_cell, 3))
    num_cell = 1
    p_closestPoints = np.zeros((num_cell, 3))
    for i in range(num_cell):
        pt_id = pt_ids.GetId(i)
        #Calculate normal vector of each point
        p_closestPoints[i,:] = np.array(surf_vtk.GetPoint(pt_id))    
        n1 = all_Normals[pt_id,:]
        vector_array[i, :] = n1
        mean_vector = np.mean(vector_array,axis=0)
        nomalVector = mean_vector / np.linalg.norm(mean_vector)
    p_closestPoint = p_closestPoints[0,:]
    return p_closestPoint, nomalVector

def IntersectionNumber(pt,ps,surf_vtk):
    #Compute intersection number of spoke and surface 
    lamda = 0.05
    spoke = pt - ps
    spoke_length = np.linalg.norm(spoke)
    spoke_dir = spoke / spoke_length
    spoke_dir_new = spoke_dir
    spoke_length_new = 0
    status = 0
    cycle_num = 0
    Last_IsInside = IsInsideCheck(pt[0],pt[1],pt[2],surf_vtk)
    while np.dot(spoke_dir_new,spoke_dir)>0:#spoke_length_new<spoke_length:
          cycle_num += 1
          spoke_length = np.linalg.norm(pt - ps) 
          spoke_dir = (pt - ps)/spoke_length
          #print('Checking intersections: cycle_num = %s' % str(cycle_num))
          pt = pt - lamda*spoke_dir
          spoke_length_new = np.linalg.norm(pt - ps) 
          spoke_dir_new = (pt - ps)/spoke_length_new
          IsInside = IsInsideCheck(pt[0],pt[1],pt[2],surf_vtk)
          #print('Decreasing spoke length: spoke_length = %s' % str(spoke_length_new))     
          if IsInside != Last_IsInside:
             status += 1
          Last_IsInside = IsInside
    return status    

def RefineSpokeLength(surf_vtk,pt_vtk,ps_vtk,eps_s):#surface, tips points, skeleton points, error
    num_pt = pt_vtk.GetNumberOfPoints()
    pt_array = np.zeros((num_pt, 3))    
    num_ps = ps_vtk.GetNumberOfPoints()
    ps_array = np.zeros((num_ps, 3))
    finished_points = 0
    #p = pv.Plotter()
    #p.add_mesh(surf_vtk, color='orange', opacity=0.4)
    for i in range(num_pt):
        finished_points += 1
        pt = np.array(pt_vtk.GetPoint(i))
        ps = np.array(ps_vtk.GetPoint(i))
        if np.linalg.norm(pt-ps)!=0:
           spoke = pt - ps
           spoke_length = np.linalg.norm(spoke)
           spoke_dir = spoke / spoke_length
           #whether pt_array[i, :] is inside the surface
           #如果有一个交点，缩短spoke使得pt刚好在物体内部（上一个步长在外部）；如果有0个交点，延长spoke使得pt刚好在物体内部（下一个步长在内部）；如果有两个交点，缩短spoke使得只有一个交点，然后继续缩短直到pt刚好在物体内部（上一个步长在外部）；如果有三个交点，缩短spoke使得只有两个交点，然后继续缩短使得只有一个交点，然后继续缩短直到pt刚好在物体内部（上一个步长在外部）
           Intersect_num = IntersectionNumber(pt,ps,surf_vtk)
           IsInside = IsInsideCheck(pt[0],pt[1],pt[2],surf_vtk)
           pt_original = pt
           if Intersect_num == 0:
              #print('Before reparing spoke length: %s' % str(spoke_length))
              #print('0 intersection.')
              while IsInside==1:
                    pt = pt+spoke_dir*eps_s
                    IsInside = IsInsideCheck(pt[0],pt[1],pt[2],surf_vtk)
                    spoke_length = np.linalg.norm(pt-ps)
                    #print('Intersect_num == 0: Reparing spoke length: %s' % str(spoke_length))
                    if spoke_length>15:
                         pt = pt_original
                         break
           elif Intersect_num == 1:
                #print('1 intersection.')
                while IsInside==0:
                      pt = pt-spoke_dir*eps_s
                      IsInside = IsInsideCheck(pt[0],pt[1],pt[2],surf_vtk)
                      spoke_length = np.linalg.norm(pt-ps)
                      #print('Intersect_num == 1: Reparing spoke length: %s' % str(spoke_length))
                      if spoke_length>15:
                         pt  = pt_original
                      break
           elif Intersect_num == 2:
                #print('2 intersections.')
                while Intersect_num==1:
                      pt = pt-spoke_dir*eps_s
                      Intersect_num = IntersectionNumber(pt,ps,surf_vtk)
                      spoke_length = np.linalg.norm(pt-ps)
                      #print('Intersect_num == 2: Reparing spoke length Step 1: %s' % str(spoke_length))
                      if spoke_length>15:
                         pt = pt_original
                      break
                while IsInside==0:
                      pt = pt-spoke_dir*eps_s
                      IsInside = IsInsideCheck(pt[0],pt[1],pt[2],surf_vtk)
                      spoke_length = np.linalg.norm(pt-ps)
                      #print('Intersect_num == 2: Reparing spoke length Step 2: %s' % str(spoke_length))
                      if spoke_length>15:
                         pt = pt_original
                      break
           elif Intersect_num == 3:
                while Intersect_num==2:
                      pt = pt-spoke_dir*eps_s
                      Intersect_num = IntersectionNumber(pt,ps,surf_vtk)
                      spoke_length = np.linalg.norm(pt-ps)
                      #print('Intersect_num == 3: Reparing spoke length Step 1: %s' % str(spoke_length))
                      if spoke_length>15:
                         pt = pt_original
                      break
                while Intersect_num==1:
                      pt = pt-spoke_dir*eps_s
                      Intersect_num = IntersectionNumber(pt,ps,surf_vtk)
                      spoke_length = np.linalg.norm(pt-ps)
                      #print('Intersect_num == 3: Reparing spoke length Step 2: %s' % str(spoke_length))
                      if spoke_length>15:
                         pt = pt_original
                      break
                while IsInside==0:
                      pt = pt-spoke_dir*eps_s
                      IsInside = IsInsideCheck(pt[0],pt[1],pt[2],surf_vtk)
                      spoke_length = np.linalg.norm(pt-ps)
                      #print('Intersect_num == 3: Reparing spoke length Step 3: %s' % str(spoke_length))
                      if spoke_length>15:
                         pt = pt_original
                      break
                #print('2 intersections.')
        pt_array[i, :] = pt     
        #print('Refine Spoke Length Finished Points: %s' % str(finished_points))
    pts = vtk.vtkPoints()
    #print(pt_array)
    for j in range(num_pt):
            pts.InsertNextPoint((pt_array[j, 0]),(pt_array[j, 1]),(pt_array[j, 2]))
    pt_vtk_LengthRefined = vtk.vtkPolyData()
    pt_vtk_LengthRefined.SetPoints(pts)  
    return pt_vtk_LengthRefined

def RefineSpokeDirection(surf_vtk,pt_vtk,ps_vtk,eps_d):    
    num_pt = pt_vtk.GetNumberOfPoints()
    pt_array = np.zeros((num_pt, 3))    
    num_ps = ps_vtk.GetNumberOfPoints()
    ps_array = np.zeros((num_ps, 3))
    #p = pv.Plotter()
    #p.add_mesh(surf_vtk, color='orange', opacity=0.4)
    #num_pt = 100
    finished_points = 0
    alpha = 0.5 #alpha越大越接近normalvector
    for i in range(num_pt):
        finished_points += 1
        pt_array[i, :] = np.array(pt_vtk.GetPoint(i))
        ps_array[i, :] = np.array(ps_vtk.GetPoint(i))
        pt = pt_array[i, :]
        ps = ps_array[i, :]
        if np.linalg.norm(pt-ps)!=0:
           spoke = pt - ps
           spoke_length = np.linalg.norm(spoke)
           spoke_dir = spoke / spoke_length
           #eps_d = 0
           cos_angle = 0
           circle_num = 0
           while (1-cos_angle)>eps_d and circle_num<50: 
                 circle_num += 1
           #perpendicular to boundary surface & equal angle with the counterpart spoke
                 #Calculate the normal vector at the intersection on the surface
                 p_closestPoint, nomalVector = CalculateNormalVectorofIntersection(pt,surf_vtk)
                 #p.add_mesh(pv.Arrow(ps, spoke, scale=spoke_length), color='blue')
                 #p.add_mesh(pv.Arrow(p_closestPoint, nomalVector, scale=1), color='green')
                 #Update the spoke direction as the normal vector's
                 spoke = pt - ps
                 spoke_length = np.linalg.norm(spoke)
                 spoke_dir = spoke / spoke_length
                 #print(nomalVector)
                 medial_vector = (alpha*spoke_dir + (1-alpha)*nomalVector)
                 mv = medial_vector/np.linalg.norm(medial_vector)
                 pt = ps + spoke_length*mv
                 spoke = pt - ps
                 spoke_length = np.linalg.norm(spoke)
                 spoke_dir = spoke / spoke_length
                 cos_angle = np.dot(spoke_dir,nomalVector)
                 #print('Cos_angle = %s' % str(cos_angle))
                 #print(pt)
                 #p.add_mesh(pv.Arrow(ps, spoke, scale=spoke_length), color='black')
        #p.add_mesh(pv.Arrow(ps, spoke, scale=spoke_length), color='red')
        pt_array[i, :] = pt
        #print('Refine Spoke Direction Finished Points: %s' % str(finished_points))
    pts = vtk.vtkPoints()
    #print(pt_array)
    for i in range(num_pt):
        pts.InsertNextPoint((pt_array[i, 0]),(pt_array[i, 1]),(pt_array[i, 2]))
    pt_vtk_DirectionRefined = vtk.vtkPolyData()
    pt_vtk_DirectionRefined.SetPoints(pts)      
    #p.show() 
    return pt_vtk_DirectionRefined

def EqualSpokeLength(pt_vtk,ps_vtk,eps_e):#tips points, skeleton points, error
    num_pt = pt_vtk.GetNumberOfPoints()
    pt_array = np.zeros((num_pt, 3))    
    num_ps = ps_vtk.GetNumberOfPoints()
    ps_array = np.zeros((num_ps, 3))
    InnerSpokeLength = np.zeros((549, 2))#((549, 2))
    finished_points = 0
    for i in range(num_pt):
        pt_array[i,:] = np.array(pt_vtk.GetPoint(i))
        ps_array[i,:] = np.array(ps_vtk.GetPoint(i))
        pt = pt_array[i,:] 
        ps = ps_array[i,:]
        if np.linalg.norm(pt-ps)!=0:
           spoke = pt - ps
           spoke_length = np.linalg.norm(spoke)
        else:
           spoke_length = 0
        if i<549:#549:
           InnerSpokeLength[i,0] = spoke_length
        elif i>=549 and i<1098:
           InnerSpokeLength[i-549,1] = spoke_length#InnerSpokeLength[i-549,1] = spoke_length
        #print(InnerSpokeLength)
        #print(pt_array)
    for i in range(549):
        finished_points += 1
        lengthDiff = InnerSpokeLength[i,0]-InnerSpokeLength[i,1]
        #print(lengthDiff)
        if abs(lengthDiff)>eps_e:
           if InnerSpokeLength[i,0]>InnerSpokeLength[i,1] :
              pt = pt_array[i,:] 
              ps = ps_array[i,:]
              if np.linalg.norm(pt_array[i+549,:]-ps_array[i+549,:])!=0:
                 spoke = pt - ps
                 spoke_dir = spoke / InnerSpokeLength[i,1]
                 pt = ps + InnerSpokeLength[i,1]*spoke_dir
              else:
                 pt = ps
              InnerSpokeLength[i,0] = np.linalg.norm(pt-ps)
              pt_array[i, :] = pt
           elif InnerSpokeLength[i,0]<InnerSpokeLength[i,1]:
              pt = pt_array[i+549,:] #[i+549,:] 
              ps = ps_array[i+549,:] #[i+549,:]
              if np.linalg.norm(pt_array[i-549,:]-ps_array[i-549,:])!=0:
                 spoke = pt - ps
                 spoke_dir = spoke / InnerSpokeLength[i,1]
                 pt = ps + InnerSpokeLength[i,0]*spoke_dir
              else:
                 pt = ps
              InnerSpokeLength[i,1] = np.linalg.norm(pt-ps)
              pt_array[i+549, :] = pt#pt_array[i+549, :] = pt       
        #print('Equalize Spoke Length Finished Points: %s of Total 549' % str(finished_points))
    #print(InnerSpokeLength)
    #print(pt_array)
    pts = vtk.vtkPoints()
    for i in range(num_pt):
        pts.InsertNextPoint((pt_array[i, 0]),(pt_array[i, 1]),(pt_array[i, 2]))
    pt_vtk_EualizedLength = vtk.vtkPolyData()
    pt_vtk_EualizedLength.SetPoints(pts)      
    return pt_vtk_EualizedLength 

def ClosestSurfPoint(ps,surf_vtk):
    cell_locator = vtk.vtkCellLocator()
    cell_locator.SetDataSet(surf_vtk)
    cell_locator.BuildLocator()
    cellId = vtk.reference(0)
    c = [0.0, 0.0, 0.0]
    subId = vtk.reference(0)
    d = vtk.reference(0.0)
    cell_locator.FindClosestPoint(ps, c, cellId, subId, d)
    pt_ids = vtk.vtkIdList()
    surf_vtk.GetCellPoints(cellId, pt_ids)
    num_cell = pt_ids.GetNumberOfIds()
    #print(num_cell)
    vector_array = np.zeros((num_cell, 3))
    num_cell = 1
    p_closestPoints = np.zeros((num_cell, 3))
    for i in range(num_cell):
        pt_id = pt_ids.GetId(i)
        #Calculate normal vector of each point
        p_closestPoints[i,:] = np.array(surf_vtk.GetPoint(pt_id))
    #print(p_closestPoints)  
    p_closestPoint = p_closestPoints[0,:]
    #print(p_closestPoint)
    return p_closestPoint
 
def RepairSkeleton(surf_vtk,pt_vtk,ps_vtk):
    num_pt = pt_vtk.GetNumberOfPoints()
    pt_array = np.zeros((num_pt, 3))    
    num_ps = ps_vtk.GetNumberOfPoints()
    ps_array = np.zeros((num_ps, 3))
    #p = pv.Plotter()
    #p.add_mesh(surf_vtk, color='orange', opacity=0.4)
    #num_pt = 100 
    for i in range(num_ps):
        pt_array[i,:] = np.array(pt_vtk.GetPoint(i))
        ps_array[i,:] = np.array(ps_vtk.GetPoint(i))
        pt = pt_array[i,:] 
        ps = ps_array[i,:]
        IsInside = IsInsideCheck(ps[0],ps[1],ps[2],surf_vtk)
        if IsInside==0:
              #print('Point %s is inside surface'% str(i))
              #替换该点为最临近的曲面点
              ps = ClosestSurfPoint(ps,surf_vtk)
              pt = ps
        pt_array[i, :] = pt
        ps_array[i, :] = ps
    pts = vtk.vtkPoints()
    pss = vtk.vtkPoints()
    for j in range(num_pt):
        pts.InsertNextPoint((pt_array[j, 0]),(pt_array[j, 1]),(pt_array[j, 2]))
        pss.InsertNextPoint((ps_array[j, 0]),(ps_array[j, 1]),(ps_array[j, 2]))
    RepairSkeleton = vtk.vtkPolyData()
    RepairTips = vtk.vtkPolyData()
    RepairSkeleton.SetPoints(pss)  
    RepairTips.SetPoints(pts) 
    return RepairSkeleton, RepairTips
 
if __name__ == '__main__':

   parent_folder = '/home/hitlab/gaona/interp/revision/srep/R_surf_CN/tp0'  
   parent_surf_folder = '/home/hitlab/gaona/interp/revision/surf/R_surf_CN/tp0'

   eps_s = 0.1
   eps_d = 0.1
   eps_e = 0.1
   r_num = 1

   finished_cases = 0

   for case_folder in os.listdir(parent_folder):
       case_dir = os.path.join(parent_folder, case_folder)

       output_dir = case_dir + "/Repaired_pt.vtk" 
       output_dir2 = case_dir + "/Repaired_ps.vtk" 

       #把regression的结果分成两部分：tips，skeleton
       #filedir = case_dir + "/output/GeodesicRegression__GeodesicFlow__hippo__tp_3__age_3.00.vtk"
       #pt_vtk, ps_vtk = ArrangeRegressionResult(filedir)

       pt = case_dir + "/" + "pt.vtk"
       pt_reader = vtk.vtkPolyDataReader()
       pt_reader.SetFileName(pt)
       pt_reader.Update()
       pt_vtk = pt_reader.GetOutput() 

       ps = case_dir + "/" + "ps.vtk"
       ps_reader = vtk.vtkPolyDataReader()
       ps_reader.SetFileName(ps)
       ps_reader.Update()
       ps_vtk = ps_reader.GetOutput() 

       surf = parent_surf_folder + "/" + case_folder + ".vtk"
       surf_reader = vtk.vtkPolyDataReader()
       surf_reader.SetFileName(surf)
       surf_reader.Update()
       surf_vtk = surf_reader.GetOutput()      

       print('Start refine: %s ...'% str(case_folder))
       pt_vtk_RefinedLength = RefineSpokeLength(surf_vtk,pt_vtk,ps_vtk,eps_s)
       print('Finished RefinedLength1')

       writer = vtk.vtkPolyDataWriter()
       writer.SetFileName(output_dir)
       writer.SetInputData(pt_vtk_RefinedLength)
       #writer.SetDataModeToAscii()
       writer.Update()

       writer = vtk.vtkPolyDataWriter()
       writer.SetFileName(output_dir2)
       writer.SetInputData(ps_vtk)
       #writer.SetDataModeToAscii()
       writer.Update()
   
       finished_cases += 1
       print('Finished cases: %s'% str(finished_cases))
'''
   #Display results
   num_pt = pt_vtk.GetNumberOfPoints()
   pt_array = np.zeros((num_pt, 3))    
   num_ps = ps_vtk.GetNumberOfPoints()
   ps_array = np.zeros((num_ps, 3))
   p = pv.Plotter()
   p.add_mesh(surf_vtk, color='orange', opacity=0.4)
   for i in range(num_pt):
       pt = np.array(pt_vtk_RefinedLength.GetPoint(i))#Replace pt file to show before and after refine filter
       ps = np.array(ps_vtk.GetPoint(i))
       spoke = pt - ps
       spoke_length = np.linalg.norm(spoke)
       p.add_mesh(pv.Arrow(ps, spoke, scale=spoke_length), color='red')
   p.show() 
'''

