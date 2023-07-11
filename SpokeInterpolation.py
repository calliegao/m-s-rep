import vtk
import xml.etree.ElementTree as ET
import os
import pyvista as pv
import numpy as np
class SrepWrapper(object):
    def __init__(self, polydata, num_rows=5, num_cols=5):
        self.data = polydata
        self.rows = num_rows
        self.cols = num_cols
class Deformer:
    def __init__(self):
        pass
    def organize_srep(self, spoke_polydata, num_rows, num_cols):
        spokePointData = spoke_polydata.GetPointData()
        numberOfArrays = spokePointData.GetNumberOfArrays()
        if numberOfArrays is 0:
            print('empty array')

        spokePoints = vtk.vtkPoints()
        spokeLines = vtk.vtkCellArray()

        arr_length = spokePointData.GetArray('spokeLength')
        arr_dirs = spokePointData.GetArray('spokeDirection')
        base_pts_array = vtk.vtkDoubleArray()
        base_pts_array.SetNumberOfComponents(3)
        base_pts_array.SetName("basePoints")
        for i in range(spoke_polydata.GetNumberOfPoints()):
            pt = [0] * 3
            spoke_polydata.GetPoint(i, pt)
            # base point of up arrows
            id0 = spokePoints.InsertNextPoint(pt)
            base_pts_array.InsertNextTuple(pt)

            # head of up arrows
            spoke_length = arr_length.GetValue(i)
            baseIdx = i * 3
            dirX = arr_dirs.GetValue(baseIdx)
            dirY = arr_dirs.GetValue(baseIdx + 1)
            dirZ = arr_dirs.GetValue(baseIdx + 2)
            pt1 = [0] * 3
            pt1[0] = pt[0] + spoke_length * dirX
            pt1[1] = pt[1] + spoke_length * dirY
            pt1[2] = pt[2] + spoke_length * dirZ
            id1 = spokePoints.InsertNextPoint(pt1)

            up_arrow = vtk.vtkLine()
            up_arrow.GetPointIds().SetId(0, id0)
            up_arrow.GetPointIds().SetId(1, id1)
            spokeLines.InsertNextCell(up_arrow)

        renderable_srep = vtk.vtkPolyData()
        renderable_srep.SetPoints(spokePoints)
        renderable_srep.SetLines(spokeLines)
        renderable_srep.GetPointData().AddArray(arr_length)

        renderable_srep.GetPointData().AddArray(arr_dirs)
        renderable_srep.GetPointData().AddArray(base_pts_array)
        # organized_srep = {"skeletal_points_vtk": spokePoints,
        #                   "radii_da": arr_length,
        #                   "dirs_da": arr_dirs,
        #                   "num_rows": num_rows,
        #                   "num_cols": num_cols}
        return SrepWrapper(renderable_srep, num_rows, num_cols)
    def readSrepFromXML(self, filename):
        """ Parse header.xml file, create models from the data, and visualize it. """
        # 1. parse header file
        tree = ET.parse(filename)
        upFileName = ''
        crestFileName = ''
        downFileName = ''
        nCols = 0
        nRows = 0
        headerFolder = os.path.dirname(filename)
        for child in tree.getroot():
            if child.tag == 'upSpoke':
                upFileName = os.path.join(headerFolder, child.text.split('/')[-1])
            elif child.tag == 'downSpoke':
                downFileName = os.path.join(headerFolder, child.text.split('/')[-1])
            elif child.tag == 'crestSpoke':
                crestFileName = os.path.join(headerFolder, child.text.split('/')[-1])
            elif child.tag == 'nRows':
                nRows = (int)(child.text)
            elif child.tag == 'nCols':
                nCols = (int)(child.text)

        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(upFileName)
        reader.Update()

        upSpokes = reader.GetOutput()

        down_reader = vtk.vtkXMLPolyDataReader()
        down_reader.SetFileName(downFileName)
        down_reader.Update()
        downSpokes = down_reader.GetOutput()

        crest_reader = vtk.vtkXMLPolyDataReader()
        crest_reader.SetFileName(crestFileName)
        crest_reader.Update()
        crestSpokes = crest_reader.GetOutput()

        up_renderable    = self.organize_srep(upSpokes, nRows, nCols)
        down_renderable  = self.organize_srep(downSpokes, nRows, nCols)
        crest_renderable = self.organize_srep(crestSpokes, nRows, nCols)
        return up_renderable, down_renderable, crest_renderable, nRows, nCols
    def load_srep(self, srep_header_file):
        up_renderable, down_renderable, crest_renderable, num_rows, num_cols = \
                                                                                        self.readSrepFromXML(srep_header_file)
        up_poly, down_poly, crest_poly = up_renderable.data, down_renderable.data, crest_renderable.data
        
        return up_poly, down_poly, crest_poly
    def deform_to(self, surf0_file, srep0, surf1_file):  
        """The surface is deformed from surf0 to surf1
        Meanwhile, the s-rep of surf0 is deformed, yielding srep1
        """
        surf0_reader = vtk.vtkPolyDataReader()
        surf0_reader.SetFileName(surf0_file)
        surf0_reader.Update()
        surf0 = surf0_reader.GetOutput()

        surf1_reader = vtk.vtkPolyDataReader()
        surf1_reader.SetFileName(surf1_file)
        surf1_reader.Update()
        surf1 = surf1_reader.GetOutput()
        
        source_pts = vtk.vtkPoints()
        target_pts = vtk.vtkPoints()
        for i in range(surf0.GetNumberOfPoints()):
            source_pts.InsertNextPoint(surf0.GetPoint(i))
            target_pts.InsertNextPoint(surf1.GetPoint(i))
            # p = pv.Plotter()
            
            # p.add_mesh(surf1, color='orange', opacity=0.4)
            # p.add_mesh(surf0, color='white', opacity=0.4)
            # p.add_mesh(np.array(surf0.GetPoint(i)), point_size=20, render_points_as_spheres=True, color='red')
            # p.add_mesh(np.array(surf1.GetPoint(i)), point_size=20, render_points_as_spheres=True, color='blue')
            # p.show()
        tps = vtk.vtkThinPlateSplineTransform()
        tps.SetSourceLandmarks(source_pts)
        tps.SetTargetLandmarks(target_pts)
        tps.SetBasisToR()
        tps.Update()
        
        deformed_pts = []
        
        up_spokes, down_spokes, crest_spokes = srep0
        ## deform up spokes
        up_spokes_poly = vtk.vtkPolyData()
        up_spokes_pts = vtk.vtkPoints()
        up_spokes_vec = vtk.vtkCellArray()

        for i in range(up_spokes.GetNumberOfPoints()//2):
            up_skel_pt = np.array(tps.TransformPoint(up_spokes.GetPoint(i * 2)))
            up_bdry_pt = np.array(tps.TransformPoint(up_spokes.GetPoint(i * 2+1)))
            deformed_pts.append(up_skel_pt)

            id0 = up_spokes_pts.InsertNextPoint(up_skel_pt)
            id1 = up_spokes_pts.InsertNextPoint(up_bdry_pt)
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, id0)
            line.GetPointIds().SetId(1, id1)
            up_spokes_vec.InsertNextCell(line)
            # print(i)
            # p = pv.Plotter()
            
            # p.add_mesh(surf1, color='orange', opacity=0.4)
            # p.add_mesh(surf0, color='white', opacity=0.4)
            # p.add_mesh(pv.Arrow(up_skel_pt, up_bdry_pt-up_skel_pt, scale=3), color='red')
            # p.add_mesh(pv.Arrow(up_spokes.GetPoint(i * 2), np.array(up_spokes.GetPoint(i * 2+1))-np.array(up_spokes.GetPoint(i * 2)), scale=3), color='blue')
            # p.show()
        up_spokes_poly.SetPoints(up_spokes_pts)
        up_spokes_poly.SetLines(up_spokes_vec)

        ## deform up spokes
        down_spokes_poly = vtk.vtkPolyData()
        down_spokes_pts = vtk.vtkPoints()
        down_spokes_vec = vtk.vtkCellArray()

        for i in range(down_spokes.GetNumberOfPoints()//2):
            down_skel_pt = deformed_pts[i]
            down_bdry_pt = np.array(tps.TransformPoint(down_spokes.GetPoint(i * 2+1)))
            id0 = down_spokes_pts.InsertNextPoint(down_skel_pt)
            id1 = down_spokes_pts.InsertNextPoint(down_bdry_pt)
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, id0)
            line.GetPointIds().SetId(1, id1)
            down_spokes_vec.InsertNextCell(line)
        down_spokes_poly.SetPoints(down_spokes_pts)
        down_spokes_poly.SetLines(down_spokes_vec)

        ## deform crest spokes
        crest_spokes_poly = vtk.vtkPolyData()
        crest_spokes_pts = vtk.vtkPoints()
        crest_spokes_vec = vtk.vtkCellArray()

        for i in range(crest_spokes.GetNumberOfPoints()//2):
            crest_skel_pt = np.array(tps.TransformPoint(crest_spokes.GetPoint(i * 2)))
            crest_bdry_pt = np.array(tps.TransformPoint(crest_spokes.GetPoint(i * 2+1)))

            id0 = crest_spokes_pts.InsertNextPoint(crest_skel_pt)
            id1 = crest_spokes_pts.InsertNextPoint(crest_bdry_pt)
            line = vtk.vtkLine()
            line.GetPointIds().SetId(0, id0)
            line.GetPointIds().SetId(1, id1)
            crest_spokes_vec.InsertNextCell(line)
        crest_spokes_poly.SetPoints(crest_spokes_pts)
        crest_spokes_poly.SetLines(crest_spokes_vec)
        return up_spokes_poly, down_spokes_poly, crest_spokes_poly, surf1
    def write_to_vtp(self, file_name, vtk_poly):
        new_poly = vtk.vtkPolyData()
        new_poly_pts = vtk.vtkPoints()
        new_poly_vec = vtk.vtkCellArray()
        radii_array = vtk.vtkDoubleArray()
        radii_array.SetNumberOfComponents(1)
        radii_array.SetName("spokeLength")

        direction_array = vtk.vtkDoubleArray()
        direction_array.SetNumberOfComponents(3)
        direction_array.SetName("spokeDirection")
        for i in range(vtk_poly.GetNumberOfPoints()//2):
            skel_pt = np.array(vtk_poly.GetPoint(i*2))
            bdry_pt = np.array(vtk_poly.GetPoint(i*2+1))
            spoke = bdry_pt - skel_pt
            spoke_len = np.linalg.norm(spoke)
            spoke_dir = spoke / spoke_len

            new_poly_pts.InsertNextPoint(skel_pt)
            direction_array.InsertNextTuple(spoke_dir)
            radii_array.InsertNextValue(spoke_len)
        new_poly.SetPoints(new_poly_pts)
        new_poly.GetPointData().AddArray(radii_array)
        new_poly.GetPointData().AddArray(direction_array)
        new_poly.Modified()
        
        poly_writer = vtk.vtkXMLPolyDataWriter()
        poly_writer.SetFileName(file_name)
        poly_writer.SetInputData(new_poly)
        poly_writer.Update()
            

if __name__ == '__main__':
    
   deformer = Deformer()            
   parent_folder = '/home/hitlab/gaona/interp/revision/srep/R_surf_CN/tp2'   
   parent_surf_folder = '/home/hitlab/gaona/interp/revision/surf/R_surf_CN/tp2'
   finished_cases = 0

   for case_folder in os.listdir(parent_folder):
       case_dir = os.path.join(parent_folder, case_folder)

       target_srep_xml = case_dir + '/header.xml'
       up_spokes, down_spokes, crest_spokes = deformer.load_srep(target_srep_xml)
       t0_surf = parent_surf_folder + '/' + case_folder + '.vtk'
       output_dir = case_dir
       file_name = 'interplated.vtk'      

       ##########################
       ###### Interpolate S-reps
       ##########################
       from interpolater import Interpolater
       surf_reader = vtk.vtkPolyDataReader()
       surf_reader.SetFileName(t0_surf)
       surf_reader.Update()
       surf_vtk = surf_reader.GetOutput()
       ## the interp_level = 1 will produce 486 spokes, while interp_level=2 will produce 1700+ spokes.
       interp = Interpolater(2)
       primary_spokes_app = vtk.vtkAppendPolyData()
       primary_spokes_app.AddInputData(up_spokes)
       primary_spokes_app.AddInputData(down_spokes)
       primary_spokes_app.AddInputData(crest_spokes)
       primary_spokes_app.Update()
       prim_spokes = primary_spokes_app.GetOutput()
    
       interp_spokes = interp.interpolate_all(prim_spokes, 16, 3)

       spoke_vtk = vtk.vtkPolyData()
       spoke_end = vtk.vtkPoints()
       spoke_vec = vtk.vtkCellArray()

       for s in interp_spokes:
           id0 = spoke_end.InsertNextPoint(s.p)
           id1 = spoke_end.InsertNextPoint(s.getB())
           spoke_line = vtk.vtkLine()
           spoke_line.GetPointIds().SetId(0, id0)
           spoke_line.GetPointIds().SetId(1, id1)
           spoke_vec.InsertNextCell(spoke_line)
       spoke_vtk.SetPoints(spoke_end)
       spoke_vtk.SetLines(spoke_vec)
    
       writer = vtk.vtkPolyDataWriter()
       writer.SetFileName(os.path.join(output_dir, file_name))
       writer.SetInputData(spoke_vtk)
       #writer.SetDataModeToAscii()
       writer.Update()

       finished_cases = finished_cases + 1
       print('Finished cases: %s'% str(finished_cases))
   print('done')


   p = pv.Plotter()
   p.add_mesh(surf_vtk, color='white', opacity=0.3)
   p.add_mesh(spoke_vtk, line_width=4, color='red')
   p.add_mesh(prim_spokes, line_width=4, color='green')
   p.show()
    
    ###########################
    ##### Generate longitudinal s-reps
    ##################
    # up_spokes_poly, down_spokes_poly, crest_spokes_poly, surf1 = deformer.deform_to(t0_surf, [up_spokes, down_spokes, crest_spokes], t1_surf)

    # deformer.write_to_vtp('/Users/zhiy/PhD/my_paper/Gao/test2/t1_up.vtp', up_spokes_poly)
    # deformer.write_to_vtp('/Users/zhiy/PhD/my_paper/Gao/test2/t1_down.vtp', down_spokes_poly)
    # deformer.write_to_vtp('/Users/zhiy/PhD/my_paper/Gao/test2/t1_crest.vtp', crest_spokes_poly)

    # # up_spokes_2, down_spokes_2, crest_spokes_2 = deformer.load_srep('/Users/zhiy/PhD/my_paper/Gao/test_data/t1_srep/header.xml')
    # # compare_appender = vtk.vtkAppendPolyData()
    # # compare_appender.AddInputData(up_spokes_2)
    # # compare_appender.AddInputData(down_spokes_2)
    # # compare_appender.AddInputData(crest_spokes_2)
    # # compare_appender.Update()

    # spharm_pts = []
    # for i in range(surf1.GetNumberOfPoints()):
    #     spharm_pts.append(np.array(surf1.GetPoint(i)))
    
    # p = pv.Plotter()
    # p.add_mesh(up_spokes_poly, color='cyan',line_width=4)
    # p.add_mesh(surf1, color='white', opacity=0.3)
    # p.add_mesh(down_spokes_poly, color='magenta', line_width=4)
    # p.add_mesh(crest_spokes_poly, color='yellow',line_width=4)
    # p.add_mesh(np.array(spharm_pts), color='white', point_size=10, render_points_as_spheres=True)

    # # p.add_mesh(compare_appender.GetOutput(), color='blue', line_width=4)
    # p.show()

