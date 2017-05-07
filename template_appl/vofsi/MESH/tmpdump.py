# -*- coding: utf-8 -*-

###
### This file is generated automatically by SALOME v7.7.1 with dump python functionality
###

import sys
import salome

salome.salome_init()
theStudy = salome.myStudy

import salome_notebook
notebook = salome_notebook.NoteBook(theStudy)
sys.path.insert( 0, r'/home/daniele/software/femus/USER_APPL/vofsi_1/MESH')

###
### GEOM component
###

import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS


geompy = geomBuilder.New(theStudy)

O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)
geomObj_1 = geompy.MakeVertex(0, 0, 0)
geomObj_2 = geompy.MakeVectorDXDYDZ(1, 0, 0)
geomObj_3 = geompy.MakeVectorDXDYDZ(0, 1, 0)
geomObj_4 = geompy.MakeVectorDXDYDZ(0, 0, 1)
geomObj_5 = geompy.MakeVertex(0, 0, 0)
geomObj_6 = geompy.MakeVertex(1, 0, 0)
geomObj_7 = geompy.MakeVertex(1, 1, 0)
geomObj_8 = geompy.MakeVertex(0, 1, 0)
O_1 = geompy.MakeVertex(0, 0, 0)
OX_1 = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY_1 = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ_1 = geompy.MakeVectorDXDYDZ(0, 0, 1)
Vertex_1 = geompy.MakeVertex(0, 0, 0)
Vertex_2 = geompy.MakeVertex(1, 0, 0)
Vertex_3 = geompy.MakeVertex(1, 1, 0)
Vertex_4 = geompy.MakeVertex(0, 1, 0)
Vertex_5 = geompy.MakeVertex(0.4, 0, 0)
Vertex_6 = geompy.MakeVertex(0.6, 0, 0)
Vertex_7 = geompy.MakeVertex(0.4, 0.6, 0)
Vertex_8 = geompy.MakeVertex(0.6, 0.6, 0)
Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_5)
Line_2 = geompy.MakeLineTwoPnt(Vertex_5, Vertex_7)
Line_3 = geompy.MakeLineTwoPnt(Vertex_7, Vertex_8)
Line_4 = geompy.MakeLineTwoPnt(Vertex_8, Vertex_6)
Line_5 = geompy.MakeLineTwoPnt(Vertex_6, Vertex_2)
Line_6 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_3)
Line_7 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_4)
Line_1_vertex_2 = geompy.GetSubShape(Line_1, [2])
Line_8 = geompy.MakeLineTwoPnt(Vertex_4, Line_1_vertex_2)
Line_1_vertex_3 = geompy.GetSubShape(Line_1, [3])
Line_4_vertex_3 = geompy.GetSubShape(Line_4, [3])
Line_9 = geompy.MakeLineTwoPnt(Line_1_vertex_3, Line_4_vertex_3)
Line_3_vertex_3 = geompy.GetSubShape(Line_3, [3])
Line_7_vertex_2 = geompy.GetSubShape(Line_7, [2])
Line_10 = geompy.MakeLineTwoPnt(Line_3_vertex_3, Line_7_vertex_2)
Line_8_vertex_2 = geompy.GetSubShape(Line_8, [2])
Line_2_vertex_3 = geompy.GetSubShape(Line_2, [3])
Line_11 = geompy.MakeLineTwoPnt(Line_8_vertex_2, Line_2_vertex_3)
Face_1 = geompy.MakeFaceWires([Line_1, Line_2, Line_8, Line_11], 1)
geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )
geompy.addToStudy( O_1, 'O' )
geompy.addToStudy( OX_1, 'OX' )
geompy.addToStudy( OY_1, 'OY' )
geompy.addToStudy( OZ_1, 'OZ' )
geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Vertex_5, 'Vertex_5' )
geompy.addToStudy( Vertex_6, 'Vertex_6' )
geompy.addToStudy( Vertex_7, 'Vertex_7' )
geompy.addToStudy( Vertex_8, 'Vertex_8' )
geompy.addToStudy( Line_1, 'Line_1' )
geompy.addToStudy( Line_2, 'Line_2' )
geompy.addToStudy( Line_3, 'Line_3' )
geompy.addToStudy( Line_4, 'Line_4' )
geompy.addToStudy( Line_5, 'Line_5' )
geompy.addToStudy( Line_6, 'Line_6' )
geompy.addToStudy( Line_7, 'Line_7' )
geompy.addToStudyInFather( Line_1, Line_1_vertex_2, 'Line_1:vertex_2' )
geompy.addToStudy( Line_8, 'Line_8' )
geompy.addToStudyInFather( Line_1, Line_1_vertex_3, 'Line_1:vertex_3' )
geompy.addToStudyInFather( Line_4, Line_4_vertex_3, 'Line_4:vertex_3' )
geompy.addToStudy( Line_9, 'Line_9' )
geompy.addToStudyInFather( Line_3, Line_3_vertex_3, 'Line_3:vertex_3' )
geompy.addToStudyInFather( Line_7, Line_7_vertex_2, 'Line_7:vertex_2' )
geompy.addToStudy( Line_10, 'Line_10' )
geompy.addToStudyInFather( Line_8, Line_8_vertex_2, 'Line_8:vertex_2' )
geompy.addToStudyInFather( Line_2, Line_2_vertex_3, 'Line_2:vertex_3' )
geompy.addToStudy( Line_11, 'Line_11' )
geompy.addToStudy( Face_1, 'Face_1' )


if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(1)
