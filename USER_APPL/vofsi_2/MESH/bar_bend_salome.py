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

Vertex_1 = geompy.MakeVertex(0, 0.0, 0)
Vertex_2= geompy.MakeVertex(1, 0, 0)
Vertex_3= geompy.MakeVertex(1, 1.0, 0)
Vertex_4 = geompy.MakeVertex(0, 1.0, 0)
Vertex_5 = geompy.MakeVertex(0.4, 0.0 ,0)
Vertex_6 = geompy.MakeVertex(0.6, 0.0, 0)
Vertex_7 = geompy.MakeVertex(0.4, 0.6, 0)
Vertex_8 = geompy.MakeVertex(0.6, 0.6, 0)

Line_1 = geompy.MakeLineTwoPnt(Vertex_1, Vertex_5)
Line_2 = geompy.MakeLineTwoPnt(Vertex_5, Vertex_7)
Line_3 = geompy.MakeLineTwoPnt(Vertex_7, Vertex_8)
Line_4 = geompy.MakeLineTwoPnt(Vertex_8, Vertex_6)
Line_5 = geompy.MakeLineTwoPnt(Vertex_6, Vertex_2)
Line_6 = geompy.MakeLineTwoPnt(Vertex_2, Vertex_3)
Line_7 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_4)
Line_8 = geompy.MakeLineTwoPnt(Vertex_4, Vertex_1)
Line_9 = geompy.MakeLineTwoPnt(Vertex_5, Vertex_6)
Line_10 = geompy.MakeLineTwoPnt(Vertex_3, Vertex_8)
Line_11 = geompy.MakeLineTwoPnt(Vertex_4, Vertex_7)






geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )



geompy.addToStudy( Vertex_1, 'Vertex_1' )
geompy.addToStudy( Vertex_2, 'Vertex_2' )
geompy.addToStudy( Vertex_3, 'Vertex_3' )
geompy.addToStudy( Vertex_4, 'Vertex_4' )
geompy.addToStudy( Vertex_5, 'Vertex_5' )
geompy.addToStudy( Vertex_6, 'Vertex_6' )
geompy.addToStudy( Vertex_7, 'Vertex_7' )
geompy.addToStudy( Vertex_8, 'Vertex_8' )






if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser(1)
