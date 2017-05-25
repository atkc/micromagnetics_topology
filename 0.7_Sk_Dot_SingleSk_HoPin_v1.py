#everything is in IS
# say not to outdated not standard units of measurementent such as CGS
#eV is cool though
import os
import sys
import re
import shutil 
import math 
import numpy as np
from itertools import *
from textwrap import *
import glob 
import tarfile
import time
import pandas as pd
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)

import queue
import threading

import plotly.plotly as py
from plotly.graph_objs import *
py.sign_in('Gattsus', 'zb8uk03bif')
#py.sign_in('DemoAccount', 'lr1c37zw81')
#import matplotlib.pyplot as plt
#import cufflinks as cf

#import pdb
#pdb.set_trace()

initial=os.getcwd()
#mumax_command='/home/anibal/go/bin/mumax3'
mumax_command='/Users/Anthony/OneDrive/PHYSICS/FYP-skyrmions/Micromagnetics/mumax3'
def iscomment(s):
   return s.startswith('#')

def vector_normalization(line, i0):# what is this for?

   vector=[]
   i_min=i0-1
  
   for i in range(i_min, (i_min+3), 1):  # columns start with 0
      vector.append(line.split()[i])

   vector=array(vector, dtype=float)

   Norm_vector=linalg.norm(vector)
   if (Norm_vector <= 0.1):
      vector=[0, 0, 0]
   else:
      vector=[x/Norm_vector for x in vector]

   return vector

def drange(start, stop, step):# what is this for?
    r = start
    while (r <= (stop+step/2.0)):
        yield r
        r += step

def countdown(t): # timer to wait before plotting again - but for fuck?
    for i in range(t,0,-1):
        print('waiting %d seconds, before trying to plot again \r' % i)
        sys.stdout.flush()
        time.sleep(1)
#start_plot_hyst
def Plot_Hyst(fig_name, x_list, y_list, plot_title, curves_names, x_Axis_Title, y_Axis_Title, x_range, y_range, plot_symbols, plot_modes, AxisFontSize, PlotFontSize, TicksFontSize): #, sim_mx3,  Ex, B_Max, M, D, K1, K2, size_X, size_Y, size_Z, Method):

   #yRange= 1.05*max( abs(i) for i in Y0_in )

   #yRange= 1.05*max(abs(y_list) ) # [func(l) for l in main_array for func in (min, max)] 
   #yRange = [func(l) for l in y_list for func in (min, max)]
   #yRange=1 
   
   print("Plotting hyst started ")

   traces=[]
   for i in np.arange(len(x_list)):
      traces.append(
         Scatter(
            x=x_list[i],
            y=y_list[i],
            mode=plot_modes[i],
            name=curves_names[i],
            marker=Marker(
               symbol=plot_symbols[i],
               size=12
            )
         )
      )

   data=Data(traces) 

   layout = Layout(

      autosize=False,
      width=1280,
      height=720,

      title=plot_title,

      titlefont=Font(
         family='Helvetica',
         size=PlotFontSize
      ),

      xaxis=XAxis(
         title=x_Axis_Title,
         titlefont=Font(
            family='Helvetica',
            size=AxisFontSize,
            color='black'
         ),
         range=x_range,
         showgrid=True,
         gridcolor='7B7676',
         zeroline=False,
         showline=True, 
         linewidth=2,
         mirror='ticks',
         showticklabels=True,
         tickfont=Font(
            family='Helvetica',
            size=TicksFontSize,
            color='black'
         ),
      ),
 
      yaxis=YAxis(
         title=y_Axis_Title,
         titlefont=Font(
            family='Helvetica',
            size=AxisFontSize,
            color='black'
         ),
         range=y_range, #hardcoded
         showgrid=True,
         gridcolor='7B7676',
         showline=True,
         linewidth=2,
         mirror='ticks',
         showticklabels=True,
         tickfont=Font(
            family='Helvetica',
            size=TicksFontSize,
            color='black'
         )
      ),

#      yaxis2=YAxis(
#         title='&tau;<sub>z</sub> [T]',
#         titlefont=Font(
#            family='Helvetica',
#            size=AxisFontSize,
#            color='black'
#         ),
#         range=[-1.05,1.05], #hardcoded
#         showgrid=False,
##         gridcolor='7B7676',
#         showline=True,
#         linewidth=2,
##         mirror='ticks',
##         showticklabels=True,
#         tickfont=Font(
#            family='Helvetica',
#            size=TicksFontSize,
#            color='black'
#         ),
#         overlaying='y',
#         side='right' 
#      ),

      showlegend=True,
#      showlegend=False,
      legend=Legend(
         x=1.05,
         y=1,
         font=Font(
            family='Helvetica',
            size=AxisFontSize,
            color='black'
         )
      )
   )    #)layout

   fig = Figure(data=data, layout=layout)
   fig_name

   py.image.save_as(fig, fig_name)

#   plot_url = py.plot(fig, filename =fig_name, world_readable=False)
   print("Plot done, saved at "+fig_name)
#end_plot_hyst



def writting_mx3(sim_mx3, Damping, Ex, M, D, K1_eff, K1, K2, B_Min, B_Max, B_Step, size_X, size_Y, size_Z, Nx, Ny, Nz, nrepeats, Method, Index, aux, working_dir, results_dir_Plots, results_dir_PlotsM, results_dir_Data, src_dir, src_data, dst_data):


   #### 0.0 dir creation
                   
   #aux_m0="v"+sim_mx3+"m0_"+Configuration_Fl+"_"+str("%.3f"%BiasFieldX)+"_"+str("%.2f"%M_Fl)+"_"+str("%.2f"%K_Fl_eff_m0)+"/"+sim_mx3+"m0.out/"
   #m0_file=os.path.join(initial,aux_m0+"m0.ovf")

   if not os.path.exists(working_dir): 
      os.makedirs(working_dir)
   print(working_dir)
   
   if not os.path.exists(working_dir): 
      print(working_dir+" doesn't exist")

   if not os.path.exists(results_dir_Data): 
      os.makedirs(results_dir_Data)
   print(results_dir_Data)

   if not os.path.exists(results_dir_Plots): 
      os.makedirs(results_dir_Plots)
   print(results_dir_Plots)

   
   Nz_list = np.arange(Nz) # it doesn't iterate over the last number; Nz = 5, Nz_list = [0 1 2 3 4]
   """
   if Nz%2 == 0:
      index_middle = Nz_list[ int(Nz/2.0)] #we get the centre of the list  Nz = 5, index_middle = 2
      print(index_middle)
   else:
      index_middle = Nz_list[int(np.ceil(Nz/2.0))] #we get the centre of the list  Nz = 5, index_middle = 2
      print(index_middle)
   """
   index_middle = Nz_list[int(np.ceil(Nz/2.0))] #we get the centre of the list  Nz = 5, index_middle = 2
   
      
   layers_names_list = []

   for aux in np.arange(Nz) : #ant_why not aux in Nz_list???
      layers_names_list.append('layer'+str(aux))

   
   layers_names_list = []

   for aux in np.arange(Nz) : #ant_why do it another time?
      layers_names_list.append('layer'+str(aux))
   """   
   layers         = ''
   Damping_layers = ''
   Ex_layers     = ''
   M_layers      = ''
   DMI_layers    = ''
   Ku1_layers    = ''
   M0_layers     = ''
   def_regions   = ''
   """   

   #layers = layers_names_list[0]+' := cuboid(size_X*Nano, size_Y*Nano, size_Z*'+str(Nz)+'*Nano)'
   layers   = layers_names_list[0]+' := cylinder(size_X*Nano, size_Z*'+str(Nz)+'*Nano)'
   geometry       = 'setgeom('+layers_names_list[0]+')'

   """      
   for aux in Nz_list:
      #layers += layers_names_list[aux]+' :=cuboid(size_X*Nano, size_Y*Nano, size_Z*Nano).transl(0,0,'+str((aux-(nrepeats-1))*size_Z)+'*Nano) \n'

      layers += layers_names_list[aux]+' :=cuboid(size_X*Nano, size_Y*Nano, size_Z*Nano).transl(0,0,'+str((aux-index_middle)*size_Z)+'*Nano) \n   '

      #def_regions  +=  'defregion('+str(aux)+', '+layers_names_list[aux]+')\n'
      def_regions  +=  'defregion('+str(aux)+',layer( '+str(aux)+'))\n   '
      #def_regions  +=  'defregion('+str(aux)+')\n'

      if aux%2== 0:
         Damping_layers += 'alpha.setregion('+str(aux)+', Damping) \n   '
         Ex_layers    += 'Aex.setregion('+str(aux)+', Exchange) \n   '
         M_layers     += 'Msat.setregion('+str(aux)+', Mag) \n   '
         DMI_layers   += 'Dind.setregion('+str(aux)+', D) \n   '
         Ku1_layers   += 'Ku1.setregion('+str(aux)+', K1) \n   '
         M0_layers    += 'm.setregion('+str(aux)+', randomMagSeed(Seed)) \n   '      
      else:
         Damping_layers += 'alpha.setregion('+str(aux)+', 0) \n   '
         Ex_layers    += 'Aex.setregion('+str(aux)+', 0) \n   '
         M_layers     += 'Msat.setregion('+str(aux)+', 0) \n   '
         DMI_layers   += 'Dind.setregion('+str(aux)+', 0) \n   '
         Ku1_layers   += 'Ku1.setregion('+str(aux)+', 0) \n   '
         M0_layers    += 'm.setregion('+str(aux)+',  uniform(0.0, 0.0, 0.0)) \n   '      
   """
   """
   geometry       = 'setgeom('+layers_names_list[0]         
   for aux in np.arange(1,Nz):
      geometry +='.add('+layers_names_list[aux]+')'
   geometry+=') \n   '
   """
   M0_layers = dedent("""
    for i:=0; i<Nx; i++{
	for j:=0; j<Ny; j++{
		for k:=0; k<Nz; k++{
		   //r := index2coord(i, j, k)
		   //x := r.X()
		   //y := r.Y()
		   //z := r.Z()	

		   //theta := randNorm()*pi + 0.0 //Z-angle, mean of 2*pi
		   theta := rand()*pi     + 0.0  //XY-angle, mean of 0.0
		   phi   := rand()*2*pi     + 0.0  //XY-angle, mean of 0.0
		   ////if ( (sqrt((i-Nx/2.0)*(i-Nx/2.0) + (i-Ny/2.0)*(i-Ny/2.0)) <=  sqrt((Nx/2.0)*(Nx/2.0) + (Ny/2.0)*(Ny/2.0)) )/2.0 and Mod(k, 2) <= 0.1 ) {
		   //if  sqrt((i-Nx/2.0)*(i-Nx/2.0) + (j-Ny/2.0)*(j-Ny/2.0)) <=  Nx/2.0  {
		   //   if  Mod(k, 2) == 0.0 {
   		   //       m.setcell(i,j,k, vector(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)))
		   //    } else {
		   //       m.setcell(i,j,k, vector(0,0 ,0))	
		   //   }   
		   //} else {
		   //  m.setcell(i,j,k, vector(0,0 ,0))	
		   //}

		       if  Mod(k, 2) == 0 {
   		          m.setcell(i,j,k, vector(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)))
		       } else {
		     m.setcell(i,j,k, vector(0,0 ,0))	
		   }

		}
	   }
        }""")


   M0_layers     = ''

   for aux in Nz_list:
      if aux%4== 0:
         #M0_layers+= 'm.setinshape(cuboid(size_X*nano, size_X*nano, size_Z*nano).transl(0, 0, '+str(aux-index_middle)+'*size_Z*nano), NeelSkyrmion(1,-1).scale(3,3,1) )\n'
         #M0_layers+= 'm.setinshape(cylinder(size_X*nano, size_Z*nano).transl(0, 0, '+str(aux-index_middle+1.0)+'*size_Z*nano), NeelSkyrmion(1,-1).scale(1,1,1) )\n' 
         M0_layers+= 'm.setinshape(cylinder(size_X*nano, size_Z*nano).transl(0, 0, '+str(aux-index_middle+1.0)+'*size_Z*nano), RandomMag().scale(1,1,1) )\n' 
      else:
         #M0_layers+= 'm.setinshape(cuboid(size_X*nano, size_X*nano, size_Z*nano).transl(0, 0, '+str(aux-index_middle)+'*size_Z*nano), uniform(0.0, 0.0, 0.0))\n'
         M0_layers+= 'm.setinshape(cylinder(size_X*nano, size_Z*nano).transl(0, 0, '+str(aux-index_middle+1.0)+'*size_Z*nano), uniform(0.0, 0.0, 0.0))\n' 
   

   """
   print(layers)
   print(geometry)
   print(Damping_layers)
   print(def_regions)
   print(Ex_layers)
   print(M_layers)
   print(DMI_layers)
   print(Ku1_layers)
   print(M0_layers)
   """
   
   mumax_commands=dedent("""

   Pico :=1e-12
   Mega :=1e6
   Nano :=1e-9
   Mili :=1e-3

   Damping  :=%f       // 
   Exchange :=%f*Pico  // in J/m^3
   Mag      :=%f*Mega  // in A/m
   D        :=%f*Mili  // in J/m^2
   K1       :=%f*Mega  // in J/m^3
   Bmin     :=%f       // BZ in T
   Bmax     :=%f       // BZ in T
   Bstep    :=%f 

   size_X   :=%f
   size_Y   :=%f
   size_Z   :=%f

   Nx   :=%.0f
   Ny   :=%.0f
   Nz   :=%.0f

   Seed:=%.0f

   SetGridsize(Ny, Ny, Nz)
   SetCellsize(size_X*Nano/Nx, size_Y*Nano/Ny, size_Z*Nano)
   %s // layers
   %s // geometry

   alpha = Damping   
   Aex   = Exchange
   Msat  = Mag
   Dind  = D     //   Dbulk  = D
   Ku1   = K1

   //m = randomMagSeed(0)

   %s // M0_layers

   anisU  = vector(0, 0, 1)

   //B_ext = vector(0, 0, B_Max) //in mT - doh not A/m

   TableAdd(B_ext)
   TableAdd(E_Total)
   tableAdd(ext_topologicalcharge)
   OutputFormat = OVF1_TEXT


   tablesave()
   //unindented-ANT
   MinimizerStop = 1e-6
   %s        // high-energy states best minimized by relax()
   save(m)
   //unindented-ANT

   
   //virgin curve
   for B:=0.0; B<=Bmax; B+=Bstep{
       B_ext = vector(0, 0, B)
       minimize()   // small changes best minimized by minimize()
       tablesave()
       save(m)
   }

   """% (Damping, Ex, M, D, K1, B_Min, B_Max, B_Step, size_X, size_Y, size_Z, Nx, Ny, Nz, Index, layers, geometry, M0_layers, Method))
   #   """% (Damping, Ex, M, D, K1, B_Max, size_X, size_Y, size_Z, Nx, Ny, Nz, Index, def_regions, Damping_layers,  Ex_layers, M_layers, DMI_layers, Ku1_layers, M0_layers)) 

   """
   alpha = Damping   
   Aex   = Exchange
   Msat  = Mag
   Dind  = D     //   Dbulk  = D
   Ku1   = K1
   %s //def_regions
   %s // Damping_layers
   %s // Ex_layers
   %s // M_layers
   %s // D_layers
   %s // Ku1_layers
   %s // M0_layers


   %s // M0_layers

   i=0
   for B:=0.0; B>=-Bmax; B-=Bstep{
       B_ext = vector(0, 0, B)
       %s
       tablesave()
       save(m)
       i= i+1
   }
   """
   

   #print(mumax_commands)
   
   """
   //virgin curve
   for B:=0.0; B<=Bmax; B+=Bstep{
       B_ext = vector(0, 0, B)
       minimize()   // small changes best minimized by minimize()
       tablesave()
   }


   //autosave(m,Pico)
   //tableautosave(Pico)
   //run(20*Nano) 
   """

   #defining the location of the .mx3 script
   executable=os.path.join(working_dir, sim_mx3+".mx3")

   #opening and saving it
   executable_file=open(executable, "w")
   executable_file.write(mumax_commands)
   executable_file.close()

   return 0
   
def plotting_m(sim_mx3, Damping, Ex, M, D, K1_eff, K1, K2, B_Min, B_Max, B_Step, size_X, size_Y, size_Z, Nx, Ny, Nz, nrepeats, Method, Index, aux, working_dir, results_dir_Plots, results_dir_PlotsM, results_dir_Data, src_dir, src_data, dst_data):


   if not os.path.exists(results_dir_Plots): 
      os.makedirs(results_dir_Plots)
   print(results_dir_Plots)
   
   print('plotting_m STARTED')
   ###2.0 plotting the magnetisation
   
   #src_file = os.path.join(src_dir, 'm000000.png')
   #src_file_list = [os.path.join(src_dir, 'm000000.png'), os.path.join(src_dir, 'm000001.png')]
   
   src_file_list = glob.iglob(os.path.join(src_dir, "*.png"))
   if len(list(src_file_list)) >= 1:
      for src_file in src_file_list:
         if os.path.isfile(src_file):
            os.remove(src_file)

   #os.system("rm "+src_dir+"*.png")
   for item in list(glob.iglob(os.path.join(src_dir, "*.ovf"))):
      for i in np.arange(Nz)[::4]:
         m_plot_command = mumax_command+'-convert -comp=2 -color=\"white,gray,black\" -resize '+str(int(5*Nx))+'x'+str(int(5*Ny))+'x1 -zrange '+str(i)+':'+str(i+1)+' -png '+item
         print('plotting layer'+str(i))
         print(m_plot_command)
         os.system(m_plot_command)
         #os.system(mumax_command+'-convert -comp=2 -color=\"blue,white,red\" -resize 512x512x1 -zrange '+str(i)+':'+str(i+1)+' -png '+item)
         #os.system(mumax_command+'-convert -comp=2 -color=\"black,white,white\" -resize 512x512x1 -zrange '+str(i)+':'+str(i+1)+' -png '+src_dir+'*.ovf')
         #os.system(mumax_command+'-convert -comp=2 -color=\"white,white,black\" -resize 512x512x1 -zrange '+str(i)+':'+str(i+1)+' -png '+src_dir+'*.ovf')

         src_file_list = glob.iglob(os.path.join(src_dir, "*.png"))
         #print(list(src_file_list))
         for src_file in src_file_list:
            dst_file = os.path.join(results_dir_Plots, aux.replace('/','_')[0:len(aux.replace('/','_'))-1]+'_nlayer'+str(i)+'_'+os.path.splitext(os.path.basename(src_file))[0]+'.png') #os.path.splitext("path_to_file")[0]

            if os.path.isfile(src_file):
               print(src_file)
               if os.path.exists(dst_file):
                  os.remove(dst_file)
               shutil.move(src_file, dst_file)
               print("M plot copied to "+dst_file)
            #os.remove(src_file)

   """
   #os.system("rm "+src_dir+"*.png")
   print('plotting layer'+str(i))
   print(mumax_command+'-convert -comp=2 -color=\"blue,white,red\" -resize 512x512x1 -zrange '+str(i)+':'+str(i+1)+' -png '+src_dir+'*.ovf')
   os.system(mumax_command+'-convert -comp=2 -color=\"blue,white,red\" -resize 512x512x1 -zrange '+str(size_Z*i)+':'+str(size_Z*(i+1))+' -png '+src_dir+'*.ovf')

   src_file_list = glob.iglob(os.path.join(src_dir, "*.png"))
   #print(list(src_file_list))
   for src_file in src_file_list:
      dst_file = os.path.join(results_dir_Plots, aux.replace('/','_')[0:len(aux.replace('/','_'))-1]+'_nlayer'+str(i)+'_'+os.path.splitext(os.path.basename(src_file))[0]+'.png') #os.path.splitext("path_to_file")[0]

   if os.path.isfile(src_file):
      if os.path.exists(dst_file):
         os.remove(dst_file)
      shutil.move(src_file, dst_file)
   print("M plot copied to "+dst_file)
   """         

   """
   ###3.1 move the magnetisation plots to the plot dir

   files = glob.iglob(os.path.join(src_dir, "*.png"))
   dst_dir=results_dir_PlotsM

   for src_file in files:
      dst_file=os.path.join(dst_dir, os.path.basename(src_file))
      if os.path.isfile(src_file):
         if os.path.exists(dst_file):
            os.remove(os.path.join(dst_dir,os.path.basename(dst_file)))
         shutil.move(src_file, dst_dir)
   print("Magnetisation plots copied to "+dst_dir)

   ###3.2 Compressing *.ovf files
   src_files = glob.glob(working_dir+'0.4_SK.out/*.ovf')
   dst_tarfile=working_dir+sim_mx3+'.out/ovf_files.tgz'

   print('compressing OVF files STARTED '+dst_tarfile)
   tar= tarfile.open(dst_tarfile, 'w:gz')
   for src_file in src_files:
      tar.add(src_file, arcname=os.path.basename(src_file)) #tar.addfile(tarfile.TarInfo("myfilename.txt"), file("/path/to/filename.txt"))
   tar.close
   print('compressing OVF files FINISHED '+dst_tarfile)


   print('REMOVING OVF '+src_dir)
   for src_file in src_files:
      os.remove(src_file)

   ###3.3 Compressing up - compress png dirs
   src_dir = results_dir_PlotsM
   src_files = glob.iglob(os.path.join(src_dir, '*.png'))
   dst_tarfile=results_dir_Plots+aux[0:len(aux)-1]+'.tgz'

   print('compressing PNG files STARTED '+src_dir)

   tar= tarfile.open(dst_tarfile, 'w:gz')
   for src_file in src_files:
      tar.add(src_file, arcname=os.path.basename(src_file)) #tar.addfile(tarfile.TarInfo("myfilename.txt"), file("/path/to/filename.txt"))
   tar.close
   print('compressing PNG files FINISHED '+dst_tarfile)
   """
   #shutil.rmtree(src_dir)
   print('plotting_m FINISHED')

   return 0



def cleaning(sim_mx3, Damping, Ex, M, D, K1, K2, B_Max, size_X, size_Y, size_Z, Nx, Ny, Nz, Method,  i, aux, working_dir, results_dir_Plots, results_dir_PlotsM, results_dir_Data, src_dir, src_data, dst_data):


   print('cleaning STARTED')
   #### 0.0 dir definitions
                   
   if not os.path.exists(working_dir): 
      print(working_dir+" doesn't exist")

   if not os.path.exists(results_dir_Data): 
            print(results_dir_Data+" doesn't exist")
   print(results_dir_Data)

   if not os.path.exists(results_dir_Plots): 
      print(results_dir_Plots+" doesn't exist")
   print(results_dir_Plots)

   if not os.path.exists(results_dir_PlotsM): 
      print(results_dir_PlotsM+" doesn't exist")
   print(results_dir_PlotsM)

   ### 2.0 Some data analysis and plotting don't remember why the conversion to .dat
   #src_file=os.path.join(working_dir,sim_mx3+".out/table.txt")
   #dst_file=os.path.join(working_dir,sim_mx3+".out/"+sim_mx3+".dat")

   if os.path.isfile(src_data):
      os.rename(src_data,dst_data) #      copyfile(src_file,dst_file) #      os.remove(src_file)
      print(src_data+" renamed into "+dst_data)
   elif (os.path.isfile(dst_data)) and (not os.path.isfile(src_data)):
      print(dst_data+" already exists")
   elif (not os.path.isfile(dst_data)) and (not os.path.isfile(src_data)):
      print(src_data+" doesn't exist")

   ### 2.1 copy the data tables into the data dir 
   #src_file=os.path.join(working_dir,sim_mx3+".out/"+sim_mx3+".dat")
   dst_file=os.path.join(working_dir,results_dir_Data+"v"+sim_mx3+"_"+str("%.3f"%Damping)+"_"+str(Ex)+"_"+str("%.2f"%B_Max)+"_"+str("%.2f"%M)+"_"+str(D)+"_"+str("%.2f"%K1)+'_'+str("%.2f"%K2)+"_"+str(size_X)+"_"+str(size_Y)+"_"+str(size_Z)+"_"+str(re.search(r"(\w+)", Method).group(1))+'.dat')

   if os.path.isfile(dst_data):
      shutil.copyfile(dst_data,dst_file) #      os.remove(src_file) #os.rename(src_file,dst_file)
      print(dst_data+" renamed into "+dst_file)
   else:
      print(dst_data+" doesn't exist")


   print('cleaning FINISHED')

   return 0


def plotting(sim_mx3, Damping, Ex, M, D, K1, K2, B_Max, size_X, size_Y, size_Z, Nx, Ny, Nz, Method, i, aux, working_dir,results_dir_Plots, results_dir_PlotsM, results_dir_Data, src_dir, src_data, dst_data):


   print('Plotting Started')
   #### 0.0 dir creation
                   
   if not os.path.exists(working_dir): 
      print(working_dir+' doesn\'t exist')

   if not os.path.exists(results_dir_Data): 
      print(results_dir_Data+' doesn\'t exist')
   print(results_dir_Data)

   if not os.path.exists(results_dir_Plots): 
      print(results_dir_Plots+' doesn\'t exist')
   print(results_dir_Plots)

   ### 2.0 Some data analysis and plotting don't remember why the conversion to .dat
   #src_file=os.path.join(working_dir,sim_mx3+'.out/'+sim_mx3+'.dat')

   if not os.path.isfile(dst_data):
      print(dst_data+' doesn\'t exist')

   """
   # 2.1 reading the using pandas, because why not

#   i=0
#   columnUnits=[]
#   with open(dst_file,'r') as f:
#      for line in f:
#         if i == 0:
#            columnNames=re.findall('[\w]+([\w])', line)  #print(columnNames)
#            columnNames=re.findall('[\w]+:[\w]*:[\w /]*(?![\w]*:)', line)  #print(columnNames)
#         elif i == 4:
#            for word in line.split():
#               columnUnits.append("["+word+"]") # 4th line #print(columnUnits[2:len(columnUnits)])
#         elif i > 0:
#            break
#         i+=1
#   f.close()
#   print(columnNames)

   #here i0=2 because we skip the first two elements, since they are a '#' and the string 'Units:'
#   zipped=(zip(columnNames,columnUnits[2:len(columnUnits)]))

#   listHeader=[]
#   for i in zipped:
#      listHeader.append(' '.join(i))

   """

   #3.0 plotting is real

   Data_Hyst = pd.read_csv(dst_data, sep='[\s]+', engine='python', comment='#', header=None, skiprows=2) #, skiprows=38, skipfooter=1,  comment='#' ,comment='#' skiprows=1, 
   Data_Hyst.columns = ['t (s)', 'mx ()', 'my ()', 'mz ()', 'B_extx (T)', 'B_exty (T)', 'B_extz (T)', 'E_total (J)', 'ext_topologicalcharge ()', 'i ()'] 

   #3.1 for M(H)
   x_list=[]
   y_list=[]
   branch_number=3
   branch_min=1
   x0_HystBranch = int(len(Data_Hyst.index)/branch_number)

   for i in np.arange(branch_min,branch_number):
      x_list.append(Data_Hyst.ix[i*x0_HystBranch:(i+1)*x0_HystBranch,'B_extz (T)'].values)
      y_list.append(Data_Hyst.ix[i*x0_HystBranch:(i+1)*x0_HystBranch, 'mz ()'].values)


   K1    = float(K1)
   K2    = float(K2)
   M     = float(M)
   K_eff = (K1*1e6 - (4.0*np.pi*1e-7)*((M*1e6)**2)/2.0)*1e-6
      
   fig_name=results_dir_Plots+'M_vs_H_v'+sim_mx3+'_'+str(Damping)+'_'+str(Ex)+'_'+str(B_Max)+'_'+str(M)+'_'+ str(D)+'_'+str("%.2f"%K1)+'_'+str("%.2f"%K_eff)+'_'+str("%.2f"%K2)+'_'+str(size_X)+'_'+str(size_Y)+'_'+str(size_Z)+'_'+str(re.search(r"(\w+)", Method).group(1))+".png"
   plot_title='B<sub>Max</sub> '+str(B_Max)+'[T]; Ex '+str(Ex)+'[pJ/m]; M '+str(M)+'[MA/m]; D '+str(D)+'[mJ/m<sup>2</sup>]; K<sub>eff</sub> '+str("%.2f"%K_eff)+'[MJ/m<sup>3</sup>]; K1 '+str("%.2f"%K1)+'[MJ/m<sup>3</sup>]; <br>K2 '+str(K2)+'[MJ/m<sup>3</sup>];  X '+str(size_X)+"[nm] Y "+str(size_Y)+"[nm] Z "+str(size_Z)+"[nm]"
   curves_names=['$\\text{B}_{min} \\rightarrow \\text{B}_{max}$',  '$\\text{B}_{max} \\rightarrow \\text{B}_{min}$']

   x_Axis_Title= '$\\text{B [T]}$'
   y_Axis_Title= '$\\text{M}_{z}$'

   plot_symbols=['dot','triangle-down','diamond']
   plot_modes=['lines+markers','lines+markers','lines+markers']

   AxisFontSize=32
   PlotFontSize=25
   TicksFontSize=24


   Plot_Hyst(fig_name, x_list, y_list, plot_title, curves_names, x_Axis_Title, y_Axis_Title, x_range, y_range, plot_symbols, plot_modes, AxisFontSize, PlotFontSize, TicksFontSize) 

   #3.2 for SKNumber(H)
   x_list=[]
   y_list=[]
   branch_number=3
   branch_min=1
   x0_HystBranch = int(len(Data_Hyst.index)/branch_number)

   for i in np.arange(branch_min,branch_number):
      aux_x=Data_Hyst.ix[i*x0_HystBranch:(i+1)*x0_HystBranch,'B_extz (T)'].values
      aux_y=Data_Hyst.ix[i*x0_HystBranch:(i+1)*x0_HystBranch, 'ext_topologicalcharge ()'].values
      x_list.append(aux_x[int(len(aux_x)/2)+1:len(aux_x)])
      y_list.append(abs(aux_y[int(len(aux_y)/2)+1:len(aux_y)]))

      #x_list.append(Data_Hyst.ix[i*x0_HystBranch:(i+1)*x0_HystBranch,'B_extz (T)'].values)
      #y_list.append(Data_Hyst.ix[i*x0_HystBranch:(i+1)*x0_HystBranch, 'ext_topologicalcharge ()'].values)


   fig_name=results_dir_Plots+'SkNumber_vs_H_v'+sim_mx3+'_'+str(Damping)+'_'+str(Ex)+'_'+str(B_Max)+'_'+str(M)+'_'+ str(D)+'_'+str("%.2f"%K1)+'_'+str("%.2f"%K_eff)+'_'+str("%.2f"%K2)+'_'+str(size_X)+'_'+str(size_Y)+'_'+str(size_Z)+'_'+str(re.search(r"(\w+)", Method).group(1))+".png"
   plot_title='B<sub>Max</sub> '+str(B_Max)+'[T]; Ex '+str(Ex)+'[pJ/m]; M '+str(M)+'[MA/m]; D '+str(D)+'[mJ/m<sup>2</sup>]; K<sub>eff</sub> '+str("%.2f"%K_eff)+'[MJ/m<sup>3</sup>]; K1 '+str("%.2f"%K1)+'[MJ/m<sup>3</sup>]; <br>K2 '+str(K2)+'[MJ/m<sup>3</sup>]; X '+str(size_X)+"[nm] Y "+str(size_Y)+"[nm] Z "+str(size_Z)+"[nm]"
   #curves_names=['$\\text{B}_{max} \\rightarrow \\text{B}_{min} $',  '$\\text{B}_{min} \\rightarrow \\text{B}_{max}$',  '$\\text{B}_{max} \\rightarrow \\text{B}_{min}$']
   curves_names=['$\\text{B}_{min} \\rightarrow \\text{B}_{max}$',  '$\\text{B}_{max} \\rightarrow \\text{B}_{min}$']

   x_Axis_Title= '$\\text{B [T]}$'
   y_Axis_Title= '$\\text{Sk}_{\#}$'

   plot_symbols=['dot','triangle-down','diamond']
   plot_modes=['lines+markers','lines+markers','lines+markers']

   AxisFontSize=32
   PlotFontSize=25
   TicksFontSize=24

   if D == 2:
      x_range=[-2,2]
      y_range=[0,900] #'Automatic'
   else :
      x_range=[-4,4]
      y_range=[0,900] #'Automatic'

   Plot_Hyst(fig_name, x_list, y_list, plot_title, curves_names, x_Axis_Title, y_Axis_Title, x_range, y_range, plot_symbols, plot_modes, AxisFontSize, PlotFontSize, TicksFontSize)

   return 0 


def copy_tgz_file(sim_mx3, Damping, Ex, M, D, K1, K2, B_Max, size_X, size_Y, size_Z, Nx, Ny, Nz, Method, i, aux, working_dir, results_dir_Plots, results_dir_PlotsM, results_dir_Data, dst_dir, src_data, dst_data):

   src_file = working_dir+'0.4_SK.out/ovf_files.tgz'

   if not os.path.exists(src_file):
      print('src_file doesn\'t exist'+src_file)
   
   if not os.path.exists(dst_dir):
      print('dst_dir doesn\'t exist'+dst_dir)
      
   name_file=aux.replace('/','_')[0:len(aux.replace('/','_'))-1]+'.tgz'
   dst_file = dst_dir+name_file
   print('src_file: '+src_file)
   print('dst_file: '+dst_file)

   if os.path.exists(dst_file):
      os.remove(dst_file)
      print('pre-existing dst_file, thus deleted')
      
   shutil.copyfile(src_file, dst_file)

   return 0


def re_plotting(exception_file):

   #we exported to an exception file because didn't want to deal with threads+I/O 
   #exception_file.seek(0)
   exception_list = []

   for line in exception_file:
      inner_list = [elt.strip() for elt in line.split('\t')] #splitting to and excluding last element ''
      exception_list.append(inner_list[0:len(inner_list)-1])
   while len(exception_list) >= 1:
      try:
         print('Length except_connection_list '+str(len(exception_list)))
         plotting(*exception_list[0])
         exception_list.pop(0)
      except Exception as error:
         print(error)
         print('sleeping before executing %s'% exception_list[0][0:9])
         countdown(5)
         pass          #continue
      #else:
      #   del exception_list[0]

   return 0

def worker(i, q, function, exception_file):

   print('running worker: '+str(i))
   while True:
      item = q.get()  #print(q.get())
      #print(item)
      if item is None:
         break

      try:
         function(*item)
      except FileNotFoundError as error_file:
         print('plotting '+str(error_file))
         #pass
      except NameError as error_name:
         print('plotting '+str(error_name))
         #exception_file.write('name/variable not found %s \n'% parameters[17])
         #pass
      except Exception as error_connection:# TimeoutError as error_connection:
         print('Error found in '+str(function)+' '+str(error_connection))
         #print(item)
         for i in item:
            exception_file.write(str(i)+'\t') # is the dst_data file, it has the datatable
         exception_file.write('\n') # is the dst_data file, it has the datatable
         pass
      #except Exception:
      #   pass

      #print(function)
      q.task_done()


def threading_function(q, num_worker_threads, worker_name, function, parameter_list, exception_file):
#q = queue.Queue()
#num_worker_threads = 4
   threads = []
   for i in range(num_worker_threads):
      t = threading.Thread(target=worker_name, args=(i,q,function, exception_file))
      t.start()#q has nothing at this moment? -> 
      #t.daemon = True
      threads.append(t)

   for parameters in parameter_list:
      q.put(parameters)

   # block until all tasks are done
   #q.join() # didn't work with this, don't know how to solve it

   # stop workers
   for i in range(num_worker_threads):
      q.put(None)
   for t in threads:
      t.join()

   return 0

###################### complete sims ######################
###################### complete sims ######################
###################### complete sims ######################


Damping=0
sim_mx3="0.7_SK"

exception_f            = os.path.join(initial, 'exception.txt')
exception_file         = open(exception_f, 'w')
except_connection_list = []
parameter_list         = []
working_dir_list       = []

n = 20        #number of magnetic layers
i_list = np.arange(0,1,1)
for Index in i_list : #Index is the seed
   #input of material parameters
   
   #sample-  D,   Ex,  size_Z, Ms, K1_eff, K2, Bmin,   Bmax, Bstep
   #   11b-[1.8, 11.83, 0.8,  0.87,  0.25, 0.0, 0.0000,0.0700,0.0050]#why is the bmax >11f?
   #   11c-[_,    _,    0.9,  1.02,  0.15, 0.0, 0.0000,0.0400,0.0050]
   #   11d-[_,    _,    1.0,  1.10,  0.05, 0.0, 0.0000,0.0400,0.0050]
   #   11e-[_,    _,    0.8,  0.84,  0.11, 0.0, 0.0000,0.0400,0.0050]
   #   11f-[1.95,11.04, 1.0,  1.05, 0.005, 0.0, 0.0000,0.0400,0.0050]
   #   11p-[1.99,11.72, 0.7,  0.826,0.225, 0.0, 0.0000,0.0400,0.0050]
   mat_name_list=['11b']
   mat_param=np.array([[1.99,11.72, 0.7,  0.826,0.225, 0.0, 0.0000,0.0400,0.0050]])
   D_var =np.array([0.8,1.0,1.2])#*D_op
   Ex_var=np.array([0.8,1.0,1.2])#*Ex_op
   #input size of dots
   dot_size_var=np.array([152, 200, 252, 300, 352, 400]) #diameter(nm) 
   [mat_len,_]=np.shape(mat_param)
   
   for mat_i in range(mat_len): #no idea what this is - MJ/m^3  #[0.0, -K1*0.1, -K1*0.2, -K1*0.3, -K1*0.4 1 Merg/cc = 10^6 10^-7/10^-6=10^5 = 0.1*10^6 J/m^3
        
        mat_name=mat_name_list[mat_i]
        D_op=mat_param[mat_i,0]
        Ex_op=mat_param[mat_i,1]
        size_Z=mat_param[mat_i,2]
        M=mat_param[mat_i,3]
        K1_eff=mat_param[mat_i,4]
        K2=mat_param[mat_i,5]
        B_Min=mat_param[mat_i,6]
        B_Max=mat_param[mat_i,7]
        B_Step=mat_param[mat_i,8]
        #sims parameters
        nrepeats=n #number of magnetic layers
        Nz=4*n-2 #number of total layers = magnetic + non-magnetic
        cell_size=4 #nm per cell
        N_cell_var=dot_size_var/cell_size #calculation of number cells
    
        
        Method='relax()' #Method='run(5e-9)'
        
        #calculation of K from Keff and Ms
        K1 = (K1_eff*1e6 + (4.0*np.pi*1e-7)*((M*1e6)**2)/2.0)*1e-6
        
        #Getting combination of Ex and D pairs
        D_var      =  np.round(D_var*D_op,2)#[1.50, 1.70, 1.90]
        Ex_var     =  np.round(Ex_var*Ex_op,2) #[1.50, 1.70, 1.90]  
        Var_comb  =  np.array(list(product(N_cell_var,D_var,Ex_var)))
        
        N_cell_list=Var_comb[:,0]
        D_list=Var_comb[:,1]
        Ex_list=Var_comb[:,2]
    
        for N_cell, Ex, D in zip(N_cell_list,Ex_list, D_list): #MJ/m^3  # 1 Merg/cc = 10^6 10^-7/10^-6=10^5 = 0.1*10^6 J/m^3
           
           Nx=int(N_cell)
           Ny=Nx
           size_X=Nx*cell_size
           size_Y=Ny*cell_size
           
           aux='v'+sim_mx3+'_'+str(nrepeats)+'_'+str('%.3f'%B_Max)+'_'+str(Ex)+'_'+str('%.2f'%D)+'/'+str('%.2f'%M)+'_'+str('%.2f'%K1_eff)+'_'+str('%.2f'%K1)+'_'+str('%.2f'%K2)+'_'+str(size_X)+'_'+str(size_Y)+'_'+str(size_Z)+'_'+str(re.search(r'(\w+)', Method).group(1))+'_'+str('%.0f'%Index)+'/'

           working_dir=os.path.join(initial,mat_name+'/'+aux)
           results_dir_Plots=os.path.join(initial,'v'+sim_mx3+'_Results/Plots/'+aux.replace('/','_')[0:len(aux.replace('/','_'))-1])
           results_dir_PlotsM=os.path.join(initial,'v'+sim_mx3+'_Results/Plots/'+aux)
           results_dir_Data=os.path.join(initial,'v'+sim_mx3+'_Results/')
                       
           src_dir=os.path.join(working_dir,sim_mx3+'.out/')
           src_data=os.path.join(working_dir,sim_mx3+'.out/table.txt')
           dst_data=os.path.join(working_dir,sim_mx3+'.out/'+sim_mx3+'.dat')

           working_dir_list.append(working_dir)
           
           parameter_list.append([sim_mx3, Damping, Ex, M, D, K1_eff, K1, K2, B_Min, B_Max, B_Step, size_X, size_Y, size_Z, Nx, Ny, Nz, nrepeats, Method, Index, aux, working_dir, results_dir_Plots, results_dir_PlotsM, results_dir_Data, src_dir, src_data, dst_data])

"""                                             
         for size_X, size_Y, Nx, Ny in zip([152, 200, 252, 300, 352], [152, 200, 252, 300, 352], [38, 50, 63, 75, 88], [38, 50, 63, 75, 88]):  
            #11p original
            D_list      = [1.59, 1.99, 2.39] #[1.90, 2.10, 2.30]

            size_Z_list = []
            M_list      = []
            K1_eff_list = []
            Ex_list     = []
            B_Min_list  = []
            B_Max_list  = []
            B_Step_list = []

            nrepeats_list = []       #number of magnetic layers
            for i in D_list:
               
               size_Z_list.append(0.70)
               M_list.append(0.826)     
               K1_eff_list.append(0.225)
               Ex_list.append(11.72)
               B_Min_list.append(0.0000)
               B_Max_list.append(0.0400)
               B_Step_list.append(0.0050)
               nrepeats_list.append(n)

            nrepeats_list = np.array(nrepeats_list)
            Nz_list       = 4*nrepeats_list-2  #number of total layers = magnetic + non-magnetic

            for B_Min, B_Max, B_Step, Ex, D, K1_eff, M, size_Z, nrepeats, Nz in zip(B_Min_list, B_Max_list, B_Step_list, Ex_list, D_list, K1_eff_list, M_list, size_Z_list, nrepeats_list, Nz_list): #MJ/m^3  # 1 Merg/cc = 10^6 10^-7/10^-6=10^5 = 0.1*10^6 J/m^3

               K1 = (K1_eff*1e6 + (4.0*np.pi*1e-7)*((M*1e6)**2)/2.0)*1e-6
               Method='relax()' #Method='run(5e-9)'

               aux='v'+sim_mx3+'_'+str(nrepeats)+'_'+str('%.3f'%B_Max)+'_'+str(Ex)+'_'+str('%.2f'%D)+'/'+str('%.2f'%M)+'_'+str('%.2f'%K1_eff)+'_'+str('%.2f'%K1)+'_'+str('%.2f'%K2)+'_'+str(size_X)+'_'+str(size_Y)+'_'+str(size_Z)+'_'+str(re.search(r'(\w+)', Method).group(1))+'_'+str('%.0f'%Index)+'/'

               working_dir=os.path.join(initial,aux)
               results_dir_Plots=os.path.join(initial,'v'+sim_mx3+'_Results/Plots/'+aux.replace('/','_')[0:len(aux.replace('/','_'))-1])
               results_dir_PlotsM=os.path.join(initial,'v'+sim_mx3+'_Results/Plots/'+aux)
               results_dir_Data=os.path.join(initial,'v'+sim_mx3+'_Results/')
                           
               src_dir=os.path.join(working_dir,sim_mx3+'.out/')
               src_data=os.path.join(working_dir,sim_mx3+'.out/table.txt')
               dst_data=os.path.join(working_dir,sim_mx3+'.out/'+sim_mx3+'.dat')

               working_dir_list.append(working_dir)
               parameter_list.append([sim_mx3, Damping, Ex, M, D, K1_eff, K1, K2, B_Min, B_Max, B_Step, size_X, size_Y, size_Z, Nx, Ny, Nz, nrepeats, Method, Index, aux, working_dir, results_dir_Plots, results_dir_PlotsM, results_dir_Data, src_dir, src_data, dst_data])

         for size_X, size_Y, Nx, Ny in zip([152, 200, 252, 300, 352], [152, 200, 252, 300, 352], [38, 50, 63, 75, 88], [38, 50, 63, 75, 88]):  
            #11p
            D_list      = [1.59, 1.99, 2.39] #[1.90, 2.10, 2.30]

            size_Z_list = []
            M_list      = []
            K1_eff_list = []
            Ex_list     = []
            B_Min_list  = []
            B_Max_list  = []
            B_Step_list = []

            nrepeats_list = []       #number of magnetic layers
            for i in D_list:
               
               size_Z_list.append(0.70)
               M_list.append(0.826)     
               K1_eff_list.append(0.225)
               Ex_list.append(14.06)
               B_Min_list.append(0.0000)
               B_Max_list.append(0.0400)
               B_Step_list.append(0.0050)
               nrepeats_list.append(n)

            nrepeats_list = np.array(nrepeats_list)
            Nz_list       = 4*nrepeats_list-2  #number of total layers = magnetic + non-magnetic

            for B_Min, B_Max, B_Step, Ex, D, K1_eff, M, size_Z, nrepeats, Nz in zip(B_Min_list, B_Max_list, B_Step_list, Ex_list, D_list, K1_eff_list, M_list, size_Z_list, nrepeats_list, Nz_list): #MJ/m^3  # 1 Merg/cc = 10^6 10^-7/10^-6=10^5 = 0.1*10^6 J/m^3

               K1 = (K1_eff*1e6 + (4.0*np.pi*1e-7)*((M*1e6)**2)/2.0)*1e-6
               Method='relax()' #Method='run(5e-9)'

               aux='v'+sim_mx3+'_'+str(nrepeats)+'_'+str('%.3f'%B_Max)+'_'+str(Ex)+'_'+str('%.2f'%D)+'/'+str('%.2f'%M)+'_'+str('%.2f'%K1_eff)+'_'+str('%.2f'%K1)+'_'+str('%.2f'%K2)+'_'+str(size_X)+'_'+str(size_Y)+'_'+str(size_Z)+'_'+str(re.search(r'(\w+)', Method).group(1))+'_'+str('%.0f'%Index)+'/'

               working_dir=os.path.join(initial,aux)
               results_dir_Plots=os.path.join(initial,'v'+sim_mx3+'_Results/Plots/'+aux.replace('/','_')[0:len(aux.replace('/','_'))-1])
               results_dir_PlotsM=os.path.join(initial,'v'+sim_mx3+'_Results/Plots/'+aux)
               results_dir_Data=os.path.join(initial,'v'+sim_mx3+'_Results/')
                           
               src_dir=os.path.join(working_dir,sim_mx3+'.out/')
               src_data=os.path.join(working_dir,sim_mx3+'.out/table.txt')
               dst_data=os.path.join(working_dir,sim_mx3+'.out/'+sim_mx3+'.dat')

               working_dir_list.append(working_dir)
               parameter_list.append([sim_mx3, Damping, Ex, M, D, K1_eff, K1, K2, B_Min, B_Max, B_Step, size_X, size_Y, size_Z, Nx, Ny, Nz, nrepeats, Method, 
               , aux, working_dir, results_dir_Plots, results_dir_PlotsM, results_dir_Data, src_dir, src_data, dst_data])
               


         for size_X, size_Y, Nx, Ny in zip([152, 200, 252, 300, 352], [152, 200, 252, 300, 352], [38, 50, 63, 75, 88], [38, 50, 63, 75, 88]):  
            #11p
            D_list      = [1.59, 1.99, 2.39] #[1.90, 2.10, 2.30]
            size_Z_list = []
            M_list      = []
            K1_eff_list = []
            Ex_list     = []
            B_Min_list  = []
            B_Max_list  = []
            B_Step_list = []

            nrepeats_list = []       #number of magnetic layers
            for i in D_list:

               size_Z_list.append(0.70)
               M_list.append(0.826)     
               K1_eff_list.append(0.225)               
               Ex_list.append(9.38)
               B_Min_list.append(0.0000)
               B_Max_list.append(0.0400)
               B_Step_list.append(0.0050)
               nrepeats_list.append(n)               

            nrepeats_list = np.array(nrepeats_list)
            Nz_list       = 4*nrepeats_list-2  #number of total layers = magnetic + non-magnetic

            for B_Min, B_Max, B_Step, Ex, D, K1_eff, M, size_Z, nrepeats, Nz in zip(B_Min_list, B_Max_list, B_Step_list, Ex_list, D_list, K1_eff_list, M_list, size_Z_list, nrepeats_list, Nz_list): #MJ/m^3  # 1 Merg/cc = 10^6 10^-7/10^-6=10^5 = 0.1*10^6 J/m^3

               K1 = (K1_eff*1e6 + (4.0*np.pi*1e-7)*((M*1e6)**2)/2.0)*1e-6
               Method='relax()' #Method='run(5e-9)'

               aux='v'+sim_mx3+'_'+str(nrepeats)+'_'+str('%.3f'%B_Max)+'_'+str(Ex)+'_'+str('%.2f'%D)+'/'+str('%.2f'%M)+'_'+str('%.2f'%K1_eff)+'_'+str('%.2f'%K1)+'_'+str('%.2f'%K2)+'_'+str(size_X)+'_'+str(size_Y)+'_'+str(size_Z)+'_'+str(re.search(r'(\w+)', Method).group(1))+'_'+str('%.0f'%Index)+'/'

               working_dir=os.path.join(initial,aux)
               results_dir_Plots=os.path.join(initial,'v'+sim_mx3+'_Results/Plots/'+aux.replace('/','_')[0:len(aux.replace('/','_'))-1])
               results_dir_PlotsM=os.path.join(initial,'v'+sim_mx3+'_Results/Plots/'+aux)
               results_dir_Data=os.path.join(initial,'v'+sim_mx3+'_Results/')
                           
               src_dir=os.path.join(working_dir,sim_mx3+'.out/')
               src_data=os.path.join(working_dir,sim_mx3+'.out/table.txt')
               dst_data=os.path.join(working_dir,sim_mx3+'.out/'+sim_mx3+'.dat')

               working_dir_list.append(working_dir)
               parameter_list.append([sim_mx3, Damping, Ex, M, D, K1_eff, K1, K2, B_Min, B_Max, B_Step, size_X, size_Y, size_Z, Nx, Ny, Nz, nrepeats, Method, Index, aux, working_dir, results_dir_Plots, results_dir_PlotsM, results_dir_Data, src_dir, src_data, dst_data])
"""
         
q = queue.Queue()
num_worker_threads = 4

threading_function(q, num_worker_threads, worker, writting_mx3, parameter_list, exception_file)

simulation_command = mumax_command #+' -gpu 0 '                           
for sim_file in working_dir_list:
   simulation_command+=' '+sim_file+sim_mx3+'.mx3'
### 1.0 Running the sims
### mumax takes care of running in parallel different sims
#print(simulation_command)
#os.system(simulation_command)    

#threading_function(q, num_worker_threads, worker, plotting_m, parameter_list, exception_file)
#threading_function(q, num_worker_threads, worker, cleaning, parameter_list, exception_file)
#threading_function(q, num_worker_threads, worker, plotting, parameter_list, exception_file)
#exception_file.close() #good manners, always

#exception_file         = open(exception_f, 'r')
#re_plotting(exception_file)
#exception_file.close() #good manners, always


def extracting_png():

   sim_mx3='0.4_SK'
   results_dir_Plots=os.path.join(initial,'v'+sim_mx3+'_Results/Plots/')

   if not os.path.exists(results_dir_Plots): 
      print(results_dir_Plots+' doesn\'t exist')
   #print(results_dir_Plots)


   src_dir=results_dir_Plots
   files = glob.iglob(os.path.join(src_dir, '*.tgz'))
   for fname in files:
      #dst_dir=os.path.splitext(fname)[0]
      #print(dst_dir)
      #os.makedirs(dst_dir)

      print('Extracting '+ fname)
      tar = tarfile.open(fname)
      tar.extractall(path=src_dir)
      tar.close()
     
#extracting_png()

"""


#    executable=open(os.path.join(initial,str(Keff_Fl),'main.sh'), 'w+')

#    os.system(str(tclsh_command)+' '+ str(oommf_command)+str( +fg +err errfile_plot +out outfile_plot avf2ppm -cwd )+working_dir+"/ -f -ipat *.omf   -opatsub .gif -filter \"ppmquant 256 | ppmtogif\" ")
#    executable.close()
#    os.chmod(working_dir+"/main.sh", 0o777)
#    os.system(working_dir+"/./main.sh")    
#executable=open(os.path.join(initial,"main.sh"), 'w+')
#executable.close()

#        executable=open(working_dir+"main.sh", 'r+')
#    shutil.copyfile(mif_file,a.text)

#   if os.path.exists(initial+"/"+sim_mif+".mif"):
#      shutil.copy(initial+"/"+sim_mif+".mif",working_dir)
#   else:
#      print("sim_mif doesn't exist")
#   for file in os.listdir(m0_dir):
#      if file.endswith(".omf"):
#   initial_m=os.path.join(m0_dir,file)


"""       
"""
 cmake -D CMAKE_BUILD_TYPE=RELEASE -D CMAKE_INSTALL_PREFIX=/usr/local -D INSTALL_C_EXAMPLES=ON -D INSTALL_PYTHON_EXAMPLES=ON -D OPENCV_EXTRA_MODULES_PATH=~/opencv_contrib/modules -D BUILD_EXAMPLES=ON -D WITH_CUDA=OFF
"""

