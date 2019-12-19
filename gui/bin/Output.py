import wx
import wx.grid
import wx.lib.scrolledpanel as scrolled
import numpy as np

from copy import deepcopy

from matplotlib.figure import Figure
from matplotlib.backends.backend_wxagg import \
    FigureCanvasWxAgg as FigCanvas, \
    NavigationToolbar2WxAgg as NavigationToolbar

from matplotlib.path import Path
import matplotlib.patches as patches

class OutputFrame(wx.Frame):
    title = 'DS-Frame2D - Output'
    def __init__(self, Lim_X1, Lim_Y1, Lim_X2, Lim_Y2, NOD, ELE, BC, SEC, FOR1, FOR2, FOR3, accuracy, PATH):
        self.Lim_X1 = Lim_X1
        self.Lim_Y1 = Lim_Y1
        self.Lim_X2 = Lim_X2
        self.Lim_Y2 = Lim_Y2
        
        self.Dim_X = self.Lim_X2 - self.Lim_X1
        self.Dim_Y = self.Lim_Y2 - self.Lim_Y1

        self.NOD  = np.asarray(NOD)
        self.ELE  = np.asarray(ELE)
        self.BC   = np.asarray(BC)
        self.SEC  = np.asarray(SEC)
        self.FOR1 = np.asarray(FOR1)
        self.FOR2 = np.asarray(FOR2)
        self.FOR3 = np.asarray(FOR3)

        self.accuracy = accuracy
        
        self.PATH = PATH
        
        wx.Frame.__init__(self, None, -1, self.title)
        
        Pth = self.PATH + '/output.txt'
        f = open(Pth)
        lines = f.readlines()
        f.close()
        for i in range(0,len(lines)):
            k = list(map(str, lines[i].split()))
            if (k[0] == 'NODAL') and (k[1] == 'FORCES'):
                P_nod = i+1
            elif (k[0] == 'NODAL') and (k[1] == 'DISPLACEMENTS'):
                D_nod = i+1
            elif (k[0] == 'ELEMENT') and (k[1] == 'FORCES'):
                P_el = i+1
        
        self.P_nodal = []
        for i in range(P_nod+1,D_nod-1):
            self.P_nodal.append(list(map(float, lines[i].split())))
        
        self.D_nodal = []
        for i in range(D_nod+1,P_el-1):
            self.D_nodal.append(list(map(float, lines[i].split())))
        
        self.P_element = []
        for i in range(P_el+1,len(lines)):
            self.P_element.append(list(map(float, lines[i].split())))
        
        self.P_nodal = np.asarray(self.P_nodal)
        self.D_nodal = np.asarray(self.D_nodal)
        self.P_element = np.asarray(self.P_element)
        
        Pth = self.PATH + '/diagrams.txt'
        f = open(Pth)
        lines = f.readlines()
        f.close()
        
        self.diagrams = np.zeros((len(lines),len(lines[0].split())))
        for i in range(0,len(lines)):
            k = list(map(float, lines[i].split()))
            self.diagrams[i,:] = k

        Pth = self.PATH + '/eigenvectors.txt'
        f = open(Pth)
        lines = f.readlines()
        f.close()
        
        self.eigenvectors = np.zeros((len(lines),len(lines[0].split())))
        for i in range(0,len(lines)):
            k = list(map(float, lines[i].split()))
            self.eigenvectors[i,:] = k

        Pth = self.PATH + '/eigen_diagrams.txt'
        f = open(Pth)
        lines = f.readlines()
        f.close()
        
        self.eigen_diagrams = np.zeros((len(lines),len(lines[0].split())))
        for i in range(0,len(lines)):
            k = list(map(float, lines[i].split()))
            self.eigen_diagrams[i,:] = k
        
        self.create__panel()

    def create__panel(self):

        self.dpi = 100

        self.panel = wx.Panel(self)
        
        self.nb = wx.Notebook(self.panel)
        
        self.DeformedTab = wx.Window(self.nb)
        self.ReactionTab = wx.Window(self.nb)
        self.AxialTab = wx.Window(self.nb)
        self.ShearTab = wx.Window(self.nb)
        self.MomentTab = wx.Window(self.nb)
        self.ModeTad = wx.Window(self.nb)
        
        self.nb.AddPage(self.DeformedTab, "Deformed Shape")
        self.nb.AddPage(self.ReactionTab, "Reaction Forces")
        self.nb.AddPage(self.AxialTab, "Axial Forces")
        self.nb.AddPage(self.ShearTab, "Shear Forces")
        self.nb.AddPage(self.MomentTab, "Bending Moments")
        self.nb.AddPage(self.ModeTad, "Mode Shapes")
        
        self.Deformed_shape(self.nb)  
        self.Reaction_forces(self.nb)
        self.Axial_forces(self.nb)
        self.Shear_forces(self.nb)
        self.Bending_moments(self.nb)
        self.Mode_shapes(self.nb)
        
        self.vbox = wx.BoxSizer(wx.VERTICAL)
        self.vbox.Add(self.nb, 1, wx.LEFT | wx.TOP | wx.GROW)
        
        self.panel.SetSizer(self.vbox)
        self.vbox.Fit(self)

    def Deformed_shape(self, parent):   
        self.figD = Figure((6.0, 5.0), dpi=self.dpi, tight_layout = True, facecolor = 'white')
        self.canvasD = FigCanvas(self.DeformedTab, -1, self.figD)
        
        self.axesD = self.figD.add_subplot(111, aspect='equal', frame_on = False)
        self.axesD.set_xlim([self.Lim_X1, self.Lim_X2])
        self.axesD.set_ylim([self.Lim_Y1, self.Lim_Y2])  

        self.canvasD.mpl_connect('pick_event', self.on_pick)

        self.cb_gridD = wx.CheckBox(self.DeformedTab, -1, "Show Grid", style=wx.ALIGN_RIGHT)
        self.cb_gridD.SetValue(True)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_gridD, self.cb_gridD)

        self.toolbar = NavigationToolbar(self.canvasD)

        self.vbox2 = wx.BoxSizer(wx.VERTICAL)
        quote = wx.StaticText(self.DeformedTab, -1, "Scale factor: ")
        self.vbox2.Add(quote, 0, wx.ALIGN_CENTRE|wx.ALL, 5)

        self.ys_l = np.zeros((self.accuracy,len(self.ELE[:,0])))
        maxval = 0
        for i in range(0,len(self.ELE[:,0])):
            self.ys_l[:,i] = self.diagrams[:,5*i + 4]#np.loadtxt('diagrams.txt', usecols=range(5*i + 4,5*i + 5))
            for j in range(0,self.accuracy):
                if abs(self.ys_l[j][i]) > maxval:
                    maxval = abs(self.ys_l[j][i])
        
        if maxval == 0:
            maxval = 1
        self.Scale_factor = int(self.Dim_X/20/maxval)

        text = wx.TextCtrl(self.DeformedTab, -1, str(self.Scale_factor), size=(80,-1))
        self.Bind(wx.EVT_TEXT, self.EvtSc, text)
        self.vbox2.Add(text, 1, wx.ALIGN_CENTRE|wx.ALL,5)

        self.gridDeformed = wx.grid.Grid(self.DeformedTab, size=(323,300))

        self.gridDeformed.EnableEditing(False)
        self.gridDeformed.SetDefaultCellAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
        self.gridDeformed.CreateGrid(len(self.NOD[:,0]),4)
        self.gridDeformed.SetRowLabelSize(3)
        self.gridDeformed.SetColLabelValue(0, 'Node ID')
        self.gridDeformed.SetColLabelValue(1, 'dx')
        self.gridDeformed.SetColLabelValue(2, 'dy')
        self.gridDeformed.SetColLabelValue(3, 'd'u'\N{greek small letter phi}')

        for i in range(0,len(self.NOD[:,0])):
            self.gridDeformed.SetCellValue(i,0,str(i+1))
            self.gridDeformed.SetCellValue(i,1,'%.2E'%self.D_nodal[i][1])
            self.gridDeformed.SetCellValue(i,2,'%.2E'%self.D_nodal[i][2])
            self.gridDeformed.SetCellValue(i,3,'%.2E'%self.D_nodal[i][3])
       
        self.vbox2.AddSpacer(30)
        self.vbox2.Add(self.gridDeformed, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
        
        self.hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox1.Add(self.toolbar, 0, wx.EXPAND)
        flags = wx.ALIGN_RIGHT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        self.hbox1.AddSpacer(20)
        self.hbox1.Add(self.cb_gridD, 0, border=3, flag=flags)
        
        self.vbox1 = wx.BoxSizer(wx.VERTICAL)
        self.vbox1.Add(self.canvasD, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox1.Add(self.hbox1, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        
        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.vbox1, 0, border=3)#, flag=flags)
        self.hbox2.Add(self.vbox2, 0, border=3)#, flag=flags)
        
        self.DeformedTab.SetSizer(self.hbox2)
        
        self.draw_by_scale()
                
    def EvtSc(self,event):
        self.Scale_factor = float(event.GetString())
        self.draw_by_scale()

    def draw_by_scale(self):
        self.NOD_D = deepcopy(self.NOD)
        for i in range(0,len(self.NOD[:,0])):
            self.NOD_D[i][1] = self.NOD_D[i][1] + self.Scale_factor*float(self.D_nodal[i][1])
            self.NOD_D[i][2] = self.NOD_D[i][2] + self.Scale_factor*float(self.D_nodal[i][2])
        
        self.draw_figure(self.NOD,'grey',True,self.canvasD,self.axesD,self.cb_gridD,True,True)
        self.draw_figure(self.NOD_D,'black',False,self.canvasD,self.axesD,self.cb_gridD,True,False)
        
        for i in range(0,len(self.ELE[:,0])):
            [phi] = self.get_element_properties(i,['phi'])
            
            x1_ = self.NOD_D[self.ELE[i][1]-1][1]
            y1_ = self.NOD_D[self.ELE[i][1]-1][2]
            
            x2_ = self.NOD_D[self.ELE[i][2]-1][1]
            y2_ = self.NOD_D[self.ELE[i][2]-1][2]
 
            theta = np.arctan((y2_ - y1_)/(x2_ - x1_))
 
            Ld = ((x2_-x1_)**2 + (y2_-y1_)**2)**0.5
            Ld = Ld*abs(np.cos(phi-theta))
            delta_x = Ld/(self.accuracy-1)

            xs = np.zeros((self.accuracy-1) + 1)
            ys = np.zeros((self.accuracy-1) + 1)
            
            cy = self.ys_l[0][i]
            
            k = 0
            for j in range (0,self.accuracy):
                xs[k] = (j-1)*delta_x
                ys[k] = (self.ys_l[j][i] - cy)*self.Scale_factor
                x_new = np.cos(phi)*xs[k] - np.sin(phi)*ys[k] + x1_
                y_new = np.sin(phi)*xs[k] + np.cos(phi)*ys[k] + y1_
                xs[k] = x_new
                ys[k] = y_new
                k = k + 1
        
            self.axesD.plot(xs,ys,'black')
            
        self.canvasD.draw()
            
    def Reaction_forces(self, parent):
        self.figR = Figure((6.0, 5.0), dpi=self.dpi, tight_layout = True, facecolor = 'white')
        self.canvasR = FigCanvas(self.ReactionTab, -1, self.figR)
        
        self.axesR = self.figR.add_subplot(111, aspect='equal', frame_on = False)
        self.axesR.set_xlim([self.Lim_X1, self.Lim_X2])
        self.axesR.set_ylim([self.Lim_Y1, self.Lim_Y2])  

        self.canvasR.mpl_connect('pick_event', self.on_pick)

        self.cb_gridR = wx.CheckBox(self.ReactionTab, -1, "Show Grid", style=wx.ALIGN_RIGHT)
        self.cb_gridR.SetValue(True)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_gridR, self.cb_gridR)

        self.toolbar = NavigationToolbar(self.canvasR)
        
        self.hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox1.Add(self.toolbar, 0, wx.EXPAND)
        flags = wx.ALIGN_RIGHT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        self.hbox1.AddSpacer(20)
        self.hbox1.Add(self.cb_gridR, 0, border=3, flag=flags)
        
        self.vbox1 = wx.BoxSizer(wx.VERTICAL)
        self.vbox1.Add(self.canvasR, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox1.Add(self.hbox1, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        
        R_forces = []
        Nod_R = []
        for i in range(0,len(self.BC[:,0])):
            j = self.BC[i,0] - 1
            Nod_R.append(j+1)
            if len(self.FOR1) != 0:
                if len(self.FOR1) != 1 or self.FOR1[0][0] == 0:
                    if j+1 in self.FOR1[:,0]:
                        for k in range(0,len(self.FOR1[:,0])):
                            if j+1 == self.FOR1[k,0]:
                                R_forces.append([self.P_nodal[k][1]-self.FOR1[k][1], 
                                                 self.P_nodal[k][2]-self.FOR1[k][2], 
                                                 self.P_nodal[k][3]-self.FOR1[k][3]])
                    else:    
                        R_forces.append([self.P_nodal[j][1], self.P_nodal[j][2], self.P_nodal[j][3]])
                                
                else:    
                    R_forces.append([self.P_nodal[j][1], self.P_nodal[j][2], self.P_nodal[j][3]])
                            
            else:    
                R_forces.append([self.P_nodal[j][1], self.P_nodal[j][2], self.P_nodal[j][3]])

            if self.BC[i][1] == 0:
                R_forces[i][0] = 0
                 
            if self.BC[i][2] == 0:
                R_forces[i][1] = 0 
         
            if self.BC[i][3] == 0:
                R_forces[i][2] = 0

        self.gridReaction = wx.grid.Grid(self.ReactionTab, size=(323,300))
         
        self.gridReaction.EnableEditing(False)
        self.gridReaction.SetDefaultCellAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
        self.gridReaction.CreateGrid(len(self.BC[:,0]),4)
        self.gridReaction.SetRowLabelSize(3)
        self.gridReaction.SetColLabelValue(0, 'Node ID')
        self.gridReaction.SetColLabelValue(1, 'Rx')
        self.gridReaction.SetColLabelValue(2, 'Ry')
        self.gridReaction.SetColLabelValue(3, 'M')
 
        for i in range(0,len(Nod_R)):
            self.gridReaction.SetCellValue(i,0,str(Nod_R[i]))
            self.gridReaction.SetCellValue(i,1,'%.2f'%R_forces[i][0])
            self.gridReaction.SetCellValue(i,2,'%.2f'%R_forces[i][1])
            self.gridReaction.SetCellValue(i,3,'%.2f'%R_forces[i][2])

        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.vbox1, 0, wx.EXPAND)
        self.hbox2.Add(self.gridReaction, 0, wx.EXPAND)
         
 
        self.ReactionTab.SetSizer(self.hbox2)
         
        self.draw_figure(self.NOD,'black',True,self.canvasR,self.axesR,self.cb_gridR,False,True)
         
        max_val = 0
        for i in range(0,len(Nod_R)):
            for j in range(0,2):
                if abs(R_forces[i][j]) > max_val:
                    max_val = abs(R_forces[i][j])
        
        if max_val == 0:
            max_val = 1
          
        toler = max_val/1000
        scale = self.Dim_X/5/max_val
        for i in range(0,len(Nod_R)):
            p =[self.NOD[Nod_R[i]-1,1], self.NOD[Nod_R[i]-1,2]]
            if abs(R_forces[i][0]) > toler:
                sign = R_forces[i][0]/abs(R_forces[i][0])
                xyT = (p[0] + sign*self.Dim_X/20 + R_forces[i][0]*scale, p[1])
                self.axesR.annotate('',     xy = (p[0], p[1]),
                                   xytext     = xyT,
                                   arrowprops = dict(arrowstyle="<-",alpha=0.7, color='b'))
   
                self.axesR.text(xyT[0]+self.Dim_X/50,xyT[1]+self.Dim_X/100,
                                    'R'+str(int(Nod_R[i]))+'x', size='small',color='b',alpha=0.7 )
                  
            if abs(R_forces[i][1]) > toler:
                sign = R_forces[i][1]/abs(R_forces[i][1])
                xyT = (p[0], p[1] + sign*self.Dim_X/20 + R_forces[i][1]*scale)
                self.axesR.annotate('',     xy = (p[0], p[1]),
                                   xytext     = xyT,
                                   arrowprops = dict(arrowstyle="<-",alpha=0.7, color='b'))
   
                self.axesR.text(xyT[0]+self.Dim_X/100,xyT[1]+self.Dim_X/100,
                                    'R'+str(int(Nod_R[i]))+'y', size='small',color='b',alpha=0.7 )

            if abs(R_forces[i][2]) > toler:
                sign = R_forces[i][1]/abs(R_forces[i][1])
                xyT = (p[0] + sign*self.Dim_X/20, p[1] + sign*self.Dim_X/20)
                self.plot_arrow_moment('blue',p,0.03*self.Dim_X,0,(0,0))
                self.axesR.text(xyT[0]+self.Dim_X/100,xyT[1]+self.Dim_X/100,
                                    'M'+str(int(Nod_R[i])), size='small',color='b',alpha=0.7 )


    def plot_arrow_moment(self, col, position, scale, angle, tail):
        verts = [(0.0, 0.0  ), (1.5 , 1.5  ), (0.0, 3.0  ), (-1.5 , 1.5  ), (-1.3706, 1.9829), (-1.5 , 1.5  ), (-1.0171, 1.6294)]
        codes = [Path.MOVETO, Path.LINETO, Path.CURVE3, Path.CURVE3, Path.LINETO, Path.LINETO, Path.LINETO]
        vertsT = verts[:]
        for i in range(0,len(verts)):
            vertsT[i] = (verts[i][0]*np.cos(angle)-verts[i][1]*np.sin(angle), verts[i][0]*np.sin(angle)+verts[i][1]*np.cos(angle))

        for i in range(0,len(verts)):
            vertsT[i] = (scale*vertsT[i][0],scale*vertsT[i][1])
    
        for i in range(0,len(verts)):
            vertsT[i] = (vertsT[i][0] + position[0],vertsT[i][1] + position[1])

        pathT = Path(vertsT, codes)
        arrow = patches.PathPatch(pathT, facecolor='none', lw=1, edgecolor=col,alpha=0.8)
        self.axesR.add_patch(arrow)
        
        self.canvasR.draw()

    
    def Axial_forces(self, parent):   
        self.figA = Figure((6.0, 5.0), dpi=self.dpi, tight_layout = True, facecolor = 'white')
        self.canvasA = FigCanvas(self.AxialTab, -1, self.figA)
        
        self.axesA = self.figA.add_subplot(111, aspect='equal', frame_on = False)
        self.axesA.set_xlim([self.Lim_X1, self.Lim_X2])
        self.axesA.set_ylim([self.Lim_Y1, self.Lim_Y2])  

        self.canvasA.mpl_connect('pick_event', self.on_pick)

        self.cb_gridA = wx.CheckBox(self.AxialTab, -1, "Show Grid", style=wx.ALIGN_RIGHT)
        self.cb_gridA.SetValue(True)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_gridA, self.cb_gridA)

        self.toolbar = NavigationToolbar(self.canvasA)
        
        self.hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox1.Add(self.toolbar, 0, wx.EXPAND)
        flags = wx.ALIGN_RIGHT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        self.hbox1.AddSpacer(20)
        self.hbox1.Add(self.cb_gridA, 0, border=3, flag=flags)
        
        self.vbox1 = wx.BoxSizer(wx.VERTICAL)
        self.vbox1.Add(self.canvasA, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox1.Add(self.hbox1, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        
        self.gridAxial = wx.grid.Grid(self.AxialTab, size=(243,300))
        
        self.gridAxial.EnableEditing(False)
        self.gridAxial.SetDefaultCellAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
        self.gridAxial.CreateGrid(len(self.ELE[:,0]),3)
        self.gridAxial.SetRowLabelSize(3)
        self.gridAxial.SetColLabelValue(0, 'Element ID')
        self.gridAxial.SetColLabelValue(1, 'Nj')
        self.gridAxial.SetColLabelValue(2, 'Nk')

        for i in range(0,len(self.ELE[:,0])):
            self.gridAxial.SetCellValue(i,0,str(i+1))
            self.gridAxial.SetCellValue(i,1,'%.2f'%self.P_element[i][1])
            self.gridAxial.SetCellValue(i,2,'%.2f'%self.P_element[i][4])

        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.vbox1, 0, wx.EXPAND)
        self.hbox2.Add(self.gridAxial, 0, wx.EXPAND)


        self.AxialTab.SetSizer(self.hbox2)

        
        self.draw_figure(self.NOD,'black',True,self.canvasA,self.axesA,self.cb_gridA,False,True)
         
        scale_h = self.get_scale(np.concatenate((self.P_element[:,1],self.P_element[:,4])), self.Dim_X/10)
        
        for i in range(0,len(self.ELE[:,0])):
            [x1,y1,phi] = self.get_element_properties(i,['x1','y1','phi'])

            xs = np.zeros(self.accuracy)
            ys = np.zeros(self.accuracy)
            x_el = np.zeros(self.accuracy)
            y_el = np.zeros(self.accuracy)
            up_down = np.zeros(self.accuracy)
            
            xs_l = self.diagrams[:,5*i]#np.loadtxt('diagrams.txt', usecols=range(5*i,5*i + 1))
            ys_l = self.diagrams[:,5*i + 1]#np.loadtxt('diagrams.txt', usecols=range(5*i + 1,5*i + 2))
            
            for j in range (1,self.accuracy):
                if ys_l[j]*ys_l[j-1] <= 0:
                    up_down[j] = 1
            
            for j in range (0,self.accuracy):
                xs[j] = xs_l[j]
                ys[j] = ys_l[j]*scale_h*-1
                x_el[j] = np.cos(phi)*xs[j] + x1
                y_el[j] = np.sin(phi)*xs[j] + y1
                x_new = np.cos(phi)*xs[j] - np.sin(phi)*ys[j] + x1
                y_new = np.sin(phi)*xs[j] + np.cos(phi)*ys[j] + y1
                xs[j] = x_new
                ys[j] = y_new

            if ys[0] >= 0:
                up = 1
            else:
                up = -1

            verts = [(xs[0],ys[0])]
            codes = [Path.MOVETO]
            xy_el = (x_el[0], y_el[0])
            xy_sta = (xs[0],ys[0])

            for j in range(1,self.accuracy):
                if (up_down[j] != 1) and (j < self.accuracy - 1):
                    verts.append((xs[j],ys[j]))
                    codes.append(Path.LINETO)
                else:
                    verts.append((x_el[j],y_el[j]))
                    verts.append( xy_el )
                    verts.append( xy_sta )
                    codes.append(Path.LINETO)
                    codes.append(Path.LINETO)
                    codes.append(Path.LINETO)
                    path = Path(verts, codes)
                    if up == 1:
                        patch = patches.PathPatch(path, edgecolor = 'blue', facecolor='blue', lw=1, alpha=0.5)
                    else:
                        patch = patches.PathPatch(path, edgecolor = 'red', facecolor='red', lw=1, alpha=0.5)
                        
                    self.axesA.add_patch(patch)
                    verts = [(xs[j],ys[j])]
                    codes = [Path.MOVETO]
                    xy_el = (x_el[j], y_el[j])
                    xy_sta = (xs[j],ys[j])
                    up = -1*up

            self.axesA.text(xs[0],ys[0],'%.2f'%ys_l[0], size='xx-small',color='k')
            self.axesA.text(xs[-1],ys[-1],'%.2f'%ys_l[-1], size='xx-small',color='k')
            ind_max = np.argmax(ys_l)
            ind_min = np.argmin(ys_l)
            if (ys_l[ind_max] != ys_l[0]) and (ys_l[ind_max] != ys_l[-1]):
                self.axesA.text(xs[ind_max],ys[ind_max],'max='+'%.2f'%ys_l[ind_max], size='xx-small',color='k')

            if (ys_l[ind_min] != ys_l[0]) and (ys_l[ind_min] != ys_l[-1]):
                self.axesA.text(xs[ind_min],ys[ind_min],'min='+'%.2f'%ys_l[ind_min], size='xx-small',color='k')
                
        self.canvasA.draw()

    def Shear_forces(self, parent):   
        self.figS = Figure((6.0, 5.0), dpi=self.dpi, tight_layout = True, facecolor = 'white')
        self.canvasS = FigCanvas(self.ShearTab, -1, self.figS)
        
        self.axesS = self.figS.add_subplot(111, aspect='equal', frame_on = False)
        self.axesS.set_xlim([self.Lim_X1, self.Lim_X2])
        self.axesS.set_ylim([self.Lim_Y1, self.Lim_Y2])  

        self.canvasS.mpl_connect('pick_event', self.on_pick)

        self.cb_gridS = wx.CheckBox(self.ShearTab, -1, "Show Grid", style=wx.ALIGN_RIGHT)
        self.cb_gridS.SetValue(True)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_gridS, self.cb_gridS)

        self.toolbar = NavigationToolbar(self.canvasS)
        
        self.hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox1.Add(self.toolbar, 0, wx.EXPAND)
        flags = wx.ALIGN_RIGHT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        self.hbox1.AddSpacer(20)
        self.hbox1.Add(self.cb_gridS, 0, border=3, flag=flags)
        
        self.vbox1 = wx.BoxSizer(wx.VERTICAL)
        self.vbox1.Add(self.canvasS, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox1.Add(self.hbox1, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        
        self.gridShear = wx.grid.Grid(self.ShearTab, size=(243,300))
        
        self.gridShear.EnableEditing(False)
        self.gridShear.SetDefaultCellAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
        self.gridShear.CreateGrid(len(self.ELE[:,0]),3)
        self.gridShear.SetRowLabelSize(3)
        self.gridShear.SetColLabelValue(0, 'Element ID')
        self.gridShear.SetColLabelValue(1, 'Qj')
        self.gridShear.SetColLabelValue(2, 'Qk')

        for i in range(0,len(self.ELE[:,0])):
            self.gridShear.SetCellValue(i,0,str(i+1))
            self.gridShear.SetCellValue(i,1,'%.2f'%self.P_element[i][2])
            self.gridShear.SetCellValue(i,2,'%.2f'%self.P_element[i][5])

        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.vbox1, 0, wx.EXPAND)
        self.hbox2.Add(self.gridShear, 0, wx.EXPAND)


        self.ShearTab.SetSizer(self.hbox2)
        
        self.draw_figure(self.NOD,'black',True,self.canvasS,self.axesS,self.cb_gridS,False,True)
        
        scale_h = self.get_scale(np.concatenate((self.P_element[:,2],self.P_element[:,5])), self.Dim_X/10)
        
        for i in range(0,len(self.ELE[:,0])):
            [x1,y1,phi] = self.get_element_properties(i,['x1','y1','phi'])

            xs = np.zeros(self.accuracy)
            ys = np.zeros(self.accuracy)
            x_el = np.zeros(self.accuracy)
            y_el = np.zeros(self.accuracy)
            up_down = np.zeros(self.accuracy)
            
            xs_l = self.diagrams[:,5*i]#np.loadtxt('diagrams.txt', usecols=range(5*i,5*i + 1))
            ys_l = self.diagrams[:,5*i + 2]#np.loadtxt('diagrams.txt', usecols=range(5*i + 2,5*i + 3))

            for j in range (1,self.accuracy):
                if ys_l[j]*ys_l[j-1] <= 0:
                    up_down[j] = 1
            
            for j in range (0,self.accuracy):
                xs[j] = xs_l[j]
                ys[j] = ys_l[j]*scale_h
                x_el[j] = np.cos(phi)*xs[j] + x1
                y_el[j] = np.sin(phi)*xs[j] + y1
                x_new = np.cos(phi)*xs[j] - np.sin(phi)*ys[j] + x1
                y_new = np.sin(phi)*xs[j] + np.cos(phi)*ys[j] + y1
                xs[j] = x_new
                ys[j] = y_new
            
            if ys[0] >= 0:
                up = 1
            else:
                up = -1
            
            verts = [(xs[0],ys[0])]
            codes = [Path.MOVETO]
            xy_el = (x_el[0], y_el[0])
            xy_sta = (xs[0],ys[0])

            for j in range(1,self.accuracy):
                if (up_down[j] != 1) and (j < self.accuracy - 1):
                    verts.append((xs[j],ys[j]))
                    codes.append(Path.LINETO)
                else:
                    verts.append((x_el[j],y_el[j]))
                    verts.append( xy_el )
                    verts.append( xy_sta )
                    codes.append(Path.LINETO)
                    codes.append(Path.LINETO)
                    codes.append(Path.LINETO)
                    path = Path(verts, codes)
                    if up == 1:
                        patch = patches.PathPatch(path, edgecolor = 'blue', facecolor='blue', lw=1, alpha=0.5)
                    else:
                        patch = patches.PathPatch(path, edgecolor = 'red', facecolor='red', lw=1, alpha=0.5)
                        
                    self.axesS.add_patch(patch)
                    verts = [(xs[j],ys[j])]
                    codes = [Path.MOVETO]
                    xy_el = (x_el[j], y_el[j])
                    xy_sta = (xs[j],ys[j])
                    up = -1*up

            ys_l = ys_l*-1
            self.axesS.text(xs[0],ys[0],'%.2f'%ys_l[0], size='xx-small',color='k')
            self.axesS.text(xs[-1],ys[-1],'%.2f'%ys_l[-1], size='xx-small',color='k')
            ind_max = np.argmax(ys_l)
            ind_min = np.argmin(ys_l)
            if (ys_l[ind_max] != ys_l[0]) and (ys_l[ind_max] != ys_l[-1]):
                self.axesS.text(xs[ind_max],ys[ind_max],'max='+'%.2f'%ys_l[ind_max], size='xx-small',color='k')

            if (ys_l[ind_min] != ys_l[0]) and (ys_l[ind_min] != ys_l[-1]):
                self.axesS.text(xs[ind_min],ys[ind_min],'min='+'%.2f'%ys_l[ind_min], size='xx-small',color='k')
                
        self.canvasS.draw()

    def Bending_moments(self, parent):   
        self.figM = Figure((6.0, 5.0), dpi=self.dpi, tight_layout = True, facecolor = 'white')
        self.canvasM = FigCanvas(self.MomentTab, -1, self.figM)
        
        self.axesM = self.figM.add_subplot(111, aspect='equal', frame_on = False)
        self.axesM.set_xlim([self.Lim_X1, self.Lim_X2])
        self.axesM.set_ylim([self.Lim_Y1, self.Lim_Y2])  

        self.canvasM.mpl_connect('pick_event', self.on_pick)

        self.cb_gridM = wx.CheckBox(self.MomentTab, -1, "Show Grid", style=wx.ALIGN_RIGHT)
        self.cb_gridM.SetValue(True)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_gridM, self.cb_gridM)

        self.toolbar = NavigationToolbar(self.canvasM)
        
        self.hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox1.Add(self.toolbar, 0, wx.EXPAND)
        flags = wx.ALIGN_RIGHT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        self.hbox1.AddSpacer(20)
        self.hbox1.Add(self.cb_gridM, 0, border=3, flag=flags)
        
        self.vbox1 = wx.BoxSizer(wx.VERTICAL)
        self.vbox1.Add(self.canvasM, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox1.Add(self.hbox1, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        
        self.gridMoment = wx.grid.Grid(self.MomentTab, size=(243,300))
        
        self.gridMoment.EnableEditing(False)
        self.gridMoment.SetDefaultCellAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
        self.gridMoment.CreateGrid(len(self.ELE[:,0]),3)
        self.gridMoment.SetRowLabelSize(3)
        self.gridMoment.SetColLabelValue(0, 'Element ID')
        self.gridMoment.SetColLabelValue(1, 'Mj')
        self.gridMoment.SetColLabelValue(2, 'Mk')

        for i in range(0,len(self.ELE[:,0])):
            self.gridMoment.SetCellValue(i,0,str(i+1))
            self.gridMoment.SetCellValue(i,1,'%.2f'%self.P_element[i][3])
            self.gridMoment.SetCellValue(i,2,'%.2f'%self.P_element[i][6])

        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.vbox1, 0, wx.EXPAND)
        self.hbox2.Add(self.gridMoment, 0, wx.EXPAND)

        self.MomentTab.SetSizer(self.hbox2)

        self.draw_figure(self.NOD,'black',True,self.canvasM,self.axesM,self.cb_gridM,False,True)
         
        scale_h = 3*self.get_scale(np.concatenate((self.P_element[:,3],self.P_element[:,6])), self.Dim_X/10)
        
        
        self.Mx = np.zeros([self.accuracy,len(self.ELE[:,0])])
        for i in range(0,len(self.ELE[:,0])):
            [x1,y1,phi] = self.get_element_properties(i,['x1','y1','phi'])

            xs = np.zeros(self.accuracy)
            ys = np.zeros(self.accuracy)
            x_el = np.zeros(self.accuracy)
            y_el = np.zeros(self.accuracy)
            up_down = np.zeros(self.accuracy)
            
            xs_l = self.diagrams[:,5*i]#np.loadtxt('diagrams.txt', usecols=range(5*i,5*i + 1))
            ys_l = self.diagrams[:,5*i + 3]#np.loadtxt('diagrams.txt', usecols=range(5*i + 3,5*i + 4))

            for j in range (1,self.accuracy):
                if ys_l[j]*ys_l[j-1] <= 0:
                    up_down[j] = 1
            
            for j in range (0,self.accuracy):
                xs[j] = xs_l[j]
                ys[j] = ys_l[j]*scale_h*-1
                x_el[j] = np.cos(phi)*xs[j] + x1
                y_el[j] = np.sin(phi)*xs[j] + y1
                x_new = np.cos(phi)*xs[j] - np.sin(phi)*ys[j] + x1
                y_new = np.sin(phi)*xs[j] + np.cos(phi)*ys[j] + y1
                xs[j] = x_new
                ys[j] = y_new
            
            if ys[0] >= 0:
                up = 1
            else:
                up = -1
            
            verts = [(xs[0],ys[0])]
            codes = [Path.MOVETO]
            xy_el = (x_el[0], y_el[0])
            xy_sta = (xs[0],ys[0])
 
            for j in range(1,self.accuracy):
                if (up_down[j] != 1) and (j < self.accuracy - 1):
                    verts.append((xs[j],ys[j]))
                    codes.append(Path.LINETO)
                else:
                    verts.append((x_el[j],y_el[j]))
                    verts.append( xy_el )
                    verts.append( xy_sta )
                    codes.append(Path.LINETO)
                    codes.append(Path.LINETO)
                    codes.append(Path.LINETO)
                    path = Path(verts, codes)
                    if up == 1:
                        patch = patches.PathPatch(path, edgecolor = 'blue', facecolor='blue', lw=1, alpha=0.5)
                    else:
                        patch = patches.PathPatch(path, edgecolor = 'red', facecolor='red', lw=1, alpha=0.5)
                         
                    self.axesM.add_patch(patch)
                    verts = [(xs[j],ys[j])]
                    codes = [Path.MOVETO]
                    xy_el = (x_el[j], y_el[j])
                    xy_sta = (xs[j],ys[j])
                    up = -1*up
            
            self.Mx[:,i] = ys_l[:]
            
            self.axesM.text(xs[0],ys[0],'%.2f'%ys_l[0], size='xx-small',color='k')
            self.axesM.text(xs[-1],ys[-1],'%.2f'%ys_l[-1], size='xx-small',color='k')
            ind_max = np.argmax(ys_l)
            ind_min = np.argmin(ys_l)
            if (ys_l[ind_max] != ys_l[0]) and (ys_l[ind_max] != ys_l[-1]):
                self.axesM.text(xs[ind_max],ys[ind_max],'max='+'%.2f'%ys_l[ind_max], size='xx-small',color='k')

            if (ys_l[ind_min] != ys_l[0]) and (ys_l[ind_min] != ys_l[-1]):
                self.axesM.text(xs[ind_min],ys[ind_min],'min='+'%.2f'%ys_l[ind_min], size='xx-small',color='k')
                
        self.canvasM.draw()

    def Mode_shapes(self, parent):
        self.dm = 3*len(self.NOD[:,0])
        f = open('eigenvectors.txt')
        line = f.readline()
        f.close()
        self.n_modes = len(line.split())
         
        self.figMS = Figure((6.0, 5.0), dpi=self.dpi, tight_layout = True, facecolor = 'white')
        self.canvasMS = FigCanvas(self.ModeTad, -1, self.figMS)
        
        self.axesMS = self.figMS.add_subplot(111, aspect='equal', frame_on = False)
        self.axesMS.set_xlim([self.Lim_X1, self.Lim_X2])
        self.axesMS.set_ylim([self.Lim_Y1, self.Lim_Y2])  

        self.canvasMS.mpl_connect('pick_event', self.on_pick)

        self.cb_gridMS = wx.CheckBox(self.ModeTad, -1, "Show Grid", style=wx.ALIGN_RIGHT)
        self.cb_gridMS.SetValue(True)
        self.Bind(wx.EVT_CHECKBOX, self.on_cb_gridMS, self.cb_gridMS)

        self.toolbar = NavigationToolbar(self.canvasMS)

        self.vbox21 = wx.BoxSizer(wx.VERTICAL)
        self.vbox22 = wx.BoxSizer(wx.VERTICAL)
        
        quote = wx.StaticText(self.ModeTad, -1, "Mode: ")
        self.vbox21.Add(quote, 0, wx.ALIGN_CENTRE|wx.ALL, 5)    
        
        quote = wx.StaticText(self.ModeTad, -1, "Scale factor: ")
        self.vbox22.Add(quote, 0, wx.ALIGN_CENTRE|wx.ALL, 5)

        self.Mode = 1

        self.ys_all = np.zeros((self.n_modes*self.accuracy,len(self.ELE[:,0])))
            
        self.ys_all = self.eigen_diagrams#np.loadtxt('eigen_diagrams.txt', usecols=range(i,i+1))

        self.Scale_factor = 1#int(self.Dim_X/10)

        self.spin = wx.SpinCtrl(self.ModeTad, -1, str(int(self.Mode)), size=(80,-1))
        self.spin.SetRange(1,self.n_modes)
        self.spin.SetValue(1)
        self.Bind(wx.EVT_TEXT, self.EvtMode, self.spin)
        self.vbox21.Add(self.spin, 1, wx.ALIGN_CENTRE|wx.ALL,5)

        text = wx.TextCtrl(self.ModeTad, -1, str(self.Scale_factor), size=(80,-1))
        self.Bind(wx.EVT_TEXT, self.EvtSc_MS, text)
        self.vbox22.Add(text, 1, wx.ALIGN_CENTRE|wx.ALL,5)

        self.gridMS = wx.grid.Grid(self.ModeTad, size=(323,300))

        self.gridMS.EnableEditing(False)
        self.gridMS.SetDefaultCellAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
        self.gridMS.CreateGrid(len(self.NOD[:,0]),4)
        self.gridMS.SetRowLabelSize(3)
        self.gridMS.SetColLabelValue(0, 'Node ID')
        self.gridMS.SetColLabelValue(1, 'dx')
        self.gridMS.SetColLabelValue(2, 'dy')
        self.gridMS.SetColLabelValue(3, 'd'u'\N{greek small letter phi}')

#         for i in range(0,len(self.NOD[:,0])):
                      
        self.update_gridMS()

        self.vbox2 = wx.BoxSizer(wx.VERTICAL)
        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)        
       
        self.hbox2.Add(self.vbox21, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
        self.hbox2.Add(self.vbox22, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
        
        self.vbox2.Add(self.hbox2, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
       
        self.vbox2.AddSpacer(30)
        self.vbox2.Add(self.gridMS, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
        
        self.hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox1.Add(self.toolbar, 0, wx.EXPAND)
        flags = wx.ALIGN_RIGHT | wx.ALL | wx.ALIGN_CENTER_VERTICAL
        self.hbox1.AddSpacer(20)
        self.hbox1.Add(self.cb_gridMS, 0, border=3, flag=flags)
        
        self.vbox1 = wx.BoxSizer(wx.VERTICAL)
        self.vbox1.Add(self.canvasMS, 1, wx.LEFT | wx.TOP | wx.GROW)
        self.vbox1.Add(self.hbox1, 0, flag = wx.ALIGN_LEFT | wx.TOP)
        
        self.hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.hbox2.Add(self.vbox1, 0, border=3)#, flag=flags)
        self.hbox2.Add(self.vbox2, 0, border=3)#, flag=flags)
        
        self.ModeTad.SetSizer(self.hbox2)
        
        self.mode_draw_by_scale()
                
    def EvtSc_MS(self,event):
        self.Scale_factor = float(event.GetString())
        self.mode_draw_by_scale()

    def EvtMode(self,event):
        self.Mode = self.spin.GetValue()
        self.update_gridMS()
        self.mode_draw_by_scale()

    def update_gridMS(self):
        self.PHI_nodal = np.zeros((len(self.NOD[:,0]),3))
        PHI = self.eigenvectors[:,self.Mode - 1]#np.loadtxt('eigenvectors.txt', usecols=range(self.Mode - 1,self.Mode))
        for i in range(0,len(self.NOD[:,0])):
            self.PHI_nodal[i,0] = PHI[(i)*3]
            self.PHI_nodal[i,1] = PHI[(i)*3 + 1]
            self.PHI_nodal[i,2] = PHI[(i)*3 + 2]
            self.gridMS.SetCellValue(i,0,str(i+1))
            self.gridMS.SetCellValue(i,1,'%.2E'%self.PHI_nodal[i][0])
            self.gridMS.SetCellValue(i,2,'%.2E'%self.PHI_nodal[i][1])
            self.gridMS.SetCellValue(i,3,'%.2E'%self.PHI_nodal[i][2])

    def mode_draw_by_scale(self):
        self.NOD_D = deepcopy(self.NOD)
        for i in range(0,len(self.NOD[:,0])):
            self.NOD_D[i][1] = self.NOD_D[i][1] + self.Scale_factor*float(self.PHI_nodal[i][0])
            self.NOD_D[i][2] = self.NOD_D[i][2] + self.Scale_factor*float(self.PHI_nodal[i][1])
        
        self.draw_figure(self.NOD,'grey',True,self.canvasMS,self.axesMS,self.cb_gridMS,True,True)
        self.draw_figure(self.NOD_D,'black',False,self.canvasMS,self.axesMS,self.cb_gridMS,True,False)

        self.ys_l = np.zeros((self.accuracy,len(self.ELE[:,0])))
        self.ys_l = self.ys_all[(self.Mode - 1)*self.accuracy:self.Mode*self.accuracy,:]

        for i in range(0,len(self.ELE[:,0])):
            [phi] = self.get_element_properties(i,['phi'])
            
            x1_ = self.NOD_D[self.ELE[i][1]-1][1]
            y1_ = self.NOD_D[self.ELE[i][1]-1][2]
            
            x2_ = self.NOD_D[self.ELE[i][2]-1][1]
            y2_ = self.NOD_D[self.ELE[i][2]-1][2]
            
            if (x2_ - x1_) != 0:
                theta = np.arctan((y2_ - y1_)/(x2_ - x1_))
            elif ((y2_ - y1_) > 0):
                theta = np.pi/2
            elif ((y2_ - y1_) < 0):
                theta = 3*np.pi/2
            else:
                theta = 0
 
            Ld = ((x2_-x1_)**2 + (y2_-y1_)**2)**0.5
            Ld = Ld*abs(np.cos(phi-theta))
            delta_x = Ld/(self.accuracy-1)

            xs = np.zeros((self.accuracy-1) + 1)
            ys = np.zeros((self.accuracy-1) + 1)
                        
            cy = self.ys_l[0][i]
            
            k = 0
            for j in range (0,self.accuracy):
                xs[k] = (j-1)*delta_x
                ys[k] = (self.ys_l[j][i] - cy)*self.Scale_factor
                x_new = np.cos(phi)*xs[k] - np.sin(phi)*ys[k] + x1_
                y_new = np.sin(phi)*xs[k] + np.cos(phi)*ys[k] + y1_
                xs[k] = x_new
                ys[k] = y_new
                k = k + 1
        
            self.axesMS.plot(xs,ys,'black')
            
        self.canvasMS.draw()
            
            
    def find_nearest(self,array,value):
        idx = (np.abs(array-value)).argmin()
        return idx

    def get_scale(self,array,value):
        max_val = 0
        for i in range(0,len(array)):
            if abs(array[i]) > max_val:
                max_val = abs(array[i])

        if max_val == 0:
            max_val = 1

        scale_h = value/max_val
        return scale_h

    def get_element_properties(self, ind, props):
        x1, y1 = self.NOD[self.ELE[ind][1]-1][1], self.NOD[self.ELE[ind][1]-1][2]
        x2, y2 = self.NOD[self.ELE[ind][2]-1][1], self.NOD[self.ELE[ind][2]-1][2]
         
        dx = x2 - x1
        dy = y2 - y1
        
        if 'L' in props:
            L = (dx**2 + dy**2)**0.5
        
        if 'phi' in props:
            if dx != 0:
                if (dy >= 0) and (dx > 0):
                    phi = np.arctan(dy/dx)
                elif (dy >= 0) and (dx < 0):
                    phi = np.pi - np.arctan(abs(dy/dx))
                elif (dy <= 0) and (dx > 0):
                    phi = 2*np.pi - np.arctan(abs(dy/dx))
                elif (dy <= 0) and (dx < 0):
                    phi = np.pi + np.arctan(abs(dy/dx))       
            else:
                if dy <= 0:
                    phi = 3*np.pi/2
                else:
                    phi = np.pi/2
        
        properties = []
        for i in range(0,len(props)):
            properties.append(eval(props[i]))
            
        return properties

    def draw_figure(self,NOD,colr,clear,canv,ax,cb_grid,plot_bc, plot_bars):
        """ Redraws the figure
        """
        # clear the axes and redraw the plot anew
        #
        if clear == True:
            ax.clear()
            ax.set_xlim([self.Lim_X1, self.Lim_X2])
            ax.set_ylim([self.Lim_Y1, self.Lim_Y2])       
            ax.grid(cb_grid.IsChecked())
            

        
        for i in range(0,len(NOD[:,0])):
            if NOD[i,0] != 0:
                ax.scatter(NOD[i][1], NOD[i][2], color = 'black', s=1)
                dist = float(self.Dim_X)/float(60)
#                 if self.labelN.GetValue() == True:
#                     ax.text(NOD[i][1]+dist,NOD[i][2]+dist, 
#                                    str(int(NOD[i][0])), size='small',color='k')
        
        if plot_bars == True:
            for i in range(0,len(self.ELE[:,0])):
                if (self.ELE[i,1] != 0) and (self.ELE[i,2] != 0):
                    X1, X2 = NOD[self.ELE[i,1]-1,1], NOD[self.ELE[i,2]-1,1]
                    Y1, Y2 = NOD[self.ELE[i,1]-1,2], NOD[self.ELE[i,2]-1,2]
                    ax.plot((X1,X2), (Y1,Y2), color = colr)
                    dist = float(self.Dim_X)/float(60)
                    p = ((X2+X1)/2-2*dist,(Y2+Y1)/2+dist)
#                 if self.labelE.GetValue() == True:
#                     ax.text(p[0],p[1], 
#                                    'E'+str(int(self.ELE[i][0])), size='small',color='k',alpha=0.5)
        if plot_bc == True:    
            for i in range(0,len(self.BC[:,0])):
                p = [NOD[self.BC[i,0]-1,1], NOD[self.BC[i,0]-1,2]]
                if self.BC[i,0] != 0:
                    if (self.BC[i,1] == 1) and (self.BC[i,2] == 1) and (self.BC[i,3] == 1):
                        self.plot_fixed(p,0.03*self.Dim_X,self.BC[i,4],colr,ax)
                    elif (self.BC[i,1] == 1) and (self.BC[i,2] == 1) and (self.BC[i,3] == 0):
                        self.plot_hinge(p,0.03*self.Dim_X,self.BC[i,4],colr,ax)
                    else:
                        if (self.BC[i,1] == 0) and (self.BC[i,2] == 1) and (self.BC[i,3] == 0):
                            self.plot_roller(p,0.03*self.Dim_X,self.BC[i,4],colr,ax)
                        elif (self.BC[i,1] == 1) and (self.BC[i,2] == 0) and (self.BC[i,3] == 0):
                            self.plot_roller(p,0.03*self.Dim_X,90,colr,ax)


        canv.draw()

    def plot_fixed(self, position, scale, angle, colr, ax):
        pi = 3.14159265359
        
        verts = [
             (-1.25,  0.0), (-1.0,  0.0), (-0.5, -0.3), (-1.0,  0.0), (-0.5,  0.0), 
             (0.0  , -0.3), (-0.5,  0.0), (0.0 ,  0.0), (0.5 , -0.3), (0.0 ,  0.0), 
             (0.5  ,  0.0), (1.0 , -0.3), (0.5 ,  0.0), (1.0 ,  0.0), (1.5 , -0.3), 
             (1.0  ,  0.0), (1.5 ,  0.0)]
         
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO,
                 Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO,
                 Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO,
                 Path.LINETO, Path.LINETO]
        
        angle = pi*angle/180
        vertsT = verts[:]
        for i in range(0,len(verts)):
            vertsT[i] = (verts[i][0]*np.cos(angle)-verts[i][1]*np.sin(angle), verts[i][0]*np.sin(angle)+verts[i][1]*np.cos(angle))

        for i in range(0,len(verts)):
            vertsT[i] = (scale*vertsT[i][0],scale*vertsT[i][1])
    
        for i in range(0,len(verts)):
            vertsT[i] = (vertsT[i][0] + position[0],vertsT[i][1] + position[1])
     
        pathT = Path(vertsT, codes)

        hingeT = patches.PathPatch(pathT, edgecolor = colr,  facecolor='none', lw=1)
        
        ax.add_patch(hingeT)

    def plot_hinge(self, position, scale, angle, colr, ax):
        pi = 3.14159265359
        
        verts = [(-1.25, -0.75), (-1.0, -0.75), (-0.5 , -1.05), (-1.0, -0.75), (-0.5, -0.75), 
                 (0.0  , -1.05), (-0.5, -0.75), (0.0  , -0.75), (0.5 , -1.05), (0.0, -0.75 ), 
                 (0.5  , -0.75), (1.0 , -1.05), (0.5  , -0.75), (1.0 , -0.75), (1.5, -1.05 ), 
                 (1.0  , -0.75), (1.5 , -0.75), (-0.75, -0.75), (0.0 , 0.0  ), (0.75, -0.75)]
        
        
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO,
                 Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO,
                 Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO,
                 Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO]
        
        angle = pi*angle/180
        vertsT = verts[:]
        for i in range(0,len(verts)):
            vertsT[i] = (verts[i][0]*np.cos(angle)-verts[i][1]*np.sin(angle), verts[i][0]*np.sin(angle)+verts[i][1]*np.cos(angle))

        for i in range(0,len(verts)):
            vertsT[i] = (scale*vertsT[i][0],scale*vertsT[i][1])
    
        for i in range(0,len(verts)):
            vertsT[i] = (vertsT[i][0] + position[0],vertsT[i][1] + position[1])
     
        pathT = Path(vertsT, codes)

        hingeT = patches.PathPatch(pathT, edgecolor = colr, facecolor='none', lw=1)
        
        ax.add_patch(hingeT)

    def plot_roller(self, position, scale, angle, colr, ax):
        pi = 3.14159265359
        
        verts = [(-1.25, -1.  ), (-1.0 , -1.  ), (-0.5 , -1.30), (-1.0 , -1.  ), (-0.5, -1. ), 
                 (0.0  , -1.30), (-0.5 , -1.  ), (0.0  , -1.  ), (0.5  , -1.30), (0.0, -1.  ), 
                 (0.5  , -1.  ), (1.0  , -1.30), (0.5  , -1.  ), (1.0  , -1.  ), (1.5, -1.30), 
                 (1.0  , -1.  ), (1.5  , -1.  ), (-0.75, -0.75), (-0.75, -0.75), (0.0, 0.0  ), 
                 (0.75 , -0.75), (-0.75, -0.75)]
    
        codes = [Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, 
                 Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO,
                 Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO,
                 Path.LINETO, Path.LINETO, Path.MOVETO, Path.LINETO, Path.LINETO,
                 Path.LINETO, Path.LINETO,]
        
        angle = pi*angle/180
        vertsT = verts[:]
        for i in range(0,len(verts)):
            vertsT[i] = (verts[i][0]*np.cos(angle)-verts[i][1]*np.sin(angle), verts[i][0]*np.sin(angle)+verts[i][1]*np.cos(angle))

        for i in range(0,len(verts)):
            vertsT[i] = (scale*vertsT[i][0],scale*vertsT[i][1])
    
        for i in range(0,len(verts)):
            vertsT[i] = (vertsT[i][0] + position[0],vertsT[i][1] + position[1])
     
        pathT = Path(vertsT, codes)

        hingeT = patches.PathPatch(pathT, edgecolor = colr, facecolor='none', lw=1)
        
        ax.add_patch(hingeT)


    def on_pick(self, event):
        box_points = event.artist.get_bbox().get_points()
        msg = "You've clicked on a bar with coords:\n %s" % box_points
        
        dlg = wx.MessageDialog(
            self, 
            msg, 
            "Click!",
            wx.OK | wx.ICON_INFORMATION)

        dlg.ShowModal() 
        dlg.Destroy()        

    def on_cb_gridD(self, event):
        self.axesD.grid(self.cb_gridD.IsChecked())
        self.canvasD.draw()
        
    def on_cb_gridR(self, event):
        self.axesR.grid(self.cb_gridR.IsChecked())
        self.canvasR.draw()

    def on_cb_gridA(self, event):
        self.axesA.grid(self.cb_gridA.IsChecked())
        self.canvasA.draw()
        
    def on_cb_gridS(self, event):
        self.axesS.grid(self.cb_gridS.IsChecked())
        self.canvasS.draw()
        
    def on_cb_gridM(self, event):
        self.axesM.grid(self.cb_gridM.IsChecked())
        self.canvasM.draw()
        
    def on_cb_gridMS(self, event):
        self.axesMS.grid(self.cb_gridMS.IsChecked())
        self.canvasMS.draw()