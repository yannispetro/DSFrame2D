import wx
import wx.grid
import numpy as np

class DimDialog(wx.Dialog):
    def __init__(self,x1,y1,x2,y2,GS):
        title = 'Dimensions'
        
        self.Lim_X1 = x1
        self.Lim_Y1 = y1
        self.Lim_X2 = x2
        self.Lim_Y2 = y2
        
        self.Grid_spacing = GS
        
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)

        sizer = wx.BoxSizer(wx.VERTICAL)

        label = wx.StaticText(self, -1, "Set Plot Limits")
        label.SetHelpText("Assign the approximate dimensions of the problem. These dimensions should be larger that the dimensions of the frame in order to display the plots appropriately")
        sizer.Add(label, 0, wx.ALIGN_CENTRE|wx.ALL, 5)
        
        gridSizer = wx.GridBagSizer(5, 5)

        gridSizer.Add(wx.StaticText(self, -1, "Lower"), (0,1))
        gridSizer.Add(wx.StaticText(self, -1, "Upper"), (0,2))
        gridSizer.Add(wx.StaticText(self, -1, "X axis: "), (1,0))

        text11 = wx.TextCtrl(self, -1, str(self.Lim_X1), size=(80,-1))
        text11.SetHelpText("This is the lower limit of the X axis.")
        self.Bind(wx.EVT_TEXT, self.EvtX1, text11)
        gridSizer.Add(text11, (1,1))

        text12 = wx.TextCtrl(self, -1, str(self.Lim_X2), size=(80,-1))
        text12.SetHelpText("This is the upper limit of the X axis.")
        self.Bind(wx.EVT_TEXT, self.EvtX2, text12)
        gridSizer.Add(text12, (1,2))

        gridSizer.Add(wx.StaticText(self, -1, "Y axis: "), (2,0))

        text21 = wx.TextCtrl(self, -1, str(self.Lim_Y1), size=(80,-1))
        text21.SetHelpText("This is the lower limit of the Y axis.")
        self.Bind(wx.EVT_TEXT, self.EvtY1, text21)
        gridSizer.Add(text21, (2,1))

        text22 = wx.TextCtrl(self, -1, str(self.Lim_Y2), size=(80,-1))
        text22.SetHelpText("This is the upper limit of the Y axis.")
        self.Bind(wx.EVT_TEXT, self.EvtY2, text22)
        gridSizer.Add(text22, (2,2))

        sizer.Add(gridSizer, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        hsizer = wx.BoxSizer(wx.HORIZONTAL)
        
        hsizer.Add(wx.StaticText(self, -1, "Grid Spacing: "), 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5) 
        textGS = wx.TextCtrl(self, -1, str(self.Grid_spacing), size=(80,-1))
        self.Bind(wx.EVT_TEXT, self.EvtGS, textGS)
        hsizer.Add(textGS, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5) 

        sizer.Add(hsizer, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def EvtX1(self, event):
        self.Lim_X1 = float(event.GetString())
 
    def EvtY1(self, event):
        self.Lim_Y1 = float(event.GetString())
        
    def EvtX2(self, event):
        self.Lim_X2 = float(event.GetString())
        
    def EvtY2(self, event):
        self.Lim_Y2 = float(event.GetString())
        
    def EvtGS(self, event):
        self.Grid_spacing = float(event.GetString())
        
class NodalForce(wx.Dialog):
    def __init__(self,FOR,PATH):
        title = 'Nodal Forces'
        
        self.FOR = FOR
        self.PATH = PATH
        
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)


        sizer = wx.BoxSizer(wx.VERTICAL)

        self.gridForces1 = wx.grid.Grid(self, size=(350,150))
        
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self.FGridUpdate)
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.FGridRight)
        self.gridForces1.Bind(wx.grid.EVT_GRID_LABEL_RIGHT_CLICK, self.FGridRight)
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.ActiveRowF)
        self.gridForces1.Bind(wx.grid.EVT_GRID_LABEL_LEFT_CLICK, self.ActiveRowF)            
        
        if len(self.FOR) != 0:
            self.gridForces1.CreateGrid(len(self.FOR),4)
        else:
            self.gridForces1.CreateGrid(0,4)
                                        
        self.gridForces1.SetDefaultCellAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
        self.gridForces1.SetRowLabelSize(3)
        self.gridForces1.SetColLabelValue(0, 'Node ID')
        self.gridForces1.SetColLabelValue(1, 'X force')
        self.gridForces1.SetColLabelValue(2, 'Y force')
        self.gridForces1.SetColLabelValue(3, 'Moment')
#         self.gridForces1.SetColSize(1, 80)
#         self.gridForces1.SetColSize(2, 120)
        
        if len(self.FOR) != 0:
            for i in range(0,len(self.FOR)):
                    if self.FOR[i][0] != 0:
                        self.gridForces1.SetCellValue(i,0,str(int(self.FOR[i][0])))
                    for j in range(1,len(self.FOR[0])):
                        if self.FOR[i][j] != 0:
                            self.gridForces1.SetCellValue(i,j,str(self.FOR[i][j]))

        sizer.Add(self.gridForces1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.buttonF1 = wx.Button(self, label="Add")
        self.Bind(wx.EVT_BUTTON, self.OnClickF, self.buttonF1)
        
        self.buttonDelF1 = wx.Button(self, label="Delete")
        self.Bind(wx.EVT_BUTTON, self.ClickDelFor, self.buttonDelF1)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.buttonF1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        hbox.Add(self.buttonDelF1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)


        Img = wx.StaticBitmap(self, bitmap=wx.Image(self.PATH + '/icons/pos_nf.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap())

        sizer.Add(hbox, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        sizer.Add(Img, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def ActiveRowF(self,event):
        self.active_rowF = event.GetRow()
        event.Skip()

    def ClickDelFor(self,event):
        self.For_to_del = self.active_rowF
        self.OnDelFor(0)

    def FGridUpdate(self, event):
        R = event.GetRow()
        C = event.GetCol()
        Value = float(self.gridForces1.GetCellValue(R, C))
        self.FOR[R][C] = Value
            
    def FGridRight(self,event):
        self.For_to_del = event.GetRow()
        FRMenu = wx.Menu()
        item = wx.MenuItem(FRMenu, wx.NewId(), "Delete Force")
        FRMenu.AppendItem(item)
        self.Bind(wx.EVT_MENU, self.OnDelFor, item)
        self.PopupMenu(FRMenu)
        
    def OnClickF(self,event):
        self.gridForces1.AppendRows(1, updateLabels = True)
        self.FOR.append([0 ,0.,0., 0.])
        
    def OnDelFor(self,event):
        self.gridForces1.DeleteRows(self.For_to_del, 1, updateLabels=True)
        self.FOR = np.delete(self.FOR, self.For_to_del, 0)
        for i in range(0,len(self.FOR[:,0])):
            if self.FOR[i][0] > self.For_to_del + 1:
                self.FOR[i][0] = self.FOR[i][0] - 1

class ElementForce(wx.Dialog):
    def __init__(self,FOR,PATH):
        title = 'Element Forces'
        
        self.FOR = FOR
        self.PATH = PATH
        
        
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)


        sizer = wx.BoxSizer(wx.VERTICAL)

        self.gridForces1 = wx.grid.Grid(self, size=(490,150))
        
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self.FGridUpdate)
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.FGridRight)
        self.gridForces1.Bind(wx.grid.EVT_GRID_LABEL_RIGHT_CLICK, self.FGridRight) 
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.ActiveRowF)
        self.gridForces1.Bind(wx.grid.EVT_GRID_LABEL_LEFT_CLICK, self.ActiveRowF)                
        
        if len(self.FOR) != 0:
            self.gridForces1.CreateGrid(len(self.FOR),5)
        else:
            self.gridForces1.CreateGrid(0,5)
        
        self.gridForces1.SetDefaultCellAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
        self.gridForces1.SetRowLabelSize(3)
        self.gridForces1.SetColLabelValue(0, 'Element ID')
        self.gridForces1.SetColLabelValue(1, 'Axial force')
        self.gridForces1.SetColLabelValue(2, 'Shear force')
        self.gridForces1.SetColLabelValue(3, 'Moment')
        self.gridForces1.SetColLabelValue(4, 'Position')
        for i in range(0,3):
            self.gridForces1.SetColSize(i, 100)
#         self.gridForces1.SetColSize(1, 80)
#         self.gridForces1.SetColSize(2, 120)
        
        if len(self.FOR) != 0:
            for i in range(0,len(self.FOR)):
                    if self.FOR[i][0] != 0:
                        self.gridForces1.SetCellValue(i,0,str(int(self.FOR[i][0])))
                    for j in range(1,len(self.FOR[0][:])):
                        if self.FOR[i][j] != 0:
                            self.gridForces1.SetCellValue(i,j,str(self.FOR[i][j]))

        sizer.Add(self.gridForces1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        
        self.buttonF1 = wx.Button(self, label="Add")
        self.Bind(wx.EVT_BUTTON, self.OnClickF, self.buttonF1)
        
        self.buttonDelF1 = wx.Button(self, label="Delete")
        self.Bind(wx.EVT_BUTTON, self.ClickDelFor, self.buttonDelF1)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.buttonF1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        hbox.Add(self.buttonDelF1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)


        Img = wx.StaticBitmap(self, bitmap=wx.Image(self.PATH + '/icons/pos_ef.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap())

        quote = wx.StaticText(self, label="The variable Position, is the position of the force relative to the length of the element: x_local / L which lies in [0,1]")
        font = wx.Font(9, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
        quote.SetFont(font)
        quote.Wrap(250)
        
        line1 = wx.StaticLine(self, -1, size=(2,-1), style=wx.LI_VERTICAL)
        
        sizer.Add(hbox, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(Img, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        hbox.Add(line1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        hbox.Add(quote, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        sizer.Add(hbox, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def ActiveRowF(self,event):
        self.active_rowF = event.GetRow()
        event.Skip()

    def ClickDelFor(self,event):
        self.For_to_del = self.active_rowF
        self.OnDelFor(0)

    def FGridUpdate(self, event):
        R = event.GetRow()
        C = event.GetCol()
        Value = float(self.gridForces1.GetCellValue(R, C))
        self.FOR[R][C] = Value
            
    def FGridRight(self,event):
        self.For_to_del = event.GetRow()
        FRMenu = wx.Menu()
        item = wx.MenuItem(FRMenu, wx.NewId(), "Delete Load")
        FRMenu.AppendItem(item)
        self.Bind(wx.EVT_MENU, self.OnDelFor, item)
        self.PopupMenu(FRMenu)
        
    def OnClickF(self,event):
        self.gridForces1.AppendRows(1, updateLabels = True)
        self.FOR.append([0 ,0.,0., 0., 0.])
        
    def OnDelFor(self,event):
        self.gridForces1.DeleteRows(self.For_to_del, 1, updateLabels=True)
        self.FOR = np.delete(self.FOR, self.For_to_del, 0)
        for i in range(0,len(self.FOR[:,0])):
            if self.FOR[i][0] > self.For_to_del + 1:
                self.FOR[i][0] = self.FOR[i][0] - 1
        
class DistributedLoad(wx.Dialog):
    def __init__(self,FOR,PATH):
        title = 'Distributed Loads'
        
        self.FOR = FOR
        self.PATH = PATH
        
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)


        sizer = wx.BoxSizer(wx.VERTICAL)

        self.gridForces1 = wx.grid.Grid(self, size=(490,150))
        
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self.FGridUpdate)
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.FGridRight)
        self.gridForces1.Bind(wx.grid.EVT_GRID_LABEL_RIGHT_CLICK, self.FGridRight)   
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.ActiveRowF)
        self.gridForces1.Bind(wx.grid.EVT_GRID_LABEL_LEFT_CLICK, self.ActiveRowF)              
        
        if len(self.FOR) != 0:
            self.gridForces1.CreateGrid(len(self.FOR),4)
        else:
            self.gridForces1.CreateGrid(0,4)
            
        self.gridForces1.SetDefaultCellAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
        self.gridForces1.SetRowLabelSize(3)
        self.gridForces1.SetColLabelValue(0, 'Element ID')
        self.gridForces1.SetColLabelValue(1, 'Axial load')
        self.gridForces1.SetColLabelValue(2, 'Shear load')
        self.gridForces1.SetColLabelValue(3, 'Moment')
#         self.gridForces1.SetColLabelValue(4, 'Position')
        for i in range(0,3):
            self.gridForces1.SetColSize(i, 100)
#         self.gridForces1.SetColSize(2, 120)
        
        if len(self.FOR) != 0:
            for i in range(0,len(self.FOR)):
                    if self.FOR[i][0] != 0:
                        self.gridForces1.SetCellValue(i,0,str(int(self.FOR[i][0])))
                    for j in range(1,len(self.FOR[0][:])-1):
                        if self.FOR[i][j] != 0:
                            self.gridForces1.SetCellValue(i,j,str(self.FOR[i][j]))

        sizer.Add(self.gridForces1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.buttonF1 = wx.Button(self, label="Add")
        self.Bind(wx.EVT_BUTTON, self.OnClickF, self.buttonF1)
        
        self.buttonDelF1 = wx.Button(self, label="Delete")
        self.Bind(wx.EVT_BUTTON, self.ClickDelFor, self.buttonDelF1)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.buttonF1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        hbox.Add(self.buttonDelF1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        Img = wx.StaticBitmap(self, bitmap=wx.Image(self.PATH + '/icons/pos_ef.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap())

        sizer.Add(hbox, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        sizer.Add(Img, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def ActiveRowF(self,event):
        self.active_rowF = event.GetRow()
        event.Skip()

    def ClickDelFor(self,event):
        self.For_to_del = self.active_rowF
        self.OnDelFor(0)

    def FGridUpdate(self, event):
        R = event.GetRow()
        C = event.GetCol()
        Value = float(self.gridForces1.GetCellValue(R, C))
        self.FOR[R][C] = Value
            
    def FGridRight(self,event):
        self.For_to_del = event.GetRow()
        FRMenu = wx.Menu()
        item = wx.MenuItem(FRMenu, wx.NewId(), "Delete Load")
        FRMenu.AppendItem(item)
        self.Bind(wx.EVT_MENU, self.OnDelFor, item)
        self.PopupMenu(FRMenu)
        
    def OnClickF(self,event):
        self.gridForces1.AppendRows(1, updateLabels = True)
        self.FOR.append([0 ,0.,0., 0., 1.])
        
    def OnDelFor(self,event):
        self.gridForces1.DeleteRows(self.For_to_del, 1, updateLabels=True)
        self.FOR = np.delete(self.FOR, self.For_to_del, 0)
        for i in range(0,len(self.FOR[:][0])):
            if self.FOR[i][0] > self.For_to_del + 1:
                self.FOR[i][0] = self.FOR[i][0] - 1        
                
                
class TemperatureLoad(wx.Dialog):
    def __init__(self,FOR,PATH):
        title = 'Temperature Loads'
        
        self.FOR = FOR
        self.PATH = PATH
        
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)


        sizer = wx.BoxSizer(wx.VERTICAL)

        self.gridForces1 = wx.grid.Grid(self, size=(370,150))
        
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self.FGridUpdate)
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.FGridRight)
        self.gridForces1.Bind(wx.grid.EVT_GRID_LABEL_RIGHT_CLICK, self.FGridRight) 
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.ActiveRowF)
        self.gridForces1.Bind(wx.grid.EVT_GRID_LABEL_LEFT_CLICK, self.ActiveRowF)                
        
        if len(self.FOR) != 0:
            self.gridForces1.CreateGrid(len(self.FOR),4)
        else:
            self.gridForces1.CreateGrid(0,4)
            
        self.gridForces1.SetDefaultCellAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
        self.gridForces1.SetRowLabelSize(3)
        self.gridForces1.SetColLabelValue(0, 'Element ID')
        self.gridForces1.SetColLabelValue(1, 'T_in ('u'\N{DEGREE SIGN}C)')
        self.gridForces1.SetColLabelValue(2, 'T_out ('u'\N{DEGREE SIGN}C)')
        self.gridForces1.SetColLabelValue(3, 'T_0 ('u'\N{DEGREE SIGN}C)')
        for i in range(0,1):
            self.gridForces1.SetColSize(i, 100)
#         self.gridForces1.SetColSize(2, 120)
        
        if len(self.FOR) != 0:
            for i in range(0,len(self.FOR)):
                    if self.FOR[i][0] != 0:
                        self.gridForces1.SetCellValue(i,0,str(int(self.FOR[i][0])))
                    for j in range(1,4):
                        if self.FOR[i][j] != 0:
                            self.gridForces1.SetCellValue(i,j,str(self.FOR[i][j]))

        sizer.Add(self.gridForces1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.buttonF1 = wx.Button(self, label="Add")
        self.Bind(wx.EVT_BUTTON, self.OnClickF, self.buttonF1)
        
        self.buttonDelF1 = wx.Button(self, label="Delete")
        self.Bind(wx.EVT_BUTTON, self.ClickDelFor, self.buttonDelF1)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.buttonF1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        hbox.Add(self.buttonDelF1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)


        Img = wx.StaticBitmap(self, bitmap=wx.Image(self.PATH + '/icons/pos_tl.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap())

        sizer.Add(hbox, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)


        sizer.Add(Img, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def ActiveRowF(self,event):
        self.active_rowF = event.GetRow()
        event.Skip()

    def ClickDelFor(self,event):
        self.For_to_del = self.active_rowF
        self.OnDelFor(0)

    def FGridUpdate(self, event):
        R = event.GetRow()
        C = event.GetCol()
        Value = float(self.gridForces1.GetCellValue(R, C))
        self.FOR[R][C] = Value
            
    def FGridRight(self,event):
        self.For_to_del = event.GetRow()
        FRMenu = wx.Menu()
        item = wx.MenuItem(FRMenu, wx.NewId(), "Delete Load")
        FRMenu.AppendItem(item)
        self.Bind(wx.EVT_MENU, self.OnDelFor, item)
        self.PopupMenu(FRMenu)
        
    def OnClickF(self,event):
        self.gridForces1.AppendRows(1, updateLabels = True)
        self.FOR = np.append(self.FOR, [[0 ,0.,0., 0.]], axis=0)
        
    def OnDelFor(self,event):
        self.gridForces1.DeleteRows(self.For_to_del, 1, updateLabels=True)
        self.FOR = np.delete(self.FOR, self.For_to_del, 0)
        for i in range(0,len(self.FOR[:,0])):
            if self.FOR[i][0] > self.For_to_del + 1:
                self.FOR[i][0] = self.FOR[i][0] - 1    
                
            
class SupportDisplacement(wx.Dialog):
    def __init__(self,FOR,PATH):
        title = 'Support Displacements'
        
        self.FOR = FOR
        self.PATH = PATH
        
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)


        sizer = wx.BoxSizer(wx.VERTICAL)

        self.gridForces1 = wx.grid.Grid(self, size=(350,150))
        
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self.FGridUpdate)
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.FGridRight)
        self.gridForces1.Bind(wx.grid.EVT_GRID_LABEL_RIGHT_CLICK, self.FGridRight)
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.ActiveRowF)
        self.gridForces1.Bind(wx.grid.EVT_GRID_LABEL_LEFT_CLICK, self.ActiveRowF)                 
        
        if len(self.FOR) != 0:
            self.gridForces1.CreateGrid(len(self.FOR),4)
        else:
            self.gridForces1.CreateGrid(0,4)
            
        self.gridForces1.SetDefaultCellAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
        self.gridForces1.SetRowLabelSize(3)
        self.gridForces1.SetColLabelValue(0, 'Node ID')
        self.gridForces1.SetColLabelValue(1, u'\N{greek small letter delta}x')
        self.gridForces1.SetColLabelValue(2, u'\N{greek small letter delta}y')
        self.gridForces1.SetColLabelValue(3, u'\N{greek small letter delta}\N{greek small letter phi}')

        for i in range(0,1):
            self.gridForces1.SetColSize(i, 100)
#         self.gridForces1.SetColSize(2, 120)
        
        if len(self.FOR) != 0:
            for i in range(0,len(self.FOR)):
                    if self.FOR[i][0] != 0:
                        self.gridForces1.SetCellValue(i,0,str(int(self.FOR[i][0])))
                    for j in range(1,4):
                        if self.FOR[i][j] != 0:
                            self.gridForces1.SetCellValue(i,j,str(self.FOR[i][j]))

        sizer.Add(self.gridForces1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.buttonF1 = wx.Button(self, label="Add")
        self.Bind(wx.EVT_BUTTON, self.OnClickF, self.buttonF1)
        
        self.buttonDelF1 = wx.Button(self, label="Delete")
        self.Bind(wx.EVT_BUTTON, self.ClickDelFor, self.buttonDelF1)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.buttonF1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        hbox.Add(self.buttonDelF1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)


        Img = wx.StaticBitmap(self, bitmap=wx.Image(self.PATH + '/icons/pos_nf.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap())

        sizer.Add(hbox, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        sizer.Add(Img, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def ActiveRowF(self,event):
        self.active_rowM = event.GetRow()
        event.Skip()

    def ClickDelFor(self,event):
        self.For_to_del = self.active_rowM
        self.OnDelFor(0)

    def FGridUpdate(self, event):
        R = event.GetRow()
        C = event.GetCol()
        Value = float(self.gridForces1.GetCellValue(R, C))
        self.FOR[R][C] = Value
            
    def FGridRight(self,event):
        self.For_to_del = event.GetRow()
        FRMenu = wx.Menu()
        item = wx.MenuItem(FRMenu, wx.NewId(), "Delete Load")
        FRMenu.AppendItem(item)
        self.Bind(wx.EVT_MENU, self.OnDelFor, item)
        self.PopupMenu(FRMenu)
        
    def OnClickF(self,event):
        self.gridForces1.AppendRows(1, updateLabels = True)
        self.FOR = np.append(self.FOR, [[0 ,0.,0., 0.]], axis=0)
        
    def OnDelFor(self,event):
        self.gridForces1.DeleteRows(self.For_to_del, 1, updateLabels=True)
        self.FOR = np.delete(self.FOR, self.For_to_del, 0)
        for i in range(0,len(self.FOR[:,0])):
            if self.FOR[i][0] > self.For_to_del + 1:
                self.FOR[i][0] = self.FOR[i][0] - 1        
                
                
class StructuralDefect(wx.Dialog):
    def __init__(self,FOR,PATH):
        title = 'Structural Defects'
        
        self.FOR = FOR
        self.PATH = PATH
        
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)


        sizer = wx.BoxSizer(wx.VERTICAL)

        self.gridForces1 = wx.grid.Grid(self, size=(490,150))
        
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self.FGridUpdate)
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.FGridRight)
        self.gridForces1.Bind(wx.grid.EVT_GRID_LABEL_RIGHT_CLICK, self.FGridRight)
        self.gridForces1.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.ActiveRowF)
        self.gridForces1.Bind(wx.grid.EVT_GRID_LABEL_LEFT_CLICK, self.ActiveRowF)                 
        
        if len(self.FOR) != 0:
            self.gridForces1.CreateGrid(len(self.FOR[:,0]),len(self.FOR[0,:]))
        else:
            self.gridForces1.CreateGrid(0,5)
            
        self.gridForces1.SetDefaultCellAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
        self.gridForces1.SetRowLabelSize(3)
        self.gridForces1.SetColLabelValue(0, 'Element ID')
        self.gridForces1.SetColLabelValue(1, 'x')
        self.gridForces1.SetColLabelValue(2, 'y')
        self.gridForces1.SetColLabelValue(3, 'phi')
        self.gridForces1.SetColLabelValue(4, 'Position')
        for i in range(0,3):
            self.gridForces1.SetColSize(i, 100)
#         self.gridForces1.SetColSize(1, 80)
#         self.gridForces1.SetColSize(2, 120)
        
        if len(self.FOR) != 0:
            for i in range(0,len(self.FOR[:,0])):
                    if self.FOR[i,0] != 0:
                        self.gridForces1.SetCellValue(i,0,str(int(self.FOR[i,0])))
                    for j in range(1,len(self.FOR[0,:])):
                        if self.FOR[i,j] != 0:
                            self.gridForces1.SetCellValue(i,j,str(self.FOR[i,j]))

        sizer.Add(self.gridForces1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.buttonF1 = wx.Button(self, label="Add")
        self.Bind(wx.EVT_BUTTON, self.OnClickF, self.buttonF1)
        
        self.buttonDelF1 = wx.Button(self, label="Delete")
        self.Bind(wx.EVT_BUTTON, self.ClickDelFor, self.buttonDelF1)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.buttonF1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)
        hbox.Add(self.buttonDelF1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)


        sizer.Add(hbox, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def ActiveRowF(self,event):
        self.active_rowM = event.GetRow()
        event.Skip()

    def ClickDelFor(self,event):
        self.For_to_del = self.active_rowM
        self.OnDelFor(0)

    def FGridUpdate(self, event):
        R = event.GetRow()
        C = event.GetCol()
        Value = float(self.gridForces1.GetCellValue(R, C))
        self.FOR[R][C] = Value
            
    def FGridRight(self,event):
        self.For_to_del = event.GetRow()
        FRMenu = wx.Menu()
        item = wx.MenuItem(FRMenu, wx.NewId(), "Delete Defect")
        FRMenu.AppendItem(item)
        self.Bind(wx.EVT_MENU, self.OnDelFor, item)
        self.PopupMenu(FRMenu)
        
    def OnClickF(self,event):
        self.gridForces1.AppendRows(1, updateLabels = True)
        self.FOR = np.append(self.FOR, [[0 ,0.,0., 0., 0.]], axis=0)
        
    def OnDelFor(self,event):
        self.gridForces1.DeleteRows(self.For_to_del, 1, updateLabels=True)
        self.FOR = np.delete(self.FOR, self.For_to_del, 0)
        for i in range(0,len(self.FOR[:,0])):
            if self.FOR[i][0] > self.For_to_del + 1:
                self.FOR[i][0] = self.FOR[i][0] - 1 
                
                
class ViewOptions(wx.Dialog):
    def __init__(self,sF,lF,lG):
        title = 'View Options'
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)

        sizer = wx.BoxSizer(wx.VERTICAL)
        self.gbs = wx.GridBagSizer(5, 5)

        hbox1 = wx.BoxSizer(wx.HORIZONTAL)
        self.quoteN = wx.StaticText(self, label="   Nodes:              ")
        self.labelN = wx.CheckBox(self, -1, "Labels", style=wx.ALIGN_RIGHT)
        self.labelN.SetValue(lG[0])
        self.Bind(wx.EVT_CHECKBOX, self.on_labelF, self.labelN)
        
        hbox1.Add(self.quoteN, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        hbox1.Add(self.labelN, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        sizer.Add(hbox1, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        
        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        
        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        self.quoteE = wx.StaticText(self, label="   Elements:          ")
        self.labelE = wx.CheckBox(self, -1, "Labels", style=wx.ALIGN_RIGHT)
        self.labelE.SetValue(lG[1])
        self.Bind(wx.EVT_CHECKBOX, self.on_labelF, self.labelE)
        
        
        hbox2.Add(self.quoteE, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        hbox2.Add(self.labelE, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        sizer.Add(hbox2, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        

        self.quoteF1 = wx.StaticText(self, label="   Nodal forces: ")
        self.quoteF2 = wx.StaticText(self, label="   Element forces: ")
        self.quoteF3 = wx.StaticText(self, label="   Distributed loads: ")
        self.quoteF4 = wx.StaticText(self, label="   Temperature loads: ")
        self.quoteF5 = wx.StaticText(self, label="   Support displacements: ")
        self.quoteF6 = wx.StaticText(self, label="   Structural defects: ")    
          
        self.showF1 = wx.CheckBox(self, -1, "Show", style=wx.ALIGN_RIGHT)
        self.showF1.SetValue(sF[0])
        self.Bind(wx.EVT_CHECKBOX, self.on_labelF, self.showF1) 
          
        self.showF2 = wx.CheckBox(self, -1, "Show", style=wx.ALIGN_RIGHT)
        self.showF2.SetValue(sF[1])
        self.Bind(wx.EVT_CHECKBOX, self.on_labelF, self.showF2) 
          
        self.showF3 = wx.CheckBox(self, -1, "Show", style=wx.ALIGN_RIGHT)
        self.showF3.SetValue(sF[2])
        self.Bind(wx.EVT_CHECKBOX, self.on_labelF, self.showF3) 
          
        self.showF4 = wx.CheckBox(self, -1, "Show", style=wx.ALIGN_RIGHT)
        self.showF4.SetValue(sF[3])
        self.Bind(wx.EVT_CHECKBOX, self.on_labelF, self.showF4) 
          
        self.showF5 = wx.CheckBox(self, -1, "Show", style=wx.ALIGN_RIGHT)
        self.showF5.SetValue(sF[4])
        self.Bind(wx.EVT_CHECKBOX, self.on_labelF, self.showF5) 
          
        self.showF6 = wx.CheckBox(self, -1, "Show", style=wx.ALIGN_RIGHT)
        self.showF6.SetValue(sF[5])
        self.Bind(wx.EVT_CHECKBOX, self.on_labelF, self.showF6)   
        
        self.labelF1 = wx.CheckBox(self, -1, "Labels", style=wx.ALIGN_RIGHT)
        self.labelF1.SetValue(lF[0])
        self.Bind(wx.EVT_CHECKBOX, self.on_labelF, self.labelF1)    
        
        self.labelF2 = wx.CheckBox(self, -1, "Labels", style=wx.ALIGN_RIGHT)
        self.labelF2.SetValue(lF[1])
        self.Bind(wx.EVT_CHECKBOX, self.on_labelF, self.labelF2)
        
        self.labelF3 = wx.CheckBox(self, -1, "Labels", style=wx.ALIGN_RIGHT)
        self.labelF3.SetValue(lF[2])
        self.Bind(wx.EVT_CHECKBOX, self.on_labelF, self.labelF3)
        
        self.labelF4 = wx.CheckBox(self, -1, "Labels", style=wx.ALIGN_RIGHT)
        self.labelF4.SetValue(lF[3])
        self.Bind(wx.EVT_CHECKBOX, self.on_labelF, self.labelF4)                        

        self.labelF5 = wx.CheckBox(self, -1, "Labels", style=wx.ALIGN_RIGHT)
        self.labelF5.SetValue(lF[4])
        self.Bind(wx.EVT_CHECKBOX, self.on_labelF, self.labelF5)
        
        self.labelF6 = wx.CheckBox(self, -1, "Labels", style=wx.ALIGN_RIGHT)
        self.labelF6.SetValue(lF[5])
        self.Bind(wx.EVT_CHECKBOX, self.on_labelF, self.labelF6)        

        self.gbs.Add( self.quoteF1,   (0,0) )
        self.gbs.Add( self.quoteF2,   (1,0) )
        self.gbs.Add( self.quoteF3,   (2,0) )
        self.gbs.Add( self.quoteF4,   (3,0) )
        self.gbs.Add( self.quoteF5,   (4,0) )
        self.gbs.Add( self.quoteF6,   (5,0) )

        self.gbs.Add( self.showF1,   (0,1) )
        self.gbs.Add( self.showF2,   (1,1) )
        self.gbs.Add( self.showF3,   (2,1) )
        self.gbs.Add( self.showF4,   (3,1) )
        self.gbs.Add( self.showF5,   (4,1) )
        self.gbs.Add( self.showF6,   (5,1) )

        self.gbs.Add( self.labelF1,   (0,2) )
        self.gbs.Add( self.labelF2,   (1,2) )
        self.gbs.Add( self.labelF3,   (2,2) )
        self.gbs.Add( self.labelF4,   (3,2) )
        self.gbs.Add( self.labelF5,   (4,2) )
        self.gbs.Add( self.labelF6,   (5,2) )

        sizer.Add(self.gbs, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        
        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def on_labelF(self,event):
        pass

    def FGridUpdate(self, event):
        R = event.GetRow()
        C = event.GetCol()
        Value = float(self.gridForces1.GetCellValue(R, C))
        self.FOR[R][C] = Value
            
    def FGridRight(self,event):
        self.For_to_del = event.GetRow()
        FRMenu = wx.Menu()
        item = wx.MenuItem(FRMenu, wx.NewId(), "Delete Force")
        FRMenu.AppendItem(item)
        self.Bind(wx.EVT_MENU, self.OnDelFor, item)
        self.PopupMenu(FRMenu)
        
    def OnClickF(self,event):
        self.gridForces1.AppendRows(1, updateLabels = True)
        self.FOR = np.append(self.FOR, [[0 ,0.,0., 0.]], axis=0)
        
    def OnDelFor(self,event):
        self.gridForces1.DeleteRows(self.For_to_del, 1, updateLabels=True)
        self.FOR = np.delete(self.FOR, self.For_to_del, 0)
        for i in range(0,len(self.FOR[:,0])):
            if self.FOR[i][0] > self.For_to_del + 1:
                self.FOR[i][0] = self.FOR[i][0] - 1
                
                
class Materials(wx.Dialog):
    def __init__(self,MAT,PATH):
        title = 'Materials'
        
        self.MAT = MAT
        self.PATH = PATH
        
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)


        sizer = wx.BoxSizer(wx.VERTICAL)

        self.gridMaterials = wx.grid.Grid(self, size=(420,150))
        
        self.gridMaterials.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self.FGridUpdate)
        self.gridMaterials.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.ActiveRow)
        self.gridMaterials.Bind(wx.grid.EVT_GRID_LABEL_LEFT_CLICK, self.ActiveRow)
        self.gridMaterials.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.FGridRight)
        self.gridMaterials.Bind(wx.grid.EVT_GRID_LABEL_RIGHT_CLICK, self.FGridRight)           
        
        if len(self.MAT) != 0:
            self.gridMaterials.CreateGrid(len(self.MAT),3)
        else:
            self.gridMaterials.CreateGrid(0,3)
            
        self.gridMaterials.SetDefaultCellAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
        self.gridMaterials.SetRowLabelSize(100)
        self.gridMaterials.SetColLabelValue(0, 'Density')
        self.gridMaterials.SetColLabelValue(1, 'E')
        self.gridMaterials.SetColLabelValue(2, 'a')
        
        corn_M = self.gridMaterials.GetGridCornerLabelWindow()
        M_ID = wx.StaticText(corn_M, label="Material ID",pos=(9,8))
        font = wx.Font(9, wx.DEFAULT, wx.NORMAL, wx.BOLD)
        M_ID.SetFont(font)        
        
        for i in range(0,3):
            self.gridMaterials.SetColSize(i, 100)
        
        if len(self.MAT) != 0:
            for i in range(0,len(self.MAT)):
                if self.MAT[i][0] != 0:
                    self.gridMaterials.SetCellValue(i,0,str(self.MAT[i][1]))
                    self.gridMaterials.SetCellValue(i,1,str(self.MAT[i][2]))
                    self.gridMaterials.SetCellValue(i,2,str(self.MAT[i][3]))
                    

        sizer.Add(self.gridMaterials, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.buttonM1 = wx.Button(self, label="Add")
        self.Bind(wx.EVT_BUTTON, self.OnClickF, self.buttonM1)

        self.buttonM2 = wx.Button(self, label="Delete")
        self.Bind(wx.EVT_BUTTON, self.OnClickDel, self.buttonM2)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.buttonM1, 0, border=3, flag=wx.ALIGN_RIGHT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        hbox.Add(self.buttonM2, 0, border=3, flag=wx.ALIGN_RIGHT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)


        sizer.Add(hbox, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)



        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def ActiveRow(self,event):
        self.active_row = event.GetRow()
        event.Skip()

    def FGridUpdate(self, event):
        R = event.GetRow()
        C = event.GetCol()
        Value = float(self.gridMaterials.GetCellValue(R, C))
        self.MAT[R][C+1] = Value
            
    def FGridRight(self,event):
        self.Mat_to_del = event.GetRow()
        FRMenu = wx.Menu()
        item = wx.MenuItem(FRMenu, wx.NewId(), "Delete Material")
        FRMenu.AppendItem(item)
        self.Bind(wx.EVT_MENU, self.OnDelMat, item)
        self.PopupMenu(FRMenu)
        
    def OnClickF(self,event):
        self.gridMaterials.AppendRows(1, updateLabels = True)
        R = len(self.MAT)
        print(self.MAT)
        self.MAT.append([R+1, 0. ,0.,0.])

    def OnClickDel(self,event):
        self.Mat_to_del = self.active_row
        self.OnDelMat(0)
        
    def OnDelMat(self,event):
        self.gridMaterials.DeleteRows(self.Mat_to_del, 1, updateLabels=True)
        del self.MAT[self.Mat_to_del]
        for i in range(0,len(self.MAT)):
            if self.MAT[i][0] > self.Mat_to_del + 1:
                self.MAT[i][0] = self.MAT[i][0] - 1 
                
                
class Sections(wx.Dialog):
    def __init__(self,SEC,PATH):
        title = 'Sections'
        
        self.SEC = SEC
        self.PATH = PATH
        
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)


        sizer = wx.BoxSizer(wx.VERTICAL)

        self.gridSections = wx.grid.Grid(self, size=(520,150))
        
        self.gridSections.Bind(wx.grid.EVT_GRID_CELL_CHANGED, self.FGridUpdate)
        self.gridSections.Bind(wx.grid.EVT_GRID_CELL_LEFT_CLICK, self.ActiveRow)
        self.gridSections.Bind(wx.grid.EVT_GRID_LABEL_LEFT_CLICK, self.ActiveRow)
        self.gridSections.Bind(wx.grid.EVT_GRID_CELL_RIGHT_CLICK, self.FGridRight)
        self.gridSections.Bind(wx.grid.EVT_GRID_LABEL_RIGHT_CLICK, self.FGridRight)           
        
        if len(self.SEC) != 0:
            self.gridSections.CreateGrid(len(self.SEC),4)
        else:
            self.gridSections.CreateGrid(0,4)
            
        self.gridSections.SetDefaultCellAlignment(wx.ALIGN_CENTRE, wx.ALIGN_CENTRE )
        self.gridSections.SetRowLabelSize(100)
        self.gridSections.SetColLabelValue(0, 'A')
        self.gridSections.SetColLabelValue(1, 'I')
        self.gridSections.SetColLabelValue(2, 'h')
        self.gridSections.SetColLabelValue(3, 'Material ID')
        
        corn_M = self.gridSections.GetGridCornerLabelWindow()
        M_ID = wx.StaticText(corn_M, label="Section ID",pos=(9,8))
        font = wx.Font(9, wx.DEFAULT, wx.NORMAL, wx.BOLD)
        M_ID.SetFont(font)        
        
        for i in range(0,4):
            self.gridSections.SetColSize(i, 100)
        
        if len(self.SEC) != 0:
            for i in range(0,len(self.SEC)):
                if self.SEC[i][0] != 0:
                    self.gridSections.SetCellValue(i,0,str(self.SEC[i][1]))
                    self.gridSections.SetCellValue(i,1,str(self.SEC[i][2]))
                    self.gridSections.SetCellValue(i,2,str(self.SEC[i][3]))
                    self.gridSections.SetCellValue(i,3,str(self.SEC[i][4]))
                    

        sizer.Add(self.gridSections, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.buttonM1 = wx.Button(self, label="Add")
        self.Bind(wx.EVT_BUTTON, self.OnClickF, self.buttonM1)

        self.buttonM2 = wx.Button(self, label="Delete")
        self.Bind(wx.EVT_BUTTON, self.OnClickDel, self.buttonM2)

        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(self.buttonM1, 0, border=3, flag=wx.ALIGN_RIGHT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)
        hbox.Add(self.buttonM2, 0, border=3, flag=wx.ALIGN_RIGHT | wx.ALL | wx.ALIGN_CENTER_VERTICAL)


        sizer.Add(hbox, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)



        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def ActiveRow(self,event):
        self.active_row = event.GetRow()
        event.Skip()

    def FGridUpdate(self, event):
        R = event.GetRow()
        C = event.GetCol()
        Value = float(self.gridSections.GetCellValue(R, C))
        self.SEC[R][C+1] = Value
            
    def FGridRight(self,event):
        self.Sec_to_del = event.GetRow()
        FRMenu = wx.Menu()
        item = wx.MenuItem(FRMenu, wx.NewId(), "Delete Section")
        FRMenu.AppendItem(item)
        self.Bind(wx.EVT_MENU, self.OnDelSec, item)
        self.PopupMenu(FRMenu)
        
    def OnClickF(self,event):
        self.gridSections.AppendRows(1, updateLabels = True)
        R = len(self.SEC)
        self.SEC.append([R+1, 0. ,0.,0., 0.])

    def OnClickDel(self,event):
        self.Sec_to_del = self.active_row
        self.OnDelSec(0)
        
    def OnDelSec(self,event):
        self.gridSections.DeleteRows(self.Sec_to_del, 1, updateLabels=True)
        del self.SEC[self.Sec_to_del]
        for i in range(0,len(self.SEC)):
            if self.SEC[i][0] > self.Sec_to_del + 1:
                self.SEC[i][0] = self.SEC[i][0] - 1 
                
class Support(wx.Dialog):
    def __init__(self,dat):
        title = 'Support'
        
        self.dat = dat
        
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)
        
        self.triple = [self.dat[1],self.dat[2],self.dat[3]]
        self.sup_angle = self.dat[4]
        sizer = wx.BoxSizer(wx.VERTICAL)
        Gsizer = wx.GridBagSizer(5,5)

        Gsizer.Add(wx.StaticText(self,-1,'X Translation:   '),(0,0))
        Gsizer.Add(wx.StaticText(self,-1,'Y Translation:   '),(1,0))
        Gsizer.Add(wx.StaticText(self,-1,'Rotation:   '),(2,0))
        
        Gsizer.Add(wx.StaticText(self,-1,'Angle('u'\N{DEGREE SIGN} ):   '),(4,0))
        
        sampleList = ['0','1']
        self.cb1 = wx.ComboBox(self, value=str(self.dat[1]), choices = sampleList,style = wx.CB_READONLY)
        self.cb2 = wx.ComboBox(self, value=str(self.dat[2]), choices = sampleList,style = wx.CB_READONLY)
        self.cb3 = wx.ComboBox(self, value=str(self.dat[3]), choices = sampleList,style = wx.CB_READONLY)

        self.tc4 = wx.TextCtrl(self, -1, str(self.dat[4]), size=(50, -1))

        self.Bind(wx.EVT_COMBOBOX, self.EvtComboBox, self.cb1)
        self.Bind(wx.EVT_COMBOBOX, self.EvtComboBox, self.cb2)
        self.Bind(wx.EVT_COMBOBOX, self.EvtComboBox, self.cb3)

        self.Bind(wx.EVT_TEXT, self.EvtText, self.tc4)


        Gsizer.Add(self.cb1, (0,4))
        Gsizer.Add(self.cb2, (1,4))
        Gsizer.Add(self.cb3, (2,4))
        
        Gsizer.Add(self.tc4, (4,4))

        sizer.Add(Gsizer, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def EvtText(self, event):
        self.sup_angle = float(event.GetString())
        event.Skip()

    def EvtComboBox(self, event):
        self.triple = [self.cb1.GetSelection(),self.cb2.GetSelection(),self.cb3.GetSelection()]
        event.Skip()
        

class EditElement(wx.Dialog):
    def __init__(self,dat):
        title = 'Edit Element'
        
        self.dat = dat
        
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        Gsizer = wx.GridBagSizer(5,5)

        if len(self.dat) == 1:
            self.dat_out = [0,0,0,0,0]
            
            Gsizer.Add(wx.StaticText(self,-1,'X'),(0,1))
            Gsizer.Add(wx.StaticText(self,-1,'Y'),(0,2))
            Gsizer.Add(wx.StaticText(self,-1,'Node 1 :   '),(1,0))
            Gsizer.Add(wx.StaticText(self,-1,'Node 2 :   '),(2,0))
    
            Gsizer.Add(wx.StaticText(self,-1,'Section ID:   '),(4,0))
    
            self.tc1 = wx.TextCtrl(self, -1, str(self.dat[0][0]), size=(50, -1))
            self.tc2 = wx.TextCtrl(self, -1, str(self.dat[0][1]), size=(50, -1))
            self.tc3 = wx.TextCtrl(self, -1, str(self.dat[0][2]), size=(50, -1))
            self.tc4 = wx.TextCtrl(self, -1, str(self.dat[0][3]), size=(50, -1))
            
            self.tc5 = wx.TextCtrl(self, -1, str(self.dat[0][4]), size=(50, -1))
    
            self.Bind(wx.EVT_TEXT, self.EvtText, self.tc1)
            self.Bind(wx.EVT_TEXT, self.EvtText, self.tc2)
            self.Bind(wx.EVT_TEXT, self.EvtText, self.tc3)
            self.Bind(wx.EVT_TEXT, self.EvtText, self.tc4)
            
            self.Bind(wx.EVT_TEXT, self.EvtText, self.tc5)
    
            Gsizer.Add(self.tc1, (1,1))
            Gsizer.Add(self.tc2, (1,2))
            Gsizer.Add(self.tc3, (2,1))
            Gsizer.Add(self.tc4, (2,2))
            
            Gsizer.Add(self.tc5, (4,1))
        
        else:
            Gsizer.Add(wx.StaticText(self,-1,'Section ID:   '),(0,0))
            self.tc5 = wx.TextCtrl(self, -1, str(self.dat[0][4]), size=(50, -1))
            self.Bind(wx.EVT_TEXT, self.EvtText, self.tc5)
            Gsizer.Add(self.tc5, (0,1))
            
            
        sizer.Add(Gsizer, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def EvtText(self, event):
        if len(self.dat) == 1:
            self.dat_out[0] = float(self.tc1.GetValue())
            self.dat_out[1] = float(self.tc2.GetValue())
            self.dat_out[2] = float(self.tc3.GetValue())
            self.dat_out[3] = float(self.tc4.GetValue())
            
            self.dat_out[4] = float(self.tc5.GetValue())
        else:
            self.dat_out = float(self.tc5.GetValue())

        event.Skip()
# 
#     def EvtComboBox(self, event):
#         self.triple = [self.cb1.GetSelection(),self.cb2.GetSelection(),self.cb3.GetSelection()]
#         event.Skip()
        
class F1(wx.Dialog):
    def __init__(self,dat, PATH):
        title = 'Loads on Node'
        
        self.dat = dat
        self.PATH = PATH
        
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)
        
        self.dat_out = [self.dat[1],self.dat[2],self.dat[3]]
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        Gsizer = wx.GridBagSizer(5,5)
            
        Gsizer.Add(wx.StaticText(self,-1,'X Force: '),(0,0))
        Gsizer.Add(wx.StaticText(self,-1,'Y Force: '),(1,0))
        Gsizer.Add(wx.StaticText(self,-1,'Moment:  '),(2,0))

        self.tc1 = wx.TextCtrl(self, -1, str(self.dat[1]), size=(50, -1))
        self.tc2 = wx.TextCtrl(self, -1, str(self.dat[2]), size=(50, -1))
        self.tc3 = wx.TextCtrl(self, -1, str(self.dat[3]), size=(50, -1))

        self.Bind(wx.EVT_TEXT, self.EvtText, self.tc1)
        self.Bind(wx.EVT_TEXT, self.EvtText, self.tc2)
        self.Bind(wx.EVT_TEXT, self.EvtText, self.tc3)


        Gsizer.Add(self.tc1, (0,4))
        Gsizer.Add(self.tc2, (1,4))
        Gsizer.Add(self.tc3, (2,4))
            
        sizer.Add(Gsizer, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        Img = wx.StaticBitmap(self, bitmap=wx.Image(self.PATH + '/icons/pos_nf.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap())

        sizer.Add(Img, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def EvtText(self, event):
        self.dat_out[0] = float(self.tc1.GetValue())
        self.dat_out[1] = float(self.tc2.GetValue())
        self.dat_out[2] = float(self.tc3.GetValue())

        event.Skip()
        
        
class F2(wx.Dialog):
    def __init__(self,dat, PATH):
        title = 'Loads on Element'
        
        self.dat = dat
        self.PATH = PATH
        
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)
        
        self.dat_out = [self.dat[1],self.dat[2],self.dat[3],self.dat[4]]
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        Gsizer = wx.GridBagSizer(5,5)
            
        Gsizer.Add(wx.StaticText(self,-1,'    X Force: '),(0,0))
        Gsizer.Add(wx.StaticText(self,-1,'    Y Force: '),(1,0))
        Gsizer.Add(wx.StaticText(self,-1,'    Moment:  '),(2,0))
        Gsizer.Add(wx.StaticText(self,-1,'    Relative Position:  '),(3,0))

        self.tc1 = wx.TextCtrl(self, -1, str(self.dat[1]), size=(50, -1))
        self.tc2 = wx.TextCtrl(self, -1, str(self.dat[2]), size=(50, -1))
        self.tc3 = wx.TextCtrl(self, -1, str(self.dat[3]), size=(50, -1))
        self.tc4 = wx.TextCtrl(self, -1, str(self.dat[4]), size=(50, -1))

        self.Bind(wx.EVT_TEXT, self.EvtText, self.tc1)
        self.Bind(wx.EVT_TEXT, self.EvtText, self.tc2)
        self.Bind(wx.EVT_TEXT, self.EvtText, self.tc3)
        self.Bind(wx.EVT_TEXT, self.EvtText, self.tc4)


        Gsizer.Add(self.tc1, (0,4))
        Gsizer.Add(self.tc2, (1,4))
        Gsizer.Add(self.tc3, (2,4))
        Gsizer.Add(self.tc4, (3,4))
            
        sizer.Add(Gsizer, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        
        sizer3 = wx.BoxSizer(wx.VERTICAL)
        sizer4 = wx.BoxSizer(wx.VERTICAL)
        
        self.Local = True
        self.radio1 = wx.RadioButton(self, -1, " Local CS ", style = wx.RB_GROUP )
        self.radio2 = wx.RadioButton(self, -1, " Global CS " )
        self.Bind(wx.EVT_RADIOBUTTON, self.OnRadioSelect, self.radio1 )
        self.Bind(wx.EVT_RADIOBUTTON, self.OnRadioSelect, self.radio2 )

        Img1 = wx.StaticBitmap(self, bitmap=wx.Image(self.PATH + '/icons/pos_ef.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        Img2 = wx.StaticBitmap(self, bitmap=wx.Image(self.PATH + '/icons/pos_ef_glob.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap())

        sizer3.Add(self.radio1, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
        sizer3.Add(Img1, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
        
        sizer4.Add(self.radio2, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
        sizer4.Add(Img2, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
        
        sizer2.Add(sizer3, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        line = wx.StaticLine(self, -1, size=(-1,20), style=wx.LI_VERTICAL)
        sizer2.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        sizer2.Add(sizer4, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        sizer.Add(sizer2, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def OnRadioSelect(self,event):
        radio_selected = event.GetEventObject()
        if radio_selected is self.radio1:
            self.Local = True
        elif radio_selected is self.radio2:
            self.Local = False

    def EvtText(self, event):
        self.dat_out[0] = float(self.tc1.GetValue())
        self.dat_out[1] = float(self.tc2.GetValue())
        self.dat_out[2] = float(self.tc3.GetValue())
        self.dat_out[3] = float(self.tc4.GetValue())

        event.Skip()
        
class F3(wx.Dialog):
    def __init__(self,dat, PATH):
        title = 'Distributed Load'
        
        self.dat = dat
        self.PATH = PATH
        
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)
        
        self.dat_out = [self.dat[1],self.dat[2],self.dat[3]]
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        Gsizer = wx.GridBagSizer(5,5)
            
        Gsizer.Add(wx.StaticText(self,-1,'    X Load: '),(0,0))
        Gsizer.Add(wx.StaticText(self,-1,'    Y Load: '),(1,0))
        Gsizer.Add(wx.StaticText(self,-1,'    Moment:  '),(2,0))

        self.tc1 = wx.TextCtrl(self, -1, str(self.dat[1]), size=(50, -1))
        self.tc2 = wx.TextCtrl(self, -1, str(self.dat[2]), size=(50, -1))
        self.tc3 = wx.TextCtrl(self, -1, str(self.dat[3]), size=(50, -1))

        self.Bind(wx.EVT_TEXT, self.EvtText, self.tc1)
        self.Bind(wx.EVT_TEXT, self.EvtText, self.tc2)
        self.Bind(wx.EVT_TEXT, self.EvtText, self.tc3)

        Gsizer.Add(self.tc1, (0,4))
        Gsizer.Add(self.tc2, (1,4))
        Gsizer.Add(self.tc3, (2,4))
            
        sizer.Add(Gsizer, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        sizer2 = wx.BoxSizer(wx.HORIZONTAL)
        
        sizer3 = wx.BoxSizer(wx.VERTICAL)
        sizer4 = wx.BoxSizer(wx.VERTICAL)
        
        self.Local = True
        self.radio1 = wx.RadioButton(self, -1, " Local CS ", style = wx.RB_GROUP )
        self.radio2 = wx.RadioButton(self, -1, " Global CS " )
        self.Bind(wx.EVT_RADIOBUTTON, self.OnRadioSelect, self.radio1 )
        self.Bind(wx.EVT_RADIOBUTTON, self.OnRadioSelect, self.radio2 )

        Img1 = wx.StaticBitmap(self, bitmap=wx.Image(self.PATH + '/icons/pos_ef.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap())
        Img2 = wx.StaticBitmap(self, bitmap=wx.Image(self.PATH + '/icons/pos_ef_glob.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap())

        sizer3.Add(self.radio1, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
        sizer3.Add(Img1, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
        
        sizer4.Add(self.radio2, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
        sizer4.Add(Img2, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
        
        sizer2.Add(sizer3, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        line = wx.StaticLine(self, -1, size=(-1,20), style=wx.LI_VERTICAL)
        sizer2.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        sizer2.Add(sizer4, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        sizer.Add(sizer2, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def OnRadioSelect(self,event):
        radio_selected = event.GetEventObject()
        if radio_selected is self.radio1:
            self.Local = True
        elif radio_selected is self.radio2:
            self.Local = False

    def EvtText(self, event):
        self.dat_out[0] = float(self.tc1.GetValue())
        self.dat_out[1] = float(self.tc2.GetValue())
        self.dat_out[2] = float(self.tc3.GetValue())

        event.Skip()
        
class F4(wx.Dialog):
    def __init__(self,dat, PATH):
        title = 'Temperature Load'
        
        self.dat = dat
        self.PATH = PATH
        
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)
        
        self.dat_out = [self.dat[1],self.dat[2],self.dat[3]]
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        Gsizer = wx.GridBagSizer(5,5)
            
        Gsizer.Add(wx.StaticText(self,-1,'    T in: '),(0,0))
        Gsizer.Add(wx.StaticText(self,-1,'    T out: '),(1,0))
        Gsizer.Add(wx.StaticText(self,-1,'    T 0:  '),(2,0))

        self.tc1 = wx.TextCtrl(self, -1, str(self.dat[1]), size=(50, -1))
        self.tc2 = wx.TextCtrl(self, -1, str(self.dat[2]), size=(50, -1))
        self.tc3 = wx.TextCtrl(self, -1, str(self.dat[3]), size=(50, -1))

        self.Bind(wx.EVT_TEXT, self.EvtText, self.tc1)
        self.Bind(wx.EVT_TEXT, self.EvtText, self.tc2)
        self.Bind(wx.EVT_TEXT, self.EvtText, self.tc3)

        Gsizer.Add(self.tc1, (0,4))
        Gsizer.Add(self.tc2, (1,4))
        Gsizer.Add(self.tc3, (2,4))
            
        sizer.Add(Gsizer, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

#         sizer2 = wx.BoxSizer(wx.HORIZONTAL)
#         
#         sizer3 = wx.BoxSizer(wx.VERTICAL)
#         sizer4 = wx.BoxSizer(wx.VERTICAL)
#         
#         self.Local = True
#         self.radio1 = wx.RadioButton(self, -1, " Local CS ", style = wx.RB_GROUP )
#         self.radio2 = wx.RadioButton(self, -1, " Global CS " )
#         self.Bind(wx.EVT_RADIOBUTTON, self.OnRadioSelect, self.radio1 )
#         self.Bind(wx.EVT_RADIOBUTTON, self.OnRadioSelect, self.radio2 )

        Img1 = wx.StaticBitmap(self, bitmap=wx.Image(self.PATH + '/icons/pos_tl.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap())
#         Img2 = wx.StaticBitmap(self, bitmap=wx.Image(self.PATH + '/icons/pos_ef_glob.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap())

#         sizer3.Add(self.radio1, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
#         sizer3.Add(Img1, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
#         
#         sizer4.Add(self.radio2, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
#         sizer4.Add(Img2, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
#         
#         sizer2.Add(sizer3, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
#         line = wx.StaticLine(self, -1, size=(-1,20), style=wx.LI_VERTICAL)
#         sizer2.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
#         sizer2.Add(sizer4, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        sizer.Add(Img1, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def OnRadioSelect(self,event):
        radio_selected = event.GetEventObject()
        if radio_selected is self.radio1:
            self.Local = True
        elif radio_selected is self.radio2:
            self.Local = False

    def EvtText(self, event):
        self.dat_out[0] = float(self.tc1.GetValue())
        self.dat_out[1] = float(self.tc2.GetValue())
        self.dat_out[2] = float(self.tc3.GetValue())

        event.Skip()
        
        
class F5(wx.Dialog):
    def __init__(self,dat, PATH):
        title = 'Temperature Load'
        
        self.dat = dat
        self.PATH = PATH
        
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)
        
        self.dat_out = [self.dat[1],self.dat[2],self.dat[3]]
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        Gsizer = wx.GridBagSizer(5,5)
            
        Gsizer.Add(wx.StaticText(self,-1,'    Delta X: '),(0,0))
        Gsizer.Add(wx.StaticText(self,-1,'    Delta Y: '),(1,0))
        Gsizer.Add(wx.StaticText(self,-1,'    Delta phi:  '),(2,0))

        self.tc1 = wx.TextCtrl(self, -1, str(self.dat[1]), size=(50, -1))
        self.tc2 = wx.TextCtrl(self, -1, str(self.dat[2]), size=(50, -1))
        self.tc3 = wx.TextCtrl(self, -1, str(self.dat[3]), size=(50, -1))

        self.Bind(wx.EVT_TEXT, self.EvtText, self.tc1)
        self.Bind(wx.EVT_TEXT, self.EvtText, self.tc2)
        self.Bind(wx.EVT_TEXT, self.EvtText, self.tc3)

        Gsizer.Add(self.tc1, (0,4))
        Gsizer.Add(self.tc2, (1,4))
        Gsizer.Add(self.tc3, (2,4))
            
        sizer.Add(Gsizer, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

#         sizer2 = wx.BoxSizer(wx.HORIZONTAL)
#         
#         sizer3 = wx.BoxSizer(wx.VERTICAL)
#         sizer4 = wx.BoxSizer(wx.VERTICAL)
#         
#         self.Local = True
#         self.radio1 = wx.RadioButton(self, -1, " Local CS ", style = wx.RB_GROUP )
#         self.radio2 = wx.RadioButton(self, -1, " Global CS " )
#         self.Bind(wx.EVT_RADIOBUTTON, self.OnRadioSelect, self.radio1 )
#         self.Bind(wx.EVT_RADIOBUTTON, self.OnRadioSelect, self.radio2 )

        Img1 = wx.StaticBitmap(self, bitmap=wx.Image(self.PATH + '/icons/pos_nf.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap())
#         Img2 = wx.StaticBitmap(self, bitmap=wx.Image(self.PATH + '/icons/pos_ef_glob.png',wx.BITMAP_TYPE_PNG).ConvertToBitmap())

#         sizer3.Add(self.radio1, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
#         sizer3.Add(Img1, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
#         
#         sizer4.Add(self.radio2, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
#         sizer4.Add(Img2, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
#         
#         sizer2.Add(sizer3, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
#         line = wx.StaticLine(self, -1, size=(-1,20), style=wx.LI_VERTICAL)
#         sizer2.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
#         sizer2.Add(sizer4, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        sizer.Add(Img1, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        line = wx.StaticLine(self, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.SetSizer(sizer)
        sizer.Fit(self)

    def OnRadioSelect(self,event):
        radio_selected = event.GetEventObject()
        if radio_selected is self.radio1:
            self.Local = True
        elif radio_selected is self.radio2:
            self.Local = False

    def EvtText(self, event):
        self.dat_out[0] = float(self.tc1.GetValue())
        self.dat_out[1] = float(self.tc2.GetValue())
        self.dat_out[2] = float(self.tc3.GetValue())

        event.Skip()
        
class Analysis_Options(wx.Dialog):
    def __init__(self,dat,n_modes):
        title = 'Analysis Options'
        
        self.dat = dat
        self.n_modes = n_modes
        
        wx.Dialog.__init__(self, None, -1, title, size=wx.DefaultSize, pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE)

        self.nb = wx.Notebook(self)
        
        self.StaticTab = wx.Window(self.nb)
        self.ModalTab = wx.Window(self.nb)
        self.StaticTab.SetBackgroundColour("white")
        self.ModalTab.SetBackgroundColour("white")
        
        self.nb.AddPage(self.StaticTab, "Static Analysis")
        self.nb.AddPage(self.ModalTab, "Modal Analysis")
        
        self.dat_out = self.dat
        
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer1 = wx.BoxSizer(wx.VERTICAL)
        
        sizer1.AddSpacer(10)
        sizer1.Add(wx.StaticText(self.StaticTab,-1,'    Solution of the linear system using: '))
        sizer1.AddSpacer(10)
        
        self.group1_ctrls = []
        radio1 = wx.RadioButton( self.StaticTab, -1, " Gaussian elimination  ", style = wx.RB_GROUP )
        radio2 = wx.RadioButton( self.StaticTab, -1, " LU decomposition  " )
        
        self.group1_ctrls.append((radio1))
        self.group1_ctrls.append((radio2))
        
        if self.dat[0] == 1:
            radio1.SetValue(True)
        elif self.dat[0] == 2:
            radio2.SetValue(True)
        
        Gsizer = wx.GridBagSizer(5,5)
        Gsizer.Add(radio1,(0,0))
        Gsizer.Add(radio2,(0,1))
        sizer1.Add(Gsizer, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        line = wx.StaticLine(self.StaticTab, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer1.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        sizer1.AddSpacer(10)           
        
        hbox = wx.BoxSizer(wx.HORIZONTAL)
        hbox.Add(wx.StaticText(self.StaticTab,-1,' Calculation of diagrams in '))
        acc = wx.TextCtrl( self.StaticTab, -1, str(int(self.dat[11])), size=(50, -1)  )
        self.Bind(wx.EVT_TEXT, self.EvtAcc, acc)
        hbox.Add(acc)
        hbox.Add(wx.StaticText(self.StaticTab,-1,' points per element.'))
        
        sizer1.Add(hbox, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
        
        sizer2 = wx.BoxSizer(wx.VERTICAL)
        
        sizer2.AddSpacer(10)
        sizer2.Add(wx.StaticText(self.ModalTab,-1,'    Solution method for the eigenvalue problem: '))
        sizer2.AddSpacer(10)
        
        self.group2_ctrls = []
        radio10 = wx.RadioButton( self.ModalTab, -1, " Jacobi  ", style = wx.RB_GROUP )
        radio20 = wx.RadioButton( self.ModalTab, -1, " Subspace iteration  " )
        radio30 = wx.RadioButton( self.ModalTab, -1, " HQRI  " )

        self.group2_ctrls.append((radio10))
        self.group2_ctrls.append((radio20))
        self.group2_ctrls.append((radio30))

        if self.dat[1] == 1:
            radio10.SetValue(True)
        elif self.dat[1] == 2:
            radio20.SetValue(True)
        elif self.dat[1] == 3:
            radio30.SetValue(True)
        
        Gsizer2 = wx.GridBagSizer(5,5)
        Gsizer2.Add(radio10,(0,0))
        Gsizer2.Add(radio20,(0,1))
        Gsizer2.Add(radio30,(1,0))
        sizer2.Add(Gsizer2, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
        
        samplelist = ['all']
        for i in range(self.n_modes-1,1,-1):
            samplelist.append(str(i))
            
        self.choice = wx.Choice(self.ModalTab, -1, (150, 150), choices = samplelist)
        self.Bind(wx.EVT_CHOICE, self.EvtChoice, self.choice)
        
        if self.dat[2] == 0:
            self.choice.SetSelection(0)
        else:
            for i in range(len(samplelist)-1,0,-1):
                if int(samplelist[i]) == self.dat[2]:
                    self.choice.SetSelection(i)

        hbox2 = wx.BoxSizer(wx.HORIZONTAL)
        quoteM = wx.StaticText(self.ModalTab,-1,' Number of modes to calculate: ')
        hbox2.Add(quoteM)
        hbox2.Add(self.choice)
        
        sizer2.Add(hbox2, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        line = wx.StaticLine(self.ModalTab, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer2.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)
        sizer2.Add(wx.StaticText(self.ModalTab,-1,'   Parameters of the Jacobi method: '))

        Gsizer3 = wx.GridBagSizer(5,5)
        quotePJ1 = wx.StaticText(self.ModalTab,-1,'        Convergence tolerance: ')
        self.textPJ1 = wx.TextCtrl( self.ModalTab, -1, str(self.dat[3]), size=(50, -1)  )
        self.Bind(wx.EVT_TEXT, self.EvtTextPJ1, self.textPJ1)
        quotePJ2 = wx.StaticText(self.ModalTab,-1,'        Maximum number of iterations: ')
        self.textPJ2 = wx.TextCtrl( self.ModalTab, -1, str(self.dat[4]), size=(50, -1)  )
        self.Bind(wx.EVT_TEXT, self.EvtTextPJ2, self.textPJ2)
        Gsizer3.Add(quotePJ1,(0,0))
        Gsizer3.Add(self.textPJ1,(0,1))
        Gsizer3.Add(quotePJ2,(1,0))
        Gsizer3.Add(self.textPJ2,(1,1))

        sizer2.Add(Gsizer3, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        line = wx.StaticLine(self.ModalTab, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer2.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)        
        sizer2.Add(wx.StaticText(self.ModalTab,-1,'   Parameters of the Subspace Iteration method: '))
        
        Gsizer4 = wx.GridBagSizer(5,5)
        quotePS1 = wx.StaticText(self.ModalTab,-1,'        Convergence tolerance: ')
        self.textPS1 = wx.TextCtrl( self.ModalTab, -1, str(self.dat[5]), size=(50, -1)  )
        self.Bind(wx.EVT_TEXT, self.EvtTextPS1, self.textPS1)
        quotePS2 = wx.StaticText(self.ModalTab,-1,'        Maximum number of iterations: ')
        self.textPS2 = wx.TextCtrl( self.ModalTab, -1, str(self.dat[6]), size=(50, -1)  )
        self.Bind(wx.EVT_TEXT, self.EvtTextPS2, self.textPS2)
        Gsizer4.Add(quotePS1,(0,0))
        Gsizer4.Add(self.textPS1,(0,1))
        Gsizer4.Add(quotePS2,(1,0))
        Gsizer4.Add(self.textPS2,(1,1))
        
        sizer2.Add(Gsizer4, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)
        
        line = wx.StaticLine(self.ModalTab, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
        sizer2.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)        
        sizer2.Add(wx.StaticText(self.ModalTab,-1,'   Parameters of the HQRI algorithm: '))
        
        Gsizer5 = wx.GridBagSizer(5,5)
        quotePQ1 = wx.StaticText(self.ModalTab,-1,'        Convergence tolerance: ')
        self.textPQ1 = wx.TextCtrl( self.ModalTab, -1, str(self.dat[7]), size=(50, -1)  )
        self.Bind(wx.EVT_TEXT, self.EvtTextPQ1, self.textPQ1)
        quotePQ2 = wx.StaticText(self.ModalTab,-1,'        Maximum number of iterations: ')
        self.textPQ2 = wx.TextCtrl( self.ModalTab, -1, str(self.dat[8]), size=(50, -1)  )
        self.Bind(wx.EVT_TEXT, self.EvtTextPQ2, self.textPQ2)
        Gsizer5.Add(quotePQ1,(0,0))
        Gsizer5.Add(self.textPQ1,(0,1))
        Gsizer5.Add(quotePQ2,(1,0))
        Gsizer5.Add(self.textPQ2,(1,1))
        
        sizer2.Add(Gsizer5, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

#        line = wx.StaticLine(self.ModalTab, -1, size=(20,-1), style=wx.LI_HORIZONTAL)
#        sizer2.Add(line, 0, wx.GROW|wx.ALIGN_CENTER_VERTICAL|wx.RIGHT|wx.TOP, 5)        
#        sizer2.Add(wx.StaticText(self.ModalTab,-1,'   Parameters of the Inverse Iteration algorithm: '))
#        
#        Gsizer6 = wx.GridBagSizer(5,5)
#        quotePI1 = wx.StaticText(self.ModalTab,-1,'        Convergence tolerance: ')
#        self.textPI1 = wx.TextCtrl( self.ModalTab, -1, str(self.dat[9]), size=(50, -1) )
#        self.Bind(wx.EVT_TEXT, self.EvtTextPQ1, self.textPQ1)
#        quotePI2 = wx.StaticText(self.ModalTab,-1,'        Maximum number of iterations: ')
#        self.textPI2 = wx.TextCtrl( self.ModalTab, -1, str(self.dat[10]), size=(50, -1) )
#        self.Bind(wx.EVT_TEXT, self.EvtTextPI2, self.textPI2)
#        Gsizer6.Add(quotePI1,(0,0))
#        Gsizer6.Add(self.textPI1,(0,1))
#        Gsizer6.Add(quotePI2,(1,0))
#        Gsizer6.Add(self.textPI2,(1,1))
#        
#        sizer2.Add(Gsizer6, 0, wx.GROW|wx.ALIGN_RIGHT|wx.ALL, 5)

        btnsizer = wx.StdDialogButtonSizer()
        
        if wx.Platform != "__WXMSW__":
            btn = wx.ContextHelpButton(self)
            btnsizer.AddButton(btn)
        
        btn = wx.Button(self, wx.ID_OK)
        btn.SetDefault()
        btnsizer.AddButton(btn)

        btn = wx.Button(self, wx.ID_CANCEL)
        btnsizer.AddButton(btn)
        btnsizer.Realize()

        sizer.Add(self.nb, 5, flag =  wx.ALL  | wx.ALIGN_CENTER_VERTICAL)
        sizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL|wx.ALL, 5)

        self.StaticTab.SetSizer(sizer1)
        self.ModalTab.SetSizer(sizer2)

        self.SetSizer(sizer)
        sizer.Fit(self)

        self.OnGroup2Select( 0 )

        for radio in self.group1_ctrls:
            self.Bind(wx.EVT_RADIOBUTTON, self.OnGroup1Select, radio )

        for radio in self.group2_ctrls:
            self.Bind(wx.EVT_RADIOBUTTON, self.OnGroup2Select, radio )

    def OnGroup1Select( self, event ):
#         radio_selected = event.GetEventObject()
        radio = self.group1_ctrls[0]
        if radio.GetValue() == True:
            self.dat_out[0] = 1
        
        radio = self.group1_ctrls[1]    
        if radio.GetValue() == True:
            self.dat_out[0] = 2
            
    def OnGroup2Select( self, event ):
#         radio_selected = event.GetEventObject()
        radio = self.group2_ctrls[0]
        if radio.GetValue() == True:
            self.dat_out[1] = 1
            self.textPJ1.Enable(True)
            self.textPJ2.Enable(True)
            self.choice.SetSelection(0)
            self.choice.Disable()
        else:
            self.textPJ1.Enable(False)
            self.textPJ2.Enable(False)
            self.choice.Enable()

        radio = self.group2_ctrls[1]
        if radio.GetValue() == True:
            self.dat_out[1] = 2
            self.textPS1.Enable(True)
            self.textPS2.Enable(True)
#            self.textPJ1.Enable(True)
#            self.textPJ2.Enable(True)
        else:
            self.textPS1.Enable(False)
            self.textPS2.Enable(False)
            
        radio = self.group2_ctrls[2]
        if radio.GetValue() == True:
            self.dat_out[1] = 3
            self.textPQ1.Enable(True)
            self.textPQ2.Enable(True)
#            self.textPI1.Enable(True)
#            self.textPI2.Enable(True)
        else:
            self.textPQ1.Enable(False)
            self.textPQ2.Enable(False)
#            self.textPI1.Enable(False)
#            self.textPI2.Enable(False)
            
    def EvtTextPJ1(self,event):
        self.dat_out[3] = float(event.GetString())
        
    def EvtTextPJ2(self,event):
        self.dat_out[4] = int(event.GetString())
        
    def EvtTextPS1(self,event):
        self.dat_out[5] = float(event.GetString())
        
    def EvtTextPS2(self,event):
        self.dat_out[6] = int(event.GetString())
        
    def EvtTextPQ1(self,event):
        self.dat_out[7] = float(event.GetString())
        
    def EvtTextPQ2(self,event):
        self.dat_out[8] = int(event.GetString())
        
    def EvtTextPI1(self,event):
        self.dat_out[9] = float(event.GetString())
        
    def EvtTextPI2(self,event):
        self.dat_out[10] = int(event.GetString())
        
    def EvtChoice(self,event):
        val = event.GetString()
        if val == 'all' :
            self.dat_out[2] = 0
        else:
            self.dat_out[2] = int(val)
    
    def EvtAcc(self,event):
        self.dat_out[11] = int(event.GetString())