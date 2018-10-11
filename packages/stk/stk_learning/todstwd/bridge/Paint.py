from tkinter import *
from tkinter import messagebox
import tkinter.font
from PIL import ImageTk,Image,ImageDraw
import os
import sys
import time

def strnum(i):
  num = ""
  if (i  < 10):
    num = "0000" + str(i)
  elif (i < 100):
    num = "000" + str(i)
  elif (i < 1000):
    num = "00" + str(i)
  elif (i < 10000):
    num = "0" + str(i)
  else:
    num = str(i)
  return num

def file_cleanup():
    os.system('rm -f my_*')
    os.system('rm -f *.jou')
    os.system('rm -f bridgeFalls.e* bridgeLoads.e*')
    os.system('rm -f *.log')

def delete_frags():
    os.system("module purge;module load sierra-devel;python delete_small_frags.py")

def get_points():
    # Max Pts = 331
    pts = 0
    fname = 'bridgeFull.log'
    if os.path.isfile(fname):
      f = open(fname,'r')
      lines = f.readlines()
      for line in lines:
          if re.search("^Death Block: 'points'",line):
              stuff = line.split()
              pts = int(stuff[-2])
    return pts

def assigned_colors():
  red = rgb2int( (255,0,0) )
  green = rgb2int( (0,255,0) )
  blue = rgb2int( (0,0,255) )
  black = rgb2int( (0,0,0) )
  brown = rgb2int( hex2rgb(brownHex()) )
  grey = rgb2int( hex2rgb(grayHex()) )
  return  [ (black,1),
           (red,2),
           (green,3),
           (blue,4),
           (grey,5),
           (brown,6)
          ]

def add_required_colors_to_pic(picName):
  im = Image.open(picName)
  (pw,ph) = im.size
  pixs = im.load()
  i = ph-1
  my_colors = assigned_colors()
  for color,blkid in my_colors:
      pixs[pw-1,i] = int2rgb(color)
      i = i - 1
  im.save(picName)

def rgb2hex(r,g,b):
    hex = "#{:02x}{:02x}{:02x}".format(r,g,b)
    return hex

def int2rgb(rgbint):
    return (rgbint // 256 // 256 % 256, rgbint // 256 % 256, rgbint % 256)

def hex2rgb(value):
    value = value.lstrip('#')
    lv = len(value)
    return tuple(int(value[i:i + lv // 3], 16) for i in range(0, lv, lv // 3))

def rgb2int( rgb ):
  return 65536 * rgb[0] + 256 * rgb[1] + rgb[2]


def greenHex():
    return rgb2hex(0,255,0)

def whiteHex():
    return rgb2hex(255,255,255)

def grayHex():
    return rgb2hex(128,128,128)

def brownHex():
    return rgb2hex(139,69,19)

def genHexPix(rgbPix,pw,ph):
    hex_pixs = []
    for x in range(pw):
        pixs = []
        for y in range(ph):
            r,g,b = rgbPix[x,y]
            pixs.append(rgb2hex(r,g,b))
        hex_pixs.append(pixs)
    return hex_pixs

def pic_colors( im ):
  (pw,ph) = im.size
  pix = im.load()
  colors = {}
  pixel_id = 1
  hs = list( reversed( range(ph) ) )
  ws = range(pw)
  for i in ws:
    for j in hs:
      pixel_color = rgb2int( pix[i,j] )
      if pixel_color in colors:
        colors[ pixel_color ].append( pixel_id )
      else:
        colors[ pixel_color ] = [ pixel_id ]
      pixel_id =  pixel_id + 1
  return colors

def assign_ids_to_block( blk_id, ids, jouFile):
  strIDs = [ str(x) for x in ids ]
  jouFile.write('block ' + str(blk_id) + ' add hex ' + ' '.join(strIDs) + '\n')

def assign_colors_to_blocks( im, jouFile ):
  my_colors = assigned_colors()
  png_colors = pic_colors( im )
  for color,blid in my_colors:
    if color in png_colors:
      assign_ids_to_block( blid, png_colors[color], jouFile )

def make_mesh_of_pic( pic_name, mesh_name ):
  jouFile = open('make_mesh.jou','w')
  im = Image.open(pic_name)
  (pw,ph) = im.size
  jouFile.write('brick x ' + str(pw) + ' y ' + str(ph) + ' z 1\n')
  jouFile.write('vol all size 1.0\n')
  jouFile.write('mesh vol all\n')
  assign_colors_to_blocks( im, jouFile )
  jouFile.write('export mesh "' + mesh_name + '" overwrite\n')
  jouFile.write('exit\n')
  jouFile.close()
  os.system('cubit -batch -nographics make_mesh.jou')
  
def write_latex(builderName,score):
  f = open('bridgeReport.tex','r')
  latexLines = f.readlines()
  f.close()
  w = open('my_bridgeReport.tex','w')
  builderName = builderName.replace('_',' ')
  for line in latexLines:
    line = line.replace('<BUILDERNAME>',str(builderName))
    line = line.replace('<SCORE>',str(score))
    w.write(line)
  w.close()
  
def make_base_file_name(builderName,score):
  fname = strnum(score) + '_' + builderName.replace(' ','_') + '_' + time.strftime('%Y%m%d-%H%M%S')
  return fname
  
class PaintApp(tkinter.Frame):
    
    def run_full_sim(self):
      os.system('module purge;module load sierra;sierra -j 16 adagio -i bridgeFull.i')

    def run_load_sim(self):
      os.system('module purge;module load sierra;sierra -j 16 adagio -i bridgeLoad.i')

    def post_full_sim(self):
      os.system('module purge;module load viz;pvpython make_paraview_pics.py')
      os.system('module purge;module load viz;/usr/netpub/ffmpeg/bin/ffmpeg -i my_pic_%04d.png my_movie.avi')
      os.system('module purge;module load viz;envideo my_movie.avi')
      # os.system('module purge;module load viz;paraview --state=bridgeFalls_Stressed.pvsm')
      write_latex( self.BuilderName.get(), self.compute_score() )
      os.system('pdflatex my_bridgeReport.tex')
      os.system('evince my_bridgeReport.pdf')
      if not os.path.isdir('saved_results'):
          os.mkdir('saved_results')
      bName = make_base_file_name(self.BuilderName.get(),self.compute_score() )
      rptName = bName + '.pdf'
      picName = bName + '.png'
      os.system('cp my_bridgeReport.pdf ./saved_results/' + rptName )
      os.system('cp my_bridge.png ./saved_results/' + picName )

    def post_load_sim(self):
      os.system('module purge;module load viz;paraview --state=bridgeLoads_Stressed.pvsm')
    
    def compute_score(self):
      pts = get_points()
      pt_frac = pts / 331 ## Full Pts all Cars across
      cost_frac = 1.0 - self.cost / 9313.60 ## 8901.10 ## Less Cost more Pts
      score = int(pt_frac * 10000) + int(cost_frac * 10000)
      return score

    def reduce_image_size(self, imgName, reducedImgName, prevReductionFactor = 1):
      im = Image.open(imgName)
      (pw,ph) = im.size
      pixs = im.load()
      pwn = int(pw / 2)
      phn = int(ph / 2)
      imn = Image.new('RGB',(pwn,phn),(255,255,255))
      pixsn = imn.load()
      xs = range(0,pw,2)
      ys = range(0,ph,2)
      print ("pwn=",pwn,"phn=",phn)
      print ("xs=",len(xs),"ys=",len(ys))
      for x in range(pwn):
        for y in range(phn):
          x1 = xs[x]
          x2 = x1 + 1
          y1 = ys[y]
          y2 = y1 + 1
          ps = [ pixs[x1,y1], pixs[x1,y2], pixs[x2,y1], pixs[x2,y2] ]
          ps = [ rgb2int(p) for p in ps ]
          pm = int2rgb( min(ps) )
          pixsn[x,y] = pm
      imn.save(reducedImgName)

    def line_draw(self,xs,ys,xe,ye,w=1,c_hex=None,c_rgb=(0,0,0),event=None):
      if c_hex == None:
          c_hex = rgb2hex(c_rgb)
      else:
          c_rgb = hex2rgb(c_hex)

      self.draw_save_image.line((xs,ys,xe,ye),c_rgb,width=w)
      if event:
        event.widget.create_line(xs,ys,xe,ye,
                                 fill = c_hex,
                                 capstyle = tkinter.PROJECTING,
                                 width=w, smooth=True)
      else:
        self.drawing_area.create_line(xs,ys,xe,ye,
                                      fill=c_hex,
                                      width=w,capstyle=tkinter.PROJECTING,
                                      smooth=False)

    def redraw_pic(self,drawAll=False):
      white = rgb2hex(255,255,255)
      if drawAll:
        self.canvas_image = self.drawing_area.create_image(0,0,anchor=NW,image=self.image_tk)
        for a in range(self.width):
          for b in range(self.height):
            self.draw_save_image.point((a,b),self.image_pix[a,b])
      else:
        for a in range(self.width):
          for b in range(self.height):
              if self.image_pix_hex[a][b] != white:
                self.line_draw(a,b,a+1,b+1,1,self.image_pix_hex[a][b],self.image_pix[a,b])

    def left_but_down(self, event=None):
        self.left_but = "down"

        #self.x_pos = event.x
        #self.y_pos = event.y

        self.x1_line_pt = event.x
        self.y1_line_pt = event.y

    def left_but_up(self, event=None):
        self.left_but = "up"

        self.x_pos = None
        self.y_pos = None

        self.x2_line_pt = event.x
        self.y2_line_pt = event.y

        # Draw the line
        if self.drawing_tool.get() == "Line":
          self.line_draw(self.x1_line_pt, self.y1_line_pt, self.x2_line_pt, self.y2_line_pt, int(self.pensize), self.penColor.get(), None, event)

        # self.redraw_pic()

    def motion(self, event=None):
        if self.drawing_tool.get() == "Pencil":
            self.pencil_draw(event)

    # ---------- DRAW PENCIL ----------
    def pencil_draw(self, event=None):
        if self.input_mode.get() == "Screen":
            self.line_draw(event.x,event.y,event.x,event.y,int(self.pensize),self.penColor.get(),None,event)
        elif self.left_but == "down":
            # Make sure x and y have a value
            if self.x_pos is not None and self.y_pos is not None:
                self.line_draw(self.x_pos,self.y_pos,event.x,event.y,int(self.pensize),self.penColor.get(),None,event)
            self.x_pos = event.x
            self.y_pos = event.y

    def set_pensize(self,width):
      self.pensize = width

    def delete_app(self):
      self.master.destroy()

    def name_has_been_set(self):
      if (self.BuilderName.get() == "<insert your name>"):
        messagebox.showerror("Error", "Please enter a Builder Name")
        return False
      return True

    def cost_has_been_computed(self):
      if (self.costVar.get() == "$0.00"):
        messagebox.showerror("Error", "Before building, let's see how much this will cost. Press compute cost!")
        return False
      return True

    def compute_cost(self):
      self.redraw_pic()

      grey = hex2rgb(grayHex())
      green = (0,255,0)
      brown = hex2rgb(brownHex())

      numGrey = 0
      numGreen = 0
      numBrown = 0

      for pixel in self.saved_image.getdata():
        if pixel == grey:
          numGrey = numGrey + 1
        elif pixel == green:
          numGreen = numGreen + 1
        elif pixel == brown:
          numBrown = numBrown + 1

      self.cost = round(numGrey*0.05 + numGreen*0.03 + numBrown*0.01, 2)
      costTxt = "$%.2f" % self.cost

      self.costVar.set(costTxt)

    def save_image(self):
      file_cleanup()
      self.compute_cost()
      self.saved_image.save('my_bridge.png')
      self.reduce_image_size('my_bridge.png','my_med_bridge.png')
      self.reduce_image_size('my_med_bridge.png','my_small_bridge.png',2)
      add_required_colors_to_pic('my_small_bridge.png')
      make_mesh_of_pic('my_small_bridge.png','my_bridge.g')
      delete_frags()

    def run_full(self):
      if (not self.name_has_been_set() or not self.cost_has_been_computed()): return
      self.save_image()
      self.run_full_sim()
      self.post_full_sim()

    def run_load(self):
      if (not self.name_has_been_set() or not self.cost_has_been_computed()): return
      self.save_image()
      self.run_load_sim()
      self.post_load_sim()

    def reset_image(self):
      self.redraw_pic(True)
      self.BuilderName.set("<insert your name>")
      self.costVar.set("$0.00")

    def create_widgets(self):
      ## Quit
      self.quit = tkinter.Button(self, text="Quit", bg="red",command=self.delete_app)
      self.quit.pack(side="bottom")

      ## Simulate Full
      self.sim_full = tkinter.Button(self, text="Complete Simulation", command=self.run_full)
      self.sim_full.pack(side="bottom")

      ## Simulate Gravity Load
      self.sim_load = tkinter.Button(self, text="Simulate Gravity Load", command=self.run_load)
      self.sim_load.pack(side="bottom")

      ## Cost
      self.costRow = Frame(self)
      self.costLabel = Label(self.costRow, width=12, text="Bridge Cost", anchor='w')
      self.costVar = StringVar()
      self.costEntry = Entry(self.costRow, width=15, textvariable=self.costVar)
      self.costEntry.bind("<Key>",lambda e: "break")
      self.costButton = tkinter.Button(self.costRow, text="Compute Cost", command=self.compute_cost)
      self.costRow.pack(side=BOTTOM, fill=X)
      self.costButton.pack(side=RIGHT)
      self.costEntry.pack(side=RIGHT, expand=YES, fill=X)
      self.costLabel.pack(side=RIGHT)
      self.costVar.set("$0.00")

      ## Builder Name
      self.entryRow = Frame(self)
      self.nameLabel = Label(self.entryRow, width=12, text="Builder Name", anchor='w')
      self.BuilderName = StringVar()
      self.nameEntry = Entry(self.entryRow, width=15, textvariable=self.BuilderName)
      self.entryRow.pack(side=BOTTOM, fill=X)
      self.nameEntry.pack(side=RIGHT, expand=YES, fill=X)
      self.nameLabel.pack(side=RIGHT)
      self.BuilderName.set("<insert your name>")

      ## Screen mode
      self.input_mode = StringVar()
      self.input_tools = {"Mouse":"Mouse", "Screen":"Touch Screen"}
      for tool, toolTxt in self.input_tools.items():
          b = Radiobutton(self, text=toolTxt, variable=self.input_mode, value=tool, indicatoron=0)
          b.pack(side=RIGHT, anchor=W)
      self.input_mode.set("Mouse")
      self.bufferLabel = Label(self, width=5, text="", anchor='w')
      self.bufferLabel.pack(side=RIGHT)

      ## Drawing tool
      self.drawing_tool = StringVar()
      self.tools = {"Line":"Line Draw", "Pencil":"Free Draw"}
      for tool, toolTxt in self.tools.items():
          b = Radiobutton(self, text=toolTxt, variable=self.drawing_tool, value=tool, indicatoron=0)
          b.pack(side=RIGHT, anchor=W)
      self.drawing_tool.set("Pencil")

      ## Pensize
      self.slider = tkinter.Scale(self, label="Pen Size", from_=12, to=24, command=self.set_pensize, orient=tkinter.HORIZONTAL)
      self.slider.set(5)
      self.slider.pack(side="right")

      ## Pen color corresponding to material
      self.penColor = StringVar()
      self.materials = {"Eraser": whiteHex(),
                        "Wood": brownHex(),
                        "Aluminum": greenHex(),
                        "Steel": grayHex(),}

      for material, color in self.materials.items():
        b = Radiobutton(self, text=material, variable=self.penColor, value=color, indicatoron=0)
        b.pack(side=RIGHT, anchor=W)
      self.penColor.set(grayHex())

      ## Cleanup
      self.reset = tkinter.Button(self, text="Clean Up", command=self.redraw_pic)
      self.reset.pack(side="right")

      ## Reset
      self.reset = tkinter.Button(self, text="Reset", command=self.reset_image)
      self.reset.pack(side="right")

    def __init__(self, master=None):
      super().__init__(master)
      self.pack()
      # Tracks whether left mouse is down
      self.left_but = "up"
      # x and y positions for drawing with pencil
      self.x_pos, self.y_pos = None, None
      # x and y positions when mouse is clicked and released
      self.x1_line_pt, self.y1_line_pt, self.x2_line_pt, self.y2_line_pt= None, None, None, None
      self.pensize = 12
      self.background_image = '3car_bridge.png'
      self.image = Image.open(self.background_image)
      self.width = self.image.size[0]
      self.height = self.image.size[1]
      self.saved_image = Image.new('RGB',(self.width,self.height),(255,255,255))
      self.draw_save_image  = ImageDraw.Draw(self.saved_image)
      self.image_pix = self.image.load()
      self.image_pix_hex = genHexPix(self.image_pix,self.width,self.height)
      self.image_tk = ImageTk.PhotoImage(self.image)
      self.drawing_area = Canvas(root,width=self.width,height=self.height)
      self.drawing_area.pack()
      self.drawing_area.bind("<Motion>", self.motion)
      self.drawing_area.bind("<ButtonPress-1>", self.left_but_down)
      self.drawing_area.bind("<ButtonRelease-1>", self.left_but_up)
      self.create_widgets()
      self.reset_image()


if __name__ == "__main__":
  root = Tk()
  paint_app = PaintApp(root)
  root.mainloop()
  file_cleanup()

