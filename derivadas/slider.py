import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button

class ConfigGUI:

    def __init__(self):
        self.setup()

    def setup(self):
        self.fig, self.ax = plt.subplots()
        self.t = np.linspace(0, 10, 1000)
        self.line, = plt.plot(self.t, np.sin(self.t), lw=2)
        self.v = 0
        slider_ax = plt.axes([0.1, 0.1, 0.8, 0.02])
        slider = Slider(slider_ax, "Offset", -5, 5, valinit=0, color='#AAAAAA')
        slider.on_changed(self.on_change)
        ax1 = plt.axes([0.2, 0.5, 0.1, 0.075])
        #ax2 = plt.axes([0.7, 0.5, 0.1, 0.075])

        b1 = Button(ax1, 'Draw')
        b1.on_clicked(self.callback)
        plt.show()

    def on_change(self, val):
        self.v = val
        print(str(self.v))

    def callback(self, event):
        print("clicked:"+str(event))
        print("v="+str(self.v))
    #    sys.stdout.flush()
        self.line.set_ydata(np.sin(self.t - self.v))
        self.fig.canvas.draw()

S()
