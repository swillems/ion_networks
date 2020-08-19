import ipywidgets as widgets
from ipywidgets import HBox, VBox
import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display

range_slider = widgets.FloatRangeSlider(
    value=[-1., +1.],
    min=-5., max=+5., step=0.1,
    description='xlim:',
    readout_format='.1f',
)
color_buttons = widgets.ToggleButtons(
    options=['blue', 'red', 'green'],
    description='Color:',
)
color_picker = widgets.ColorPicker(
    concise=True,
    description='Background color:',
    value='#efefef',
)

def plot2(b=None):
    xlim = range_slider.value
    freq = 2
    color = color_buttons.value
    t = np.linspace(xlim[0], xlim[1], 1000)
    f, ax = plt.subplots(1, 1, figsize=(8, 6))
    ax.plot(t, np.sin(2 * np.pi * freq * t),
            color=color)


tab1 = VBox(children=[range_slider,
                      ])
# tab2 = VBox(children=[color_buttons,
#                       HBox(children=[title_textbox,
#                                      color_picker,
#                                      grid_button]),
#                                      ])

out = widgets.Output()
tab = widgets.Tab(children=[tab1])
tab.set_title(0, 'plot')
tab.set_title(1, 'styling')
VBox(children=[tab, out])
plt.show()
display()
