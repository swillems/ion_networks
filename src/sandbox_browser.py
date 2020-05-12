#!python

# builtin
import contextlib
# external
import numpy as np
import matplotlib.pyplot as plt
import PySimpleGUI as sg
import matplotlib
# local
import ms_run_files

plt.rcParams['toolbar'] = 'toolmanager'
# matplotlib.use('TkAgg')


class NetworkTool(matplotlib.backend_tools.ToolBase):
    '''List all the tools controlled by the `ToolManager`'''
    # keyboard shortcut
    # default_keymap = 'I'
    description = 'Select ion-network'
    image = "../lib/browser_images/network.png"

    def trigger(self, *args, **kwargs):
        print(self.figure)
        self.figure.select_network()


# class PointerTool(matplotlib.backend_tools.ToolToggleBase):
#     """Base class for `ToolZoom` and `ToolPan`"""
#     description = 'Select nodes'
#     image = "../lib/browser_images/pointer.png"
#     radio_group = 'default'
#
#     def __init__(self, *args):
#         super().__init__(*args)
#         self._button_pressed = None
#         self._xypress = None
#         self._idPress = None
#         self._idRelease = None
#         self._idScroll = None
#
#     def enable(self, event):
#         """Connect press/release events and lock the canvas"""
#         self.figure.canvas.widgetlock(self)
#         self._idPress = self.figure.canvas.mpl_connect('pick_event', onpick1)
#
#     def disable(self, event):
#         """Release the canvas and disconnect press/release events"""
#         self._cancel_action()
#         self.figure.canvas.widgetlock.release(self)
#         self.figure.canvas.mpl_disconnect(self._idPress)
#
#     def onpick1(event):
#         print(event.ind)

class PointerTool(matplotlib.backend_tools.ToolToggleBase):
    """Base class for `ToolZoom` and `ToolPan`"""
    description = 'Select nodes'
    image = "../lib/browser_images/pointer.png"
    radio_group = 'default'

    def __init__(self, *args):
        super().__init__(*args)
        self._idSelect = None

    def enable(self, event):
        # super().enable(event)
        self.figure.canvas.widgetlock(self)
        self._idSelect = self.figure.canvas.mpl_connect(
            'pick_event',
            self.select
        )

    def disable(self, event):
        # super().disable(event)
        self._cancel_action()
        self.figure.canvas.widgetlock.release(self)
        self.figure.canvas.mpl_disconnect(self._idSelect)

    def select(event):
        print(event.ind)


def select_network(self, *args, **kwargs):
    file_name = sg.popup_get_file(
        'Please select an ion-network',
        file_types=(('Ion-network', '*.inet.hdf'),)
    )
    if file_name is None:
        return
    try:
        ion_network = ms_run_files.Network(file_name)
    except (OSError, ValueError):
        sg.popup_error('This is not a valid ion_network')
        return
    try:
        evidence = ms_run_files.Evidence(file_name)
    except (OSError, ValueError):
        sg.popup_error('This ion_network has no valid evidence')
        return
    if hasattr(self, "ion_network") and (self.ion_network == ion_network):
        pass
    else:
        self.ion_network = ion_network
        self.evidence = evidence


matplotlib.pyplot.Figure.select_network = select_network
fig = plt.figure("Ion-network browser")

fig.add_subplot(111)
fig.axes[0].scatter(np.random.random(10), np.random.random(10), picker=True)

# for tool in list(fig.canvas.manager.toolmanager.tools):
#     print(tool)
#     fig.canvas.manager.toolmanager.remove_tool(tool)

fig.canvas.manager.toolmanager.remove_tool('forward')
fig.canvas.manager.toolmanager.remove_tool('back')
fig.canvas.manager.toolmanager.remove_tool('home')

fig.canvas.manager.toolmanager.add_tool('Select nodes', PointerTool)
fig.canvas.manager.toolbar.add_tool('Select nodes', 'navigation', 1)

fig.canvas.manager.toolmanager.add_tool('Select ion-network', NetworkTool)
fig.canvas.manager.toolbar.add_tool('Select ion-network', 'navigation', 1)

plt.show()
