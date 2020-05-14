#!python

import matplotlib.backend_tools


class PointerTool(matplotlib.backend_tools.ToolToggleBase):
    description = 'Select nodes with left mouse, deselect with right, scroll to undo/redo'
    image = "../lib/browser_images/pointer_25x25.png"
    radio_group = 'default'
    default_keymap = 's'

    def __init__(self, *args):
        super().__init__(*args)
        self._pick_connection = None
        self._zoom_connection = None
        self._stack = []
        self._stack_pointer = -1
        self._current_selection = set()

    def enable(self, event):
        self._pick_connection = self.figure.canvas.mpl_connect(
            'pick_event',
            self.pick_event
        )
        self._zoom_connection = self.figure.canvas.mpl_connect(
            'scroll_event',
            self.scroll_event
        )

    def disable(self, event):
        self.figure.canvas.mpl_disconnect(self._zoom_connection)
        self.figure.canvas.mpl_disconnect(self._pick_connection)

    def scroll_event(self, event):
        if event.button == "up":
            self.update_stack_pointer(direction=1)
        elif event.button == "down":
            self.update_stack_pointer(direction=-1)
        self.custom_event()

    def custom_event(self):
        scatter = self.figure.axes[0].collections[0]
        edge_colors = scatter.get_facecolors()
        try:
            edge_colors[list(self._current_selection)] = [0,0,0,1]
        except IndexError:
            edge_colors = np.repeat(edge_colors, scatter._offsets.shape[0], axis=0)
            edge_colors[list(self._current_selection)] = [0,0,0,1]
        scatter.set_edgecolors(edge_colors)
        plt.draw()
        # print(self._stack, self._stack_pointer, self._current_selection)

    def pick_event(self, event):
        for index in event.ind:
            if event.mouseevent.button == 1:
                self.update_stack(index + 1)
            if event.mouseevent.button == 3:
                self.update_stack(-index - 1)
        self.custom_event()

    def update_stack(self, index):
        if (index > 0) and (index - 1 in self._current_selection):
            return
        if (index < 0) and (-index - 1 not in self._current_selection):
            return
        self._stack = self._stack[:self._stack_pointer + 1]
        self._stack.append(index)
        self._stack_pointer += 1
        self.update_current_selection(index)

    def update_current_selection(self, index):
        if index > 0:
            self._current_selection.add(index - 1)
        elif index < 0:
            self._current_selection.remove(-index - 1)

    def update_stack_pointer(self, direction=None, target=None):
        if target is None:
            target = self._stack_pointer + direction
            target = max(0, target)
            target = min(len(self._stack) - 1, target)
        direction = 1 if target > self._stack_pointer else -1
        for i in range(abs(target - self._stack_pointer)):
            index = self._stack[self._stack_pointer]
            self.update_current_selection(direction * index)
            self._stack_pointer += direction

if __name__ == "__main__":
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib
    matplotlib.use('TkAgg')
    plt.rcParams['toolbar'] = 'toolmanager'
    fig = plt.figure("Ion-network browser")
    fig.add_subplot(111)
    xs, ys = np.random.random(10), np.random.random(10)
    fig.axes[0].scatter(xs, ys, picker=True)
    fig.canvas.manager.toolmanager.add_tool('Select nodes', PointerTool)
    fig.canvas.manager.toolbar.add_tool('Select nodes', 'navigation', 1)
    plt.show()
else:
    pass

#
# from sandbox_browser import PointerTool
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib
# # matplotlib.use('TkAgg')
# plt.rcParams['toolbar'] = 'toolmanager'
# fig = plt.figure("Ion-network browser")
# fig.add_subplot(111)
# xs, ys = np.random.random(10), np.random.random(10)
# s = fig.axes[0].scatter(xs, ys, picker=True)
# fig.canvas.manager.toolmanager.add_tool('Select nodes', PointerTool)
# fig.canvas.manager.toolbar.add_tool('Select nodes', 'navigation', 1)
# plt.show(block=False)
# tool = fig.canvas.manager.toolmanager._tools['Select nodes']
# def bind(instance, func, as_name=None):
#     """
#     Bind the function *func* to *instance*, with either provided name *as_name*
#     or the existing name of *func*. The provided *func* should accept the
#     instance as the first argument, i.e. "self".
#     """
#     if as_name is None:
#         as_name = func.__name__
#     bound_method = func.__get__(instance, instance.__class__)
#     setattr(instance, as_name, bound_method)
#     return bound_method
#
# def test(self):
#     edges = ["None"] * s._offsets.shape[0]
#     for i in self._current_selection:
#         edges[i] = "black"
#     s.set_edgecolors(edges)
#     print("test", self._stack, self._stack_pointer, self._current_selection)
#     plt.draw()
#
# tool.custom_event = bind(tool, test)
