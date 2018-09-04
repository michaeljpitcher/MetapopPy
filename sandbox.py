from metapoppy_group_by_patch import *

class NAEvent(Event):
    REACTION_PARAMETER = 1.1
    def __init__(self):
        Event.__init__(self)

a = NAEvent()