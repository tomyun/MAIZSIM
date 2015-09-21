from .organ import Organ

#FIXME not implemented yet
class Root(Organ):
    def setup(self):
        #FIXME make it general property of all organs
        self.initiated = False
