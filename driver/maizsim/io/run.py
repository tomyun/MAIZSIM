from .base import LegacyFile

class Run(LegacyFile):
    @property
    def specs(self):
        return [
            [('weather', str)],
            [('time', str)],
            [('biology', str)],
            [('climate', str)],
            [('nitrogen', str)],
            [('solute', str)],
            [('soil', str)],
            [('management', str)],
            [('water', str)],
            [('water_boundary', str)],
            [('initials', str)],
            [('variety', str)],
            [('geometry', str)],
            [('node_geometry', str)],
            [('element_geometry', str)],
            [('mass_balance', str)],
            [('plant_graphics', str)],
            [('leaf_graphics', str)],
            [('node_graphics', str)],
            [('element_graphics', str)],
            [('surface_graphics', str)],
            [('flux_graphics', str)],
            [('mass_balance_out', str)],
        ]
