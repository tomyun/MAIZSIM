from .trait import Trait

class Ratio(Trait):
    @property
    def carbon_to_mass(self):
        # Weight.C_to_CH2O_ratio
        # 40% C, See Kim et al. (2007) EEB
        return 0.40

    @property
    def shoot_to_root(self):
        return 0.7

    @property
    def root_to_shoot(self):
        return 1 - self.shoot_to_root

    @property
    def leaf_to_stem(self):
        return 0.9

    @property
    def stem_to_leaf(self):
        return 1 - self.leaf_to_stem

    @property
    def initial_leaf(self):
        #TODO how to handle primordia?
        return self.shoot_to_root * self.leaf_to_stem / self.p.primordia
