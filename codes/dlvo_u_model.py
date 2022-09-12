import numpy as np


class Doc:
    """scripting colloid interaction in the LAMMPS colloid piar_style
    and yukawa/colloid.
    The colloid equation is a combination of three interactions:
        - colloid-colloid interaction energy
        - colloid-solvent interaction
        - solvent-solvent
    And the yukawa/colloid:
        - coulombic interaction between two colloid particles,
        screened due to the presence of an electrolyte
    This script simulates all the above interactions.
    Changing the values of variables shows behavior change in the main
    interaction.
    """


class Hamaker:
    """Hamaker energy prefactor"""
    A_CC: float = 39.5  # 4*np.pi**2 = 39.478417 colloid-colloid
    A_SS: float = 144.0  # solvent-solvent (if epsilon=1, so that 144/36=4)
    A_CS: float = 5688.0  # np.sqrt(A_CC*A_SS)  colloid-solvent


class Colloid:
    """modeling Colloid interaction between particles"""
    def __init__(self) -> None:
        pass

    def get_colloid(self) -> None:
        """call all the methods to calculate the interaction"""
        pass

    def colloid_colloid(self) -> None:
        """get colloid-colloid interaction"""

    def colloid_solvent(self) -> None:
        """get colloid_solvent interaction"""

    def solvent_solvent(self) -> None:
        """get solvent_solvent interaction"""
