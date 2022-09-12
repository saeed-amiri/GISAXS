import typing
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


class Sigma:
    """constant values of the sigma
    SIGMA is the size of the solvent particle or the constituent
    particles integrated over in the colloidal particle and should
    typically be set as follows:
    """
    SIGMA_CC: float = 1.0  # colloid-colloid
    SIGMA_SS: float = 1.0  # colloid-solvent,arithmetic mixing colloid&solvent
    SIGMA_CS: float = 1.0  # solvent-solvent or size of  the solvent particle


class Colloid:
    """modeling Colloid interaction between particles
    input:
        A: float,  # energy unit, energy prefactor form Hamaker calss
        sigma: float,  # distance unit, size of the solvent particles
        from Sigma class
        d1: float,  # distance unit, particles diameters
        d2: float,  # distance unit, particles diameters
        r_cut: a flaot or an array of float  # distance units, cutoff
    """
    def __init__(self,
                 d1: float,  # distance unit, particles diameters
                 d2: float,  # distance unit, particles diameters
                 r_cut: typing.Any  # distance units, cutoff, float/array
                 ) -> None:
        self.get_colloid(d1, d2, r_cut)

    def get_colloid(self,
                    d1: float,  # distance unit, particles diameters
                    d2: float,  # distance unit, particles diameters
                    r_cut: typing.Any  # distance units, cutoff, float/array
                    ) -> None:
        """call all the methods to calculate the interaction"""
        a1: float = d1 / 2  # Radius of 1st particle
        a2: float = d2 / 2  # Radius of 2nd particle
        self.U_cc = self.colloid_colloid(a1, a2, r_cut)

    def colloid_colloid(self,
                        a1: float,  # distance unit, particles radius
                        a2: float,  # distance unit, particles radius
                        r_cut: typing.Any  # distance units, cutoff,float/array
                        ) -> typing.Any:  # float or np.array
        """get colloid-colloid interaction"""
        A = Hamaker.A_CC
        sigma = Sigma.SIGMA_CC
        divs1: float  # 1st divisor of the equation in U_A
        divs2: float  # 2nd divisor of the equation in U_A
        a12: float = 2*a1*a2  # for simplification
        divs1 = r_cut**2-(a1+a2)**2
        divs2 = r_cut**2-(a1-a2)**2
        U_a: float = -(A/6)*(a12/divs1 + a12/divs2 + np.log(divs1/divs2))
        divd1: float  # 1st divider of the equation in U_R
        divd2: float  # 2nd divider of the equation in U_R
        r2: float = r_cut*r_cut  # for simplification
        divd1 = 6*(a1**2+7*a1*a2+a2**2)
        divd2 = 6*(a1**2-7*a1*a2+a2**2)
        U_r = (A/37800)*(sigma**6)*(
            (r2-7*r_cut*(a1+a2)+divd1) / (r_cut-a1-a2)**7 +
            (r2+7*r_cut*(a1+a2)+divd1) / (r_cut+a1+a2)**7 -
            (r2+7*r_cut*(a1-a2)+divd2) / (r_cut+a1-a2)**7 -
            (r2-7*r_cut*(a1-a2)+divd2) / (r_cut-a1+a2)**7)
        return U_a + U_r

    def colloid_solvent(self) -> None:
        """get colloid_solvent interaction"""

    def solvent_solvent(self) -> None:
        """get solvent_solvent interaction"""


if __name__ == '__main__':
    r = [i*2/10 for i in range(10)]
    r = np.array(r)
    coll = Colloid(1, 10, r)
    print(coll.U_cc)
