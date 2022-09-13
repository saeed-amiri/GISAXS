import typing
import numpy as np
from colors_text import TextColor as bcolors
import matplotlib.pylab as plt


class Doc:
    """scripting colloid interaction in the LAMMPS colloid piar_style
    and yukawa/colloid.
    https://docs.lammps.org/pair_colloid.html
    "Style colloid computes pairwise interactions between large
    colloidal particles and small solvent particles using 3 formulas.
    A colloidal particle has a size > sigma; a solvent particle is the
    usual Lennard-Jones particle of size sigma."
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

    "d1 and d2 are particle diameters, so that d1 = 2*a1 and d2 = 2*a2
    in the formulas above. Both d1 and d2 must be values >= 0.
    If d1 > 0 and d2 > 0, then the pair interacts via the colloid-colloid
    formula above.
    If d1 = 0 and d2 = 0, then the pair interacts via the solvent-solvent
    formula. I.e. a d value of 0 is a Lennard-Jones particle of size sigma.
    If either d1 = 0 or d2 = 0 and the other is larger, then the pair
    interacts via the colloid-solvent formula."
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
        a: float  # particle size
        if a1 > 0 and a2 > 0:
            self.U = self.colloid_colloid(a1, a2, r_cut)
        elif (a1 == 0 or a2 == 0) and not (a1 == 0 and a2 == 0):
            if a1 > 0:
                a = a1
            elif a2 > 0:
                a = a2
            self.U = self.colloid_solvent(a, r_cut)
        elif a1 == 0 and a2 == 0:
            self.U = self.solvent_solvent(r_cut)
        else:
            exit(f'{bcolors.FAIL}{self.__class__.__name__}:\n'
                 f'Wrong values for diameters{bcolors.ENDC}\n'
                 f'{bcolors.OKGREEN}{Doc.__doc__}{bcolors.ENDC}\n')

    def colloid_colloid(self,
                        a1: float,  # distance unit, particles radius
                        a2: float,  # distance unit, particles radius
                        r_cut: typing.Any  # distance units, cutoff,float/array
                        ) -> typing.Any:  # float or np.array
        """get colloid-colloid interaction
        "This equation results from describing each colloidal particle
        as an integrated collection of Lennard-Jones particles of size
        sigma and is derived in (Everaers)."
        """
        print(f'{bcolors.OKCYAN}Colloid-Colloid interaction:{bcolors.ENDC}\n'
              f'\t{bcolors.OKBLUE}d1={2*a1}, d2={2*a2}{bcolors.ENDC}')
        A: float = Hamaker.A_CC
        sigma: float = Sigma.SIGMA_CC
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
        print(f'{bcolors.OKCYAN}\tmin: U_R={np.min(np.abs(U_r)):.3e}, '
              f'U_A={np.min(np.abs(U_a)):.3e}\n'
              f'\tmax: U_R={np.max(np.abs(U_r)):.3e}, '
              f'U_A={np.max(np.abs(U_a)):.3e}\n{bcolors.ENDC}')
        return U_a + U_r

    def colloid_solvent(self,
                        a: float,  # distance unit, colloidal particles radius
                        r_cut: typing.Any  # distance units, cutoff,float/array
                        ) -> typing.Any:  # float or np.array
        """get colloid_solvent interaction
        "This formula is derived from the colloid-colloid interaction,
        letting one of the particle sizes go to zero."
        """
        print(f'{bcolors.OKCYAN}Colloid-Solvent interaction{bcolors.ENDC}\n')
        A: float = Hamaker.A_CS
        sigma: float = Sigma.SIGMA_CS
        # For simplification:
        r2: float = r_cut*r_cut
        r4: float = r2*r2
        r6: float = r2*r4
        a2: float = a*a
        a3: float = a*a2
        a4: float = a*a3
        a6: float = a2*a4
        s2: float = sigma*sigma
        s3: float = sigma*s2
        s6: float = s3*s3
        U = (2*a3*s3*A) / (9*(a2-r2)**3)*(
            1-(
                (5*a6+45*a4*r2+63*a2*r4+15*r6)*s6 /
                (15*(a-r_cut)**6*(a+r_cut)**6)
            )
        )
        return U

    def solvent_solvent(self,
                        r_cut: typing.Any  # distance units, cutoff,float/array
                        ) -> typing.Any:
        """get solvent_solvent interaction
        The solvent-solvent interaction energy is given by the usual
        Lennard-Jones formula
        """
        print(f'{bcolors.OKCYAN}Solvent-Solvent interaction{bcolors.ENDC}\n')
        A: float = Hamaker.A_SS
        sigma: float = Sigma.SIGMA_SS
        U: typing.Any = (A/36)*(
            (sigma/r_cut)**12-(sigma/r_cut)**6
        )
        return U


if __name__ == '__main__':
    for d2 in range(1, 2):
        # d2=0
        for d1 in range(d2+1, 10):
            r = [i/10 for i in range((d1+d2+1)*10, 260)]
            r = np.array(r)
            coll = Colloid(d1=d1, d2=d2, r_cut=r)
            u = np.nan_to_num(coll.U)
            plt.plot(r, u, label=f'd1={d1}, d2={d2}')
            plt.legend()
    plt.show()
