"""
Python module to optimize the interaction cutoff in a 
CP2k simulation.

Notes
-----
Currently the auto feature of the optimizer does not work and so should not be called.
"""

import numpy as np
import matplotlib.pyplot as plt
import ase
from ase.calculators.calculator import Calculator
from ase.calculators.cp2k import CP2K
from typing import Union
import sys


class CutoffOptimizer:
    """
    Class for the optimization of the integration grid cutoff for a cp2k simulation.

    Attributes
    ----------
    calculator : Calculator
            The ase calculator object to be added to the configuration and then used for analysis.
    loop_range : np.array
            cutoff values to loop over and test. All are integers.
    energy_array : list
            list in which the energy values at each cutoff are stored for later or on-the-fly analysis.
    force_array : list
            list in which the forces are stored for later or on-the-fly analysis.
    start : float
                Starting cutoff value to use in the analysis
    stop : Union[float, str]
            Stop value to use in the analysis. If a float, the cutoff values will be taken between start and top, if
            the string auto is given, an automatic optimization method is implemented to find the optimal cutoff.
    calc_parameters : dict
            Parameters to be used by ase to run the calculation. As these are numerous and varied, we do not ask for
            them individually. The parameters should be given as a dict with the name as the key and the value of
            the name as the key value, e.g.
                {'xc': 'PBE', 'basis_set': 'DZVP-MOLOPT-SR-GTH', 'command': 'cp2k.popt'}
            Note that you need not add the cutoff parameters here as these will be set.
    atoms : ase.Atoms
            Ase atoms object on which the analysis should be performed. We ask for an ase atoms object over a
            coordinate file in order to avoid issues arising from file readers.
    tolerance : float
            The delta value at which the cutoff should be considered converged.
    """

    def __init__(self, start: float = 100, stop: Union[float, str] = 1600, atoms: ase.Atoms = None,
                 calc_parameters: dict = None, tolerance: float = 0.002):
        """
        Constructor for the cutoff optimizer class

        Parameters
        ----------
        start : float
                Starting cutoff value to use in the analysis
        stop : Union[float, str]
                Stop value to use in the analysis. If a float, the cutoff values will be taken between start and top, if
                the string auto is given, an automatic optimization method is implemented to find the optimal cutoff.
        calc_parameters : dict
                Parameters to be used by ase to run the calculation. As these are numerous and varied, we do not ask for
                them individually. The parameters should be given as a dict with the name as the key and the value of
                the name as the key value, e.g.
                    {'xc': 'PBE', 'basis_set': 'DZVP-MOLOPT-SR-GTH', 'command': 'cp2k.popt'}
                Note that you need not add the cutoff parameters here as these will be set.
        atoms : ase.Atoms
                Ase atoms object on which the analysis should be performed. We ask for an ase atoms object over a
                coordinate file in order to avoid issues arising from file readers.
        tolerance : float
                The delta value at which the cutoff should be considered converged.
        """
        self.start = start
        self.stop = stop
        self.atoms = atoms
        self.calc_parameters = calc_parameters
        self.tolerance = tolerance

        self.calculator: Calculator
        self.loop_range: list
        self.force_array: list = []
        self.energy_array: list = []

        # Run checks
        if type(stop) is str:
            self._temp_operation_check()
        if calc_parameters is None:
            print("You must give adequate calculator parameters")
            sys.exit(1)

    @staticmethod
    def _temp_operation_check():
        """
        If the auto option is called, this method will be called and the process ended.
        Returns
        -------
        Prints a message to the screen end exits the code.
        """
        print("Sorry, the 'auto' implementation for optimizing cutoff range is currently not supported. Please check"
              "back in soon for updates. In the mean time, we suggest setting a cutoff maximum of 1500 as the code can"
              "tell you later if a greater value is required.")

        sys.exit(1)

    def _get_static_range(self):
        """
        Get the loop range for a system not using the auto method.
        Returns
        -------

        """
        self.loop_range = np.linspace(self.start, self.stop, 10, dtype=int)

    def _build_ase_calculator(self):
        """
        Build an ase calculator object
        Returns
        -------

        """
        self.calculator = CP2K(**self.calc_parameters)
        self.calculator.max_scf = 1  # set this to 1 for efficiency

    def _run_scf_minimization(self):
        """
        Run the SCF minimization to get energy and forces.
        Returns
        -------

        """
        self.energy_array.append(self.atoms.get_potential_energies())
        self.force_array.append(self._reduce_forces(self.atoms.get_forces()))

    def _perform_loop_operation(self, item: int):
        """
        Run energy/force calculation on atoms. Separated for parallelization.

        Parameters
        ----------
        item : int
                Cutoff value to be used in the calculation

        Returns
        -------
        Updates the class state.
        """
        self.calculator.cutoff = item  # set the cutoff
        self.atoms.calc = self.calculator  # attach the calculator
        self._run_scf_minimization()  # run the calculation

    @staticmethod
    def _reduce_forces(forces: list):
        """
        Reduce the n_atoms x 3 output of forces to a single

        Parameters
        ----------
        forces : list
                n_atoms x 3 array of forces
        Returns
        -------

        """
        return np.mean(np.linalg.norm(forces, axis=1), axis=0)

    def _check_progress(self):
        """
        Check the progress of the analysis at run time and end it if convergence is reached
        Returns
        -------

        """
        if len(self.energy_array) < 2:
            return 0
        else:
            en_diff = np.diff(self.energy_array)
            fr_diff = np.diff(self.force_array)

            if en_diff[-1] and fr_diff[-1] < self.tolerance:
                return 0
            else:
                return 1

    def _return_results(self):
        """
        Return the results of the analysis.

        Returns
        -------
        Prints output to screen and saves images.
        """
        print(f"Energy cutoff: {self.energy_array}")
        print(f"Force cutoff: {self.force_array}")

        plt.plot(self.loop_range[:len(self.energy_array)], self.energy_array, 'o-', label="Energy convergence")
        plt.plot(self.loop_range[:len(self.force_array)], self.force_array, 'o-', label="Force convergence")
        plt.xlabel("Cutoff value (Ry)")
        plt.ylabel("Energy and Force Values")
        plt.grid()
        plt.legend()
        plt.savefig("Convergence_Results.pdf", dpi=800)

    def run_optimizer(self):
        """
        Run the optimizer
        Returns
        -------
        Plots images of the convergence and prints the ideal values.

        Notes
        -----
        TODO - Parallelize loop
        """
        self._get_static_range()  # get the correct loop range

        # Run all calculations: TODO - Parallelize.
        for item in self.loop_range:
            self._perform_loop_operation(item)
            if self._check_progress() == 1:
                continue
            else:
                break


