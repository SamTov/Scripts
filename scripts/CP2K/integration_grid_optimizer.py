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
from pycp2k import CP2K
from typing import Union, TextIO
import sys
import glob
import re

class IntegrationGridOptimizer:
    """
    Class for the optimization of the integration grid for a cp2k simulation.

    Attributes
    ----------
    calculator : Calculator
            The ase calculator object to be added to the configuration and then used for analysis.
    loop_range : dict
            cutoff values to loop over and test for the cutoff, ngrid, and rel cutoff e.g.
            {'cutoff': np.linspace(150, 500, 10), 'ngrid': np.linspace(3, 10, 5, dtype=int), ...}
    energy_array : list
            list in which the energy values at each cutoff are stored for later or on-the-fly analysis.
    force_array : list
            list in which the forces are stored for later or on-the-fly analysis.
    start : float
                Starting cutoff value to use in the analysis
    stop : Union[float, str]
            Stop value to use in the analysis. If a float, the cutoff values will be taken between start and top, if
            the string auto is given, an automatic optimization method is implemented to find the optimal cutoff.
    template_file : str
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
    m_grid : CP2K().CP2K_INPUT.FORCE_EVAL_add().DFT.MGRID
            MGRID updated object for easier updating of properties.
    """

    def __init__(self, start=None, stop: Union[dict, str] = None, atoms: ase.Atoms = None, template_file: str = None,
            tolerance: float = 0.002, cwd: str = './', project_name: str = 'optimization'):
        """
        Constructor for the cutoff optimizer class

        Parameters
        ----------
        start : dict
                Starting values to use in the analysis for the cutoff, rel_cutoff and ngrid optimization
        stop : Union[dict, str]
                Stop value to use in the analysis. If a float, the cutoff values will be taken between start and top, if
                the string auto is given, an automatic optimization method is implemented to find the optimal cutoff.
        template_file : str
                A cp2k input file to be read in as a template and modified during the analysis
                Note that you need not add the cutoff, rel cutoff, and monkhorst parameters need not be set here as they
                will be added at the appropriate stage of the analysis.
        atoms : ase.Atoms
                Ase atoms object on which the analysis should be performed. We ask for an ase atoms object over a
                coordinate file in order to avoid issues arising from file readers. If this is left empty, it is
                assumed that the cell is defined in the cp2k input file
        tolerance : float
                The delta value at which the cutoff should be considered converged.
        cwd : str
                Directory in which the analysis should run and files should be stored.
        project_name : str
                Name of the project.
        """
        if start is None:
            start = {'Cutoff': 100, 'Ngrids': 5, 'Rel_cutoff': 60}
        if stop is None:
            stop = {'Cutoff': 1500, 'Ngrids': 8, 'Rel_cutoff': 120}
        self.start = start
        self.stop = stop
        self.atoms = atoms
        self.template_file = template_file
        self.tolerance = tolerance
        self.cwd = cwd
        self.project_name = project_name

        self.calculator = CP2K()
        self.calculator.project_name = self.project_name
        self.calculator.working_directory = self.cwd
        self.calculator.parse(self.template_file)
        self.force_eval = self.calculator.CP2K_INPUT.FORCE_EVAL_list[0]
        self.force_eval.PRINT.FORCES.Section_parameters = "ON"
        self.m_grid = self.force_eval.DFT.MGRID
        self.loop_range = {}
        self.force_array = {}
        self.energy_array = {}
        self.optimized_cutoff: float
        self.optimized_rel_cutoff: float
        self.optimized_n_grids: float
        
        # Run checks
        if type(stop) is str:
            self._temp_operation_check()
        if template_file is None:
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

    def _build_cell(self):
        """
        Add the atoms to the cell if required
        Returns
        -------

        """
        if self.atoms is not None:
            sub_sys = self.force_eval.SUBSYS
            self.calculator.create_cell(sub_sys, self.atoms)
            self.calculator.create_coord(sub_sys, self.atoms)

    def _get_static_range(self):
        """
        Get the loop range for a system not using the auto method. This is currently the only valid method.
        Returns
        -------
        Updates the class
        """
        self.loop_range['Cutoff'] = np.linspace(self.start['Cutoff'], self.stop['Cutoff'], 10, dtype=int)
        self.loop_range['Rel_cutoff'] = np.linspace(self.start['Rel_cutoff'], self.stop['Rel_cutoff'], 10, dtype=int)
        self.loop_range['Ngrids'] = np.linspace(self.start['Ngrids'], self.stop['Ngrids'], 4, dtype=int)

    def _get_force_file(self):
        """
        search for the forces file in the current directory

        Returns
        -------
        force_file : str
                Returns the path to the forces file so that it can be read. If none found, will print an error and end
                the program
        """
        force_file = f"{self.project_name}.out"
        files = glob.glob(f"{self.cwd}/*.xyz")
        if force_file not in files:
            print("There is no force file to be analyzed, something is likely wrong with your cp2k input")
            sys.exit(1)

        return force_file

    @staticmethod
    def _read_energy(f_object: TextIO):
        """
        Read the energy value from a force file.

        Parameters
        ----------
        f_object : TextIO
                File object from which to read the values
        Returns
        -------
        Energy value for that configuration : float
        """
        energy = None
        pattern = "Total FORCE_EVAL"
        for line in f_object:
            if re.search(pattern, line):
                energy = float(line.split()[-1])

        if energy is None:
            print("No energy value found - Are you sure the input file was correct? Check the optimization.out file "
                  "for more details.")
            sys.exit(1)
        else:
            f_object.seek(0)  # return to the start of the file.
            return energy

    @staticmethod
    def _read_forces(f_object: TextIO, n_atoms: int):
        """
        Read in the forces from the text file.

        Parameters
        ----------
        f_object : TextIO
                file object from which to read data
        n_atoms : int
                Number of atoms in the file.

        Returns
        -------
        mean force value as a float : float
        """
        force = 0.0
        pattern_1 = "# Atom   Kind   Element          X              Y              Z"
        pattern_2 = "SUM OF ATOMIC FORCES"

        for i, line, in enumerate(f_object):
            if re.search(pattern_1, line):
                start = i + 1
            elif re.search(pattern_2, line):
                stop = i - 1
        f_object.seek(0)
        for i, line in enumerate(f_object):
            if i < start:
                continue
            elif start <= i <= stop:
                force += np.linalg.norm(np.array(line.split())[3:].astype(float))
            else:
                break

        return force/n_atoms

    def _fetch_properties(self):
        """
        Collect the forces and the global energy from the dumped frc file

        Returns
        -------
        force and energy values : tuple
        """
        force_file = self._get_force_file()
        f_object = open(force_file, 'r')
        n_atoms = int(f_object.readline().split()[0])  # get the number of atoms
        energy = self._read_energy(f_object=f_object)
        force = self._read_forces(f_object=f_object, n_atoms=n_atoms)
        f_object.close()

        return energy, force

    def _update_property(self, updater: dict):
        """
        Update the value of a specific property in the cp2k calculator

        Parameters
        ----------
        updater : dict
                Properties to update as key/value pairs, e.g.
                {'Cutoff': 500, 'Ngrids': 6, 'Rel_cutoff': 60}

        Returns
        -------
        Updates the state of the self.calculator object
        """
        for item in updater:
            print(item)
            self.m_grid.__dict__[item] = updater[item]

    def _perform_energy_force_calculation(self):
        """
        Perform a calculation on a set of atoms.

        Returns
        -------
        updates the class state
        """
        self.calculator.write_input_file()
        self.calculator.run()

    def _update_property_arrays(self, optimized_parameter: str):
        """
        Update the force and energy arrays after a calculation

        Parameters
        ----------
        optimized_parameter : str
                Name of the property to be updated

        Returns
        -------
        Updates the class state
        """
        energy, force = self._fetch_properties()  # collect the energy and force values
        self.energy_array[optimized_parameter].append(energy)
        self.force_array[optimized_parameter].append(force)

    def _perform_loop_operation(self, optimized_parameter: dict):
        """
        Run energy/force calculation on atoms. Separated for parallelization.

        Parameters
        ----------
        optimized_parameter : dict
                A dictionary with the property to be optimized as well as the new value for this property.
                e.g. {'Cutoff': 550}

        Returns
        -------
        Updates the class state.
        """
        self._update_property(optimized_parameter)
        self._perform_energy_force_calculation()
        self._update_property_arrays(list(optimized_parameter.keys())[0])

    def _check_progress(self, optimized_parameter: str):
        """
        Check the progress of the analysis at run time and end it if convergence is reached

        Parameters
        ----------
        optimized_parameter : str
                Property being optimized and therefore updated

        Returns
        -------
        outcome: int
                If 1, continue with the analysis, if 0, the analysis is complete and no more computations should be
                performed
        Will also update the class state in the event of convergence being reached.
        """

        if len(self.energy_array[optimized_parameter]) < 2:
            return 0
        else:
            en_diff = np.diff(self.energy_array[optimized_parameter])
            fr_diff = np.diff(self.force_array[optimized_parameter])

            if en_diff[-1] and fr_diff[-1] < self.tolerance:
                return 0
            else:
                return 1

    def _load_input_file(self):
        """
        Load the cp2k input file and prepare some variables for calls.

        Returns
        -------
        updates the class state
        """
        self.calculator.parse(self.template_file)

    def _run_cutoff_optimization(self):
        """
        Run the cutoff optimization
        """
        for item in self.loop_range['Cutoff']:
            self._perform_loop_operation({'Cutoff': item})
            if self._check_progress('Cutoff') == 1:
                continue
            else:
                self._update_property({'Cutoff': self.energy_array['Cutoff'][-2]})
                break

    def _run_rel_cutoff_optimization(self):
        """
        Run the rel cutoff optimization.

        Returns
        -------

        """
        for item in self.loop_range['Rel_cutoff']:
            self._perform_loop_operation({'Rel_cutoff': item})
            if self._check_progress('Rel_cutoff') == 1:
                continue
            else:
                self._update_property({'Rel_cutoff': self.energy_array['Rel_cutoff'][-2]})
                break

    def _run_n_grids_optimization(self):
        """
        Run the ngrid optimization

        Returns
        -------

        """
        for item in self.loop_range['Ngrids']:
            self._perform_loop_operation({'Ngrids': item})
            if self._check_progress('nNrids') == 1:
                continue
            else:
                self._update_property({'Ngrids': self.energy_array['Ngrids'][-2]})
                break

    def _set_defaults_parameters(self):
        """
        Set the initial rel_cutoff and Ngrids parameters before starting optimization

        Returns
        -------
        updates the calculator attribute
        """
        self._update_property({'Rel_cutoff': self.loop_range['Rel_cutoff'][0],
                               'Ngrids': self.loop_range['Ngrids'][0]})

    def _return_results(self):
        """
        Return the results of the analysis.

        Returns
        -------
        Prints output to screen and saves images.
        """
        for item in self.energy_array:
            print(f"Energy cutoff: {self.energy_array[item]}")
            print(f"Force cutoff: {self.force_array[item]}")

            plt.plot(self.loop_range[:int(len(self.energy_array[item])-1)], self.energy_array[item],
                     'o-', label="Energy convergence")
            plt.plot(self.loop_range[:int(len(self.force_array[item])-1)], self.force_array[item],
                     'o-', label="Force convergence")
            plt.xlabel("Cutoff value (Ry)")
            plt.ylabel("Energy and Force Values")
            plt.grid()
            plt.legend()
            plt.savefig(f"Convergence_Results_{item}.pdf", dpi=800)

    def run_optimizer(self):
        """
        Run the optimizer
        Returns
        -------
        Plots images of the convergence and prints the ideal values.

        Notes
        -----
        TODO - Parallelize loops for all analysis
        """
        self._get_static_range()             # get the correct loop range
        self._set_defaults_parameters()      # set default parameters
        self._build_cell()                   # build the atomic cell if it is given
        self._run_cutoff_optimization()      # optimize the cutoff
        self._run_rel_cutoff_optimization()  # optimize the rel cutoff
        self._run_n_grids_optimization()     # optimize the n grids attribute
