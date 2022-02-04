"""
Example parameter file for the running of a multi-frequency, multi-epoch
radiative transfer calculation and subsequent synthetic imaging with a variety
of telescopes/telescope-configurations.

Use would be with a RaJePy.classes.ModelRun class instance e.g.:

pipeline = RaJePy.classes.ModelRun('/full/path/to/example-pipeline-params.py')
"""

import os
import numpy as np

params = {'min_el':    20.,    # Minimum elevation for synthetic observations, deg
          'dcys':      {"model_dcy": os.sep.join([os.getcwd(), 'output_run_4'])},  # Output root directory
          # Continuum observations
          'continuum': {'times':  np.array([0.0]),  # Model times, yr
                        'freqs':  np.array([5., 15., 5.]) * 1e9,  # Frequencies of observations, Hz
                        't_obs':  np.array([358, 4800, 59400]),  # Total on-source times, s
                        'tscps':  np.array([('VLA', 'A'), ('VLA', 'A'), ('EMERLIN', '0')]),  # List of 2-tuples of (telescope, configuration)
                        't_ints': np.array([2, 2, 8]),  # Visibility integration times, s
                        'bws':    np.array([4e9, 6e9, .5e9]),  # Observational bandwidth, Hz
                        'chanws': np.array([1e8, 1e8, 1e8])},  # Channel widths, Hz
          # Radio recombination line observations
          'rrls':      {'times':  np.array([]),  # Model times, yr
                        'lines':  np.array([]),  # RRL lines to observe (Element+n+dn)
                        't_obs':  np.array([]),  # Total on-source times, s
                        'tscps':  np.array([]),  # List of 2-tuples of (telescope, configuration)
                        't_ints': np.array([]),  # Visibility integration times, s
                        'bws':    np.array([]),  # Observational bandwidth, Hz
                        'chanws': np.array([])},  # Channel widths, Hz
          }
