import os
import string
import numpy as np
import pandas as pd
import Simulation as sml
from MergingPolicies import random_merge, largest_first_merge, weighted_distribution_by_size
MAIN_DATA_DIRECTORY_PATH = '/Users/yishaiazabary/Desktop/My Career/OriAvinoamSimulation/Data'


class SimulationPerExperiment:
    def __init__(self, exp_prefix: str):
        self.exp_prefix = exp_prefix
        self.exp_path = os.sep.join([MAIN_DATA_DIRECTORY_PATH,exp_prefix[:-10], '{0}.csv'.format(exp_prefix)])
        self.exp_df = pd.read_csv(self.exp_path)
        self.max_number_of_cells = np.max(self.exp_df.loc[:, 'CTRn=1'].values)
        self.number_of_frames = len(self.exp_df.loc[:, 'CTRn=1'].values)

    def calc_experiment_fusions_in_frame(self):
        n_fusions_per_frame = [0]
        cols = ['CTRn=2', 'CTRn=3', 'CTRn=4']
        initial_ctrs = {
            cols[0]: self.exp_df.loc[0, 'CTRn=2'],
            cols[1]: self.exp_df.loc[0, 'CTRn=3'],
            cols[2]: self.exp_df.loc[0, 'CTRn=4']
        }
        for frame_idx, next_frame in enumerate(self.exp_df.loc[1:, ['CTRn=2', 'CTRn=3','CTRn=4']].values):
            fusions_in_frame = 0
            for colIDX in range(len(cols)):
                if initial_ctrs[cols[colIDX]] != next_frame[colIDX]:
                    fusions_in_frame += abs(int(initial_ctrs[cols[colIDX]] - next_frame[colIDX]))
                    initial_ctrs[cols[colIDX]] = next_frame[colIDX]
            n_fusions_per_frame.append(fusions_in_frame)
        return n_fusions_per_frame


n_replicates = 1000
for exp_prefix in ['200203_S17_processed', '200203_S19_processed', '200203_S22_processed', '200203_S24_processed', '200604_S06_processed','200604_S09_processed']:

# exp_prefix = '200203_S17_processed'
    spe = SimulationPerExperiment(exp_prefix)
    n_fusions_per_frame = spe.calc_experiment_fusions_in_frame()
    max_number_of_cells = spe.max_number_of_cells
    number_of_frames = spe. number_of_frames
    for i in range(n_replicates):
        print('{0}: replica No. {1}'.format(exp_prefix, i))
        simulation = sml.Simulation(initial_cell_number=max_number_of_cells,
                                    merging_policy=weighted_distribution_by_size,
                                    simulation_name='{}_simulation'.format(exp_prefix))
        simulation.run_simulation_fusions_by_fusions_list(n_fusions_per_frame=n_fusions_per_frame)
        # simulation.convert_simulation_to_data_csv(path=os.sep.join([MAIN_DATA_DIRECTORY_PATH, exp_prefix[:-10], '1000_randoms', '{0}_{1}_{2}.csv'.format(exp_prefix, random_merge.__name__, i)]))
        simulation.convert_simulation_to_data_csv(path=os.sep.join([MAIN_DATA_DIRECTORY_PATH, exp_prefix[:-10], '1000_weighted', '{0}_{1}_{2}.csv'.format(exp_prefix, weighted_distribution_by_size.__name__, i)]))
