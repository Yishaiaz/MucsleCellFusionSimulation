import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
MAIN_DATA_DIRECTORY_PATH = '/Users/yishaiazabary/Desktop/My Career/OriAvinoamSimulation/Data'


class BootstrappingPerExperiment:
    def __init__(self, exp_prefix, type_of_simulation, col_to_calc_by, number_of_replications, frame_num):
        self.exp_prefix = exp_prefix
        self.type_of_simulation = type_of_simulation
        self.col_to_calc_by = col_to_calc_by
        type_of_simulation = '{0}_{1}'.format(int(number_of_replications), type_of_simulation)
        self.path_to_simulations_directory = os.sep.join([MAIN_DATA_DIRECTORY_PATH, exp_prefix, type_of_simulation])
        self.number_of_replications = number_of_replications
        self.frame_num = frame_num
        self.path_to_true_exp_file = os.sep.join([MAIN_DATA_DIRECTORY_PATH, exp_prefix, ''.join([exp_prefix, '_processed.csv'])])
        self.all_simulations_data_points = np.ndarray((self.number_of_replications, self.frame_num))
        self.simulations_dist_averages = np.ndarray(self.frame_num)
        self.p_values_per_frame = np.zeros(frame_num)
        entire_exp_df = pd.read_csv(self.path_to_true_exp_file)
        sum_of_totals = entire_exp_df.loc[:, 'CTRn=1'].values + \
                        entire_exp_df.loc[:, 'CTRn=2'].values + \
                        entire_exp_df.loc[:, 'CTRn=3'].values + \
                        entire_exp_df.loc[:, 'CTRn=4'].values
        self.exp_data_points_normalized = entire_exp_df.loc[:, self.col_to_calc_by].values / sum_of_totals

    def read_simulations_col_data_points(self):
        files = (file for file in os.listdir(self.path_to_simulations_directory)
                 if os.path.isfile(os.path.join(self.path_to_simulations_directory, file)))
        idx = 0
        for filename in files:
            if filename.endswith('.csv'):
                path_to_single_file = os.sep.join([self.path_to_simulations_directory, filename])
                entire_df =  pd.read_csv(path_to_single_file)
                single_file_col_data =entire_df.loc[:, self.col_to_calc_by].values
                sum_of_totals = entire_df.loc[:, 'CTRn=1'].values + \
                                entire_df.loc[:, 'CTRn=2'].values + \
                                entire_df.loc[:, 'CTRn=3'].values + \
                                entire_df.loc[:, 'CTRn=4'].values
                single_file_data_normalized = single_file_col_data / sum_of_totals
                self.all_simulations_data_points[idx] = single_file_data_normalized
                idx += 1

    def calc_p_values_per_frame(self):
        for frame_idx, exp_frame_dist in enumerate(self.exp_data_points_normalized):
            frame_data_from_simulations = self.all_simulations_data_points[:, frame_idx]
            self.simulations_dist_averages[frame_idx] = np.average(frame_data_from_simulations)
            num_of_equal_or_above = len(frame_data_from_simulations[frame_data_from_simulations >= exp_frame_dist])
            self.p_values_per_frame[frame_idx] = 0.001 if num_of_equal_or_above / self.number_of_replications <= 0.001 else num_of_equal_or_above / self.number_of_replications

    def plot_col_data_with_p_values(self):
        fig, ax = plt.subplots()
        lines_and_legends = {}

        # plot simulation dist averages
        coding = 'gx'
        sim_data = self.simulations_dist_averages
        ax.plot(sim_data, coding)
        lines_and_legends['{} Simulation Distribution Averages'.format(self.type_of_simulation.capitalize())] = coding
        # plot true exp dist
        coding = 'r+'
        exp_data = self.exp_data_points_normalized
        ax.plot(exp_data, coding)
        lines_and_legends['True Experiment Distribution'] = coding
        # plot p values
        coding = 'b*'
        pval = self.p_values_per_frame
        ax.plot(pval, coding)
        lines_and_legends['P-Values<='] = coding

        ax.set_title('EXP:{0}\nBootstrap P-Values of Simulation type:{1}\nReplications:{2}'.format(self.exp_prefix, self.type_of_simulation, self.number_of_replications))
        ax.set_ylabel('fraction')
        ax.set_xlabel('time (h)')

        artist_patches = []
        for leg_label, coding in lines_and_legends.items():
            artist_patches.append(mlines.Line2D([], [],
                                                color=coding[:1],
                                                marker=coding[1:],
                                                label=leg_label))
        # annotate p_values
        for idx, p_value in enumerate(self.p_values_per_frame):
            text = '{:.3f}≥'.format(p_value) if p_value != 1 else str(1)
            ax.annotate(text, xy=(idx, p_value), xytext=(idx+0.05, p_value), fontsize=6, rotation=40)

        ax.legend(handles=artist_patches, loc=6)
        plt.show()
        path_to_plot = os.sep.join([MAIN_DATA_DIRECTORY_PATH, self.exp_prefix, 'plots', '{0}_{1}_{2}.png'.format(exp_prefix, self.type_of_simulation, 'pvalues')])
        fig.savefig(path_to_plot, dpi=300)


for exp_prefix in ['200203_S17', '200203_S19', '200203_S22', '200203_S24', '200604_S06', '200604_S09']:
    # exp_prefix = '200604_S09'
    type_of_simulation = 'weighted'
    col_to_calc_by = 'CTRn=4'
    number_of_replications = 1000
    frame_num = 16

    bspe = BootstrappingPerExperiment(exp_prefix, type_of_simulation, col_to_calc_by, number_of_replications, frame_num)
    bspe.read_simulations_col_data_points()
    bspe.calc_p_values_per_frame()
    bspe.plot_col_data_with_p_values()
