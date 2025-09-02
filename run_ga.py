# 1- crossover
# 2- mutation


#Individual:  species * (ecoregions*4 + 7 + 8 ) + Seeding 
import sys
sys.path.append('.')
from generate_experiment_files import *
from collect_run_data import *
import numpy as np



def simulated_annealing(RNG, FL5_SPECIES, SPECIES_MAX_AGE,  PLOTS, gt_df,plots_df, plot_measurements_df, plot_min_measurements_df, MAX_ITER = 1000, T = 100000, alpha = 0.99992 ):
    best = -np.inf
    r2 = -np.inf
    exp_seed = RNG.integers(10000)
    exp_RNG = np.random.default_rng(seed=exp_seed)
    rs, sps, spps = None, None, None
    for i in range(MAX_ITER):
        prefix = f'exp_outputs/exp_{i:03d}'
        print(subprocess.run(['rm', '-rf',prefix]))
        print(subprocess.run(['cp', '-r', './template',prefix]))
        if rs is None:

            prefix, new_rs, new_sps, new_spps = generate_experiment(exp_RNG, FL5_SPECIES, SPECIES_MAX_AGE,  PLOTS, PREFIX = prefix)
        else:
            new_rs, new_sps, new_spps = mutate_biomass_succession_individual(RNG, rs,sps, spps, SPECIES_MAX_AGE)
            prefix= generate_experiment_for_individual(exp_RNG, FL5_SPECIES, SPECIES_MAX_AGE,  PLOTS, PREFIX = prefix, rs= new_rs, sps=new_sps, spps = new_spps)
        print(subprocess.run(['bash', '-c', f'cd {prefix} && ./run_simulation.sh']))
        df = read_csvs(prefix)
        outputs= compare_data(gt_df,plots_df, plot_measurements_df, plot_min_measurements_df, df)
        r2 =outputs[0]['r2_restricted']
        r2_unrestricted = outputs[0]['r2_unrestricted']
        best = max(best, outputs[0]['r2_restricted'])
        print("#"*80)
        print("#"*80)
        print(f"{prefix}: seed:{exp_seed} r2_linear_fit: {r2_unrestricted} r2:{r2} BEST:{best}")
        print("#"*80)
        print("#"*80)
        diff = r2 - best
        e = np.exp(diff / T)
        E = RNG.random()
        b = diff > 0 or e >= E
        print(f"Diff: {diff:e}, T:{T:e}, e:{e:e}, E:{E:e}, b:{b}")
        if b:
            rs, sps, spps = new_rs, new_sps, new_spps
        T *= alpha
        if T  == 0:
            return


def random_search(RNG, FL5_SPECIES, SPECIES_MAX_AGE,  PLOTS, gt_df,plots_df, plot_measurements_df, plot_min_measurements_df ):
    
    best = -np.inf
    for i in range(30):
        exp_seed = RNG.integers(10000)
        exp_RNG = np.random.default_rng(seed=exp_seed)
        prefix = f'exp_outputs/exp_{i:03d}'
        print(subprocess.run(['rm', '-rf',prefix]))
        print(subprocess.run(['cp', '-r', './template',prefix]))
        prefix = generate_experiment(exp_RNG, FL5_SPECIES, SPECIES_MAX_AGE,  PLOTS, PREFIX = prefix)
        print(subprocess.run(['bash', '-c', f'cd {prefix} && ./run_simulation.sh']))
        df = read_csvs(prefix)
        outputs= compare_data(gt_df,plots_df, plot_measurements_df, plot_min_measurements_df, df)
        r2 =outputs[0]['r2_restricted']
        r2_unrestricted = outputs[0]['r2_unrestricted']
        best = max(best, outputs[0]['r2_restricted'])
        print("#"*80)
        print("#"*80)
        print(f"{prefix}: seed:{exp_seed} r2_linear_fit: {r2_unrestricted} r2:{r2} BEST:{best}")
        print("#"*80)
        print("#"*80)

if __name__ == '__main__':

    import sys

    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 1337
    RNG = np.random.default_rng(seed=seed)
    print(f"GLOBAL SEED = {seed}")
    PLOTS = load_plots_data('./data_fl5_plot_genus_sp_ba_age_agb_20.csv')
    SPECIES_MAX_AGE = load_max_age('./max_treeage.csv')
    ##rs = get_core_species_params()
    ##rs = replace_species(rs,get_fl5_species())
    FL5_SPECIES = get_fl5_species()

    gt_df,plots_df, plot_measurements_df, plot_min_measurements_df = load_gt('./data_fl5_plot_genus_sp_ba_age_agb_20.csv')
    #random_search(RNG, FL5_SPECIES, SPECIES_MAX_AGE,  PLOTS, gt_df,plots_df, plot_measurements_df, plot_min_measurements_df )
    simulated_annealing(RNG, FL5_SPECIES, SPECIES_MAX_AGE,  PLOTS, gt_df,plots_df, plot_measurements_df, plot_min_measurements_df )

