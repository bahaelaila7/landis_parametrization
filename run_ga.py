# 1- crossover
# 2- mutation


#Individual:  species * (ecoregions*4 + 7 + 8 ) + Seeding 
import sys
sys.path.append('.')
from generate_experiment_files import *
from collect_run_data import *
import numpy as np



def simulated_annealing(RNG, FL5_SPECIES, SPECIES_MAX_AGE,  PLOTS, gt_df,plots_df, plot_measurements_df, plot_min_measurements_df, MAX_ITER = 1000, INITIAL_T = 10000, alpha = 0.9992 ):
    best = -np.inf
    r2 = -np.inf
    exp_seed = RNG.integers(10000)
    exp_RNG = np.random.default_rng(seed=exp_seed)
    rs, sps, spps = None, None, None
    mutate_type = 'random'
    gauss_rate = 0.1
    T = INITIAL_T
    for i in range(MAX_ITER):
        prefix = f'exp_outputs/exp_{i:03d}'
        print(subprocess.run(['rm', '-rf',prefix]))
        print(subprocess.run(['cp', '-r', './template',prefix]))
        if rs is None:

            prefix, new_rs, new_sps, new_spps = generate_experiment(exp_RNG, FL5_SPECIES, SPECIES_MAX_AGE,  PLOTS, PREFIX = prefix)
        else:
            new_rs, new_sps, new_spps = mutate_biomass_succession_individual(RNG, rs,sps, spps, SPECIES_MAX_AGE, mutate_type=mutate_type, gauss_rate=gauss_rate)
            prefix= generate_experiment_for_individual(exp_RNG, FL5_SPECIES, SPECIES_MAX_AGE,  PLOTS, PREFIX = prefix, rs= new_rs, sps=new_sps, spps = new_spps)
        print(subprocess.run(['bash', '-c', f'cd {prefix} && ./run_simulation.sh']))
        df = read_csvs(prefix)
        outputs= compare_data(gt_df,plots_df, plot_measurements_df, plot_min_measurements_df, df)
        new_r2 =outputs[0]['r2_restricted']
        new_r2_unrestricted = outputs[0]['r2_unrestricted']
        if new_r2 > best:
            best = new_r2
            with open('exp_outputs/best_params.txt', 'w') as f:
                f.write(f'best_r2={best}\n')
                f.write(f'r2_linear={new_r2_unrestricted}\n')
                f.write(f'pval_unrestricted={outputs[0]['pval_unrestricted']}\n')
                f.write(f'pval={outputs[0]['pval']}\n')
                f.write(f'seed={exp_seed}\n')
                f.write(f'mutate_type={mutate_type}\n')
                f.write(f'gauss_rate={gauss_rate}\n')
                f.write(f'INITIAL_T={INITIAL_T}\n')
                f.write(f'alpha={alpha}\n')
                f.write(f'prefix={prefix}')
                for r in new_rs:
                    f.write(str(r) + '\n')
                for sp in new_sps:
                    f.write(str(sp) + '\n')
                for spp in new_spps:
                    f.write(str(spp) + '\n')

        print("#"*80)
        print("#"*80)
        print(f"{prefix}: seed:{exp_seed} r2_linear_fit: {new_r2_unrestricted} r2:{new_r2} BEST:{best}")
        print("#"*80)
        print("#"*80)
        diff = new_r2 - r2
        e = np.exp(diff / T)
        E = RNG.random()
        b = diff > 0 or e >= E
        print(f"Diff: {diff:e}, T:{T:e}, e:{e:e}, E:{E:e}, b:{b}")
        if b:
            r2, rs, sps, spps = new_r2, new_rs, new_sps, new_spps
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

