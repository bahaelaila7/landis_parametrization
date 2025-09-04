import sys
import rasterio as rx
from collections import Counter, defaultdict
import pandas as pd
import pandasql as pdsql
import statsmodels.api as sm
from scipy.stats import pearsonr
from generate_experiment_files import *
sys.path.append('.')

def match_map_codes(prefix):
    #match initial communities with output communities
    with rx.open(f'{prefix}/initial-communities.tif') as src:
        in_data = src.read(1)
    with rx.open(f'{prefix}/output-community-0.tif') as src:
        out_data = src.read(1)
    assert in_data.shape == out_data.shape
    map_codes_dict = defaultdict(lambda: defaultdict(int))
    for i,o in zip(in_data.flatten(),out_data.flatten()):
        map_codes_dict[i][o] += 1
    rows = [(i, o, c) for i, os in map_codes_dict.items() for o,c in os.items()]
    map_df = pd.DataFrame(rows, columns = ["MapCodeIn", "MapCodeOut", "Multiplier"])

    return map_df




def read_csvs(prefix, timestep=5, timehorizon = 50):
    plot_mapcode_df = pd.read_csv(f'{prefix}/plot_mapcode_mapping.csv')
    map_df = match_map_codes(prefix)
    map_df.to_csv(f'{prefix}/output_map_codes_mapping.csv')
    # in_code ,out_code, multiplier


    total_df = None
    for yr in range(0,timehorizon+1, timestep):
        filename = f'{prefix}/community-input-file-{yr}.csv'
        print(f"reading {filename}")
        df = pd.read_csv(filename)
        df['SIM_YEAR'] = yr
        if total_df is None:
            total_df = df
        else:
            total_df = pd.concat([total_df, df])
    # reclassify tree by ageclass, same as initial communities
    # MapCode SpeciesName  CohortAge  CohortBiomass  CohortANPP  SIM_YEAR
    new_total_df = pdsql.sqldf(""" SELECT SIM_YEAR, MapCode MapCodeOut, SpeciesName, ageclass, sum(CohortBiomass) CohortBiomass
                               FROM (SELECT *, (CASE WHEN CohortAge >= 96 THEN 100 ELSE MAX(1,(((CohortAge-1) / 5)+1)) * 5 END) ageclass 
                                     FROM total_df
                                     )
                               GROUP BY SIM_YEAR, MapCode, SpeciesName, ageclass""")

    new_total_df.to_csv(f"{prefix}/sim.csv")
    new_total_mapped_df = pdsql.sqldf(""" SELECT df.SIM_YEAR, map_df.MapCodeIn MapCode, df.SpeciesName, df.ageclass, sum(CohortBiomass * map_df.Multiplier) CohortBiomass
                                        FROM new_total_df df
                                        LEFT OUTER JOIN map_df ON map_df.MapCodeOut = df.MapCodeOut
                                        GROUP BY df.SIM_YEAR, map_df.MapCodeIn, df.SpeciesName, df.ageclass""")
    new_total_mapped_df.to_csv(f"{prefix}/sim_mapped.csv")

    new_total_mapped_plot_df = pdsql.sqldf( """ SELECT p.STATECD, p.UNITCD, p.COUNTYCD, p.PLOT, d.*
                                FROM new_total_mapped_df d
                                LEFT OUTER JOIN plot_mapcode_df p ON d.MapCode = p.MapCode""")
    new_total_mapped_plot_df.to_csv(f"{prefix}/sim_mapped_plot.csv")

    return new_total_mapped_plot_df

def compare_data(prefix, gt_df, plot_measurements_df, plot_min_measurements_df, total_df):
    gt = gt_df
    #print(total_df)
    #print(len(total_df))
    total_p_df = pdsql.sqldf( """ SELECT  p.*,  min_df.min_measyear, (min_df.min_measyear + p.SIM_YEAR) real_year
                                FROM total_df p 
                                LEFT OUTER JOIN plot_min_measurements_df min_df
                                    ON p.STATECD = min_df.STATECD
                                    AND p.UNITCD = min_df.UNITCD
                                    AND p.COUNTYCD = min_df.COUNTYCD
                                    AND p.PLOT = min_df.PLOT""")
    total_p_df.to_csv(f"{prefix}/sim_min_year.csv")

    total_pp_df = pdsql.sqldf("""
                        SELECT sim.*, meas_df.STATECD mSTATECD, meas_df.UNITCD mUNITCD, meas_df.COUNTYCD mCOUNTYCD, meas_df.PLOT mPLOT, meas_df.MEASYEAR mMEASYEAR
                        FROM total_p_df sim
                        --- Only take years that have FIA measurement for the plot
                        JOIN plot_measurements_df meas_df ON
                        sim.STATECD = meas_df.STATECD AND
                        sim.UNITCD = meas_df.UNITCD AND
                        sim.COUNTYCD = meas_df.COUNTYCD AND
                        sim.PLOT = meas_df.PLOT AND
                        meas_df.MEASYEAR = sim.real_year

    """)
    total_pp_df.to_csv(f"{prefix}/sim_min_year_filtered.csv")
    gt_df.to_csv(f"{prefix}/gt_df.csv")

    h = pdsql.sqldf("""SELECT 
                                sim.MapCode,sim.STATECD, sim.UNITCD, sim.COUNTYCD, sim.PLOT, sim.real_year, sim.SpeciesName, sim.ageclass, sim.SIM_YEAR, sim.min_measyear, Coalesce(sim.CohortBiomass,0) CohortBiomass,
                                d.STATECD dSTATECD, d.UNITCD dUNITCD, d.COUNTYCD dCOUNTYCD, d.PLOT dPLOT, d.SIM_YEAR dSIM_YEAR, d.measyear, d.species_symbol, d.ageclass dageclass, COALESCE(d.drybio_ag,0) drybio_ag
                        FROM total_pp_df sim
                        FULL OUTER JOIN gt_df d ON
                        sim.STATECD = d.STATECD AND
                        sim.UNITCD = d.UNITCD AND
                        sim.COUNTYCD = d.COUNTYCD AND
                        sim.PLOT = d.PLOT AND
                        sim.real_year = d.MEASYEAR AND
                        sim.SpeciesName = d.SPECIES_SYMBOL AND 
                        sim.ageclass = d.ageclass 
                        --WHERE COALESCE(sim.SIM_YEAR, d.SIM_YEAR) > 0
                        ORDER BY sim.SIM_YEAR, dSIM_YEAR,
                                sim.STATECD, dSTATECD,
                                sim.UNITCD, dUNITCD,
                                sim.COUNTYCD, dCOUNTYCD,
                                sim.PLOT, dPLOT,
                                sim.SpeciesName, species_symbol,
                                sim.ageclass, dageclass

                        """)
    '''
    h = pdsql.sqldf("""SELECT 
                                sim.MapCode,
                                Coalesce(SIM.STATECD,d.STATECD) STATECD,
                                Coalesce(SIM.UNITCD,d.UNITCD) UNITCD,
                                Coalesce(SIM.COUNTYCD,d.COUNTYCD) COUNTYCD,
                                Coalesce(SIM.PLOT,d.PLOT) PLOT,
                                Coalesce(SIM.SpeciesName,d.SPECIES_SYMBOL) species_symbol,
                                Coalesce(SIM.ageclass,d.ageclass) ageclass,
                                Coalesce(SIM.SIM_YEAR, d.SIM_YEAR) SIM_YEAR,
                                sim.min_measyear,
                                sim.real_year,
                                d.MEASYEAR,
                                Coalesce(SIM.CohortBiomass, 0) CohortBiomass,
                                Coalesce(d.drybio_ag,0) drybio_ag
                        FROM total_pp_df sim
                        FULL OUTER JOIN gt_df d ON
                        sim.STATECD = d.STATECD AND
                        sim.UNITCD = d.UNITCD AND
                        sim.COUNTYCD = d.COUNTYCD AND
                        sim.PLOT = d.PLOT AND
                        sim.real_year = d.MEASYEAR AND
                        sim.SpeciesName = d.SPECIES_SYMBOL AND 
                        sim.ageclass = d.ageclass
                        --sim.SIM_YEAR  = d.SIM_YEAR 

                        """)
    '''
    #print(dff)
    '''
    h = pdsql.sqldf("""SELECT 
                                Coalesce(SIM.STATECD,d.STATECD) STATECD,
                                Coalesce(SIM.UNITCD,d.UNITCD) UNITCD,
                                Coalesce(SIM.COUNTYCD,d.COUNTYCD) COUNTYCD,
                                Coalesce(SIM.PLOT,d.PLOT) PLOT,
                                Coalesce(SIM.SpeciesName,d.SPECIES_SYMBOL) species_symbol,
                                Coalesce(SIM.ageclass,d.ageclass) ageclass,
                                Coalesce(SIM.SIM_YEAR, d.SIM_YEAR) SIM_YEAR,
                                min_df.min_measyear,
                                d.MEASYEAR,
                                Coalesce(SIM.CohortBiomass, 0) CohortBiomass,
                                Coalesce(d.drybio_ag,0) drybio_ag
                        FROM total_p_df sim
                        --- Only take years that have FIA measurement for the plot
                        JOIN plot_min_measurements_df min_df ON
                        sim.STATECD = min_df.STATECD AND
                        sim.UNITCD = min_df.UNITCD AND
                        sim.COUNTYCD = min_df.COUNTYCD AND
                        sim.PLOT = min_df.PLOT
                        JOIN plot_measurements_df meas_df ON
                        sim.STATECD = meas_df.STATECD AND
                        sim.UNITCD = meas_df.UNITCD AND
                        sim.COUNTYCD = meas_df.COUNTYCD AND
                        sim.PLOT = meas_df.PLOT AND
                        meas_df.MEASYEAR = sim.SIM_YEAR + min_df.min_measyear 

                        FULL OUTER JOIN gt_df d ON
                        sim.STATECD = d.STATECD AND
                        sim.UNITCD = d.UNITCD AND
                        sim.COUNTYCD = d.COUNTYCD AND
                        sim.PLOT = d.PLOT AND
                        sim.SpeciesName = d.SPECIES_SYMBOL AND 
                        sim.ageclass = d.ageclass AND
                        sim.SIM_YEAR  = d.SIM_YEAR 
                        ORDER BY SIM_YEAR,
                                 MapCode,
                                 STATECD,
                                 UNITCD,
                                 COUNTYCD,
                                 PLOT,
                                 SpeciesName
                        """)
    '''
    h.to_csv(f'{prefix}/cohort_matching.csv')


    ##the hypothesis is that CohortBiomass is a good estimator of drybio_ag
    ##That means R^2 is going to be between CohortBiomass and drybio_ag directly
    ## NOT drybio_ag~CohortBiomass. Which is still easy
    ## but a p-value is a bit tricky since we do not have a statistical model.
    ## but we can invent one! drybio_ag = 0 + 1 * drybio_ag + err, Where err ~ N(0,1)
    ## But now the null hypothesis will be whether intercept = 0 and slope =1
    ## or there is a better calibration.
    ## the best calibration would be by OLS, and then F-test between the errors of both models
    ## 


    def r2pval(g):
        from scipy.special import fdtr

        N = len(g['drybio_ag'])
        # Direct (restricted): drybio_ag = CohortBiomass + err
        se_restricted = (g['drybio_ag'] - g['CohortBiomass']).pow(2)
        sse_restricted = se_restricted.sum()
        variance = (g['drybio_ag'] - g['drybio_ag'].mean()).pow(2).sum()
        if N <= 2:
            return pd.Series({'F':np.inf,'pval':1.0, 'sse_restricted':sse_restricted,'var':variance, 'r2_restricted': 1.0 - sse_restricted/(variance + 1e-14), 'sse_unrestricted':np.inf, 'r2_unrestricted':-np.inf, 'pval_unrestricted':1.0})

        #print("#"*120)
        # OLS (unrestricted): drybio_ag = intercept + slope* CohortBiomass + err
        X = sm.add_constant(g['CohortBiomass'])
        #print('--', len(X))
        model = sm.OLS(g['drybio_ag'], X).fit()
        #print('^^', len(X))
        se_unrestricted = model.resid.pow(2)
        sse_unrestricted = se_unrestricted.sum()


        df_restricted = len(X)
        df_unrestricted = len(X) - 2 # slope and intercept
        df_diff =  df_restricted - df_unrestricted  # 2


        #print('---', sse_unrestricted, df_unrestricted, N, len(g), len(X))
        sse_per_degree_of_freedom = sse_unrestricted/df_unrestricted
        #print('^^^', sse_unrestricted, df_unrestricted, N, len(g), len(X))

        #print(sse_restricted, sse_unrestricted, sse_restricted - sse_unrestricted)
        deviation_in_sse_per_degree_of_freedom_diff = (sse_restricted - sse_unrestricted) / df_diff
        if sse_per_degree_of_freedom == 0:
            F = np.inf
        else:
            F = deviation_in_sse_per_degree_of_freedom_diff / sse_per_degree_of_freedom
        pval = 1- fdtr(df_diff, df_unrestricted, F)
        #print("@"*120)
        r2_unrestricted = model.rsquared if model.centered_tss != 0.0 else -np.inf
        #print("@"*110, model.centered_tss)
        #print(deviation_in_sse_per_degree_of_freedom_diff)
        #print(sse_per_degree_of_freedom)
        #print(F)

        return pd.Series({'F':F,'pval':pval, 'sse_restricted':sse_restricted,'var':variance, 'r2_restricted': 1.0 - sse_restricted/(variance + 1e-14), 'sse_unrestricted':sse_unrestricted, 'r2_unrestricted':r2_unrestricted, 'pval_unrestricted':model.pvalues['CohortBiomass']})


    def stratified_apply(g,func,stratify_by_cols = []):
        if stratify_by_cols:
            return g.groupby(stratify_by_cols).apply(func,include_groups=False )
        else:
            return func(g)


    overall_output = stratified_apply(h,r2pval)
    #print(output)
    by_species_output = stratified_apply(h,r2pval, ['species_symbol'])
    by_year_output = stratified_apply(h,r2pval, ['SIM_YEAR'])
    by_species_year_output = stratified_apply(h,r2pval, ['species_symbol','SIM_YEAR'])
    #print(output)
    return (overall_output, by_species_output, by_year_output, by_species_year_output)

def load_gt(plots_csv):
    gt = pd.read_csv(plots_csv)
    gt = gt.sort_values(by=['invyr','measyear','measmon','measday','statecd','unitcd','countycd','plot','species_symbol','ageclass','drybio_ag'])
    print(gt)
    print(gt.shape)
    plot_measurements_df = pdsql.sqldf("""SELECT distinct STATECD,UNITCD,COUNTYCD, PLOT, MEASYEAR FROM gt""")
    plot_min_measurements_df = pdsql.sqldf("""SELECT STATECD,UNITCD,COUNTYCD, PLOT, min(MEASYEAR) min_measyear FROM gt
                                            GROUP BY STATECD,UNITCD,COUNTYCD, PLOT""")
    sql = """SELECT gt.*,(gt.MEASYEAR - mgt.min_measyear) SIM_YEAR
            FROM gt
            JOIN plot_min_measurements_df mgt
                ON gt.STATECD = mgt.STATECD
                AND gt.UNITCD = mgt.UNITCD
                AND gt.COUNTYCD = mgt.COUNTYCD
                AND gt.PLOT = mgt.PLOT
                """
    gt_df = pdsql.sqldf(sql)
    return (gt_df, plot_measurements_df, plot_min_measurements_df)

def run_comparison(plots_csv, prefix):
    df = read_csvs(prefix)
    gt_df, plot_measurements_df, plot_min_measurements_df = load_gt(plots_csv)
    return compare_data(prefix, gt_df, plot_measurements_df, plot_min_measurements_df, df)

if __name__ == '__main__':
    import sys
    prefix = sys.argv[1]
    from pprint import pprint
    pprint(run_comparison('./data_fl5_plot_genus_sp_ba_age_agb_20.csv', prefix))




    #year, plot, ecoregion, species, cohortAge, cohortBiomass

