import sys
import pandas as pd
import pandasql as pdsql
import statsmodels.api as sm
from scipy.stats import pearsonr
from generate_experiment_files import *
sys.path.append('.')

def read_csvs(prefix, timestep=5, timehorizon = 50):
    total_df = None
    for yr in range(0,timehorizon, timestep):
        df = pd.read_csv(f'{prefix}/community-input-file-{yr}.csv')
        df['SIM_YEAR'] = yr
        if total_df is None:
            total_df = df
        else:
            total_df = pd.concat([total_df, df])
    # reclassify tree by ageclass, same as initial communities
    # MapCode SpeciesName  CohortAge  CohortBiomass  CohortANPP  SIM_YEAR
    new_total_df = pdsql.sqldf(""" SELECT SIM_YEAR, MapCode, SpeciesName, ageclass, sum(CohortBiomass) CohortBiomass
                               FROM (SELECT *, (CASE WHEN CohortAge >= 96 THEN 100 ELSE MAX(1,(((CohortAge-1) / 5)+1)) * 5 END) ageclass 
                                     FROM total_df
                                     )
                               GROUP BY SIM_YEAR, MapCode, SpeciesName, ageclass""")

    return new_total_df

def compare_data(prefix, gt_df,plots_df, plot_measurements_df, plot_min_measurements_df, total_df):
    gt = gt_df
    #print(total_df)
    #print(len(total_df))
    total_p_df = pdsql.sqldf( """ SELECT p.STATECD, p.UNITCD, p.COUNTYCD, p.PLOT,d.* FROM total_df d
                                JOIN plots_df p ON d.MapCode = p.MapCode""")

    #print(dff)
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

                        """)
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
    s = set()
    plots = []

    i = 0
    for _,row in gt.iterrows():
        plot_id = (row['statecd'],row['unitcd'],row['countycd'], row['plot'])
        if plot_id not in s:
            plots.append((i, *plot_id))
            s.add(plot_id)
            i+=1
    plots_df = pd.DataFrame(plots , columns=['MapCode','STATECD','UNITCD','COUNTYCD','PLOT'])
    #print(df)
    plot_measurements_df = pdsql.sqldf("""SELECT distinct STATECD,UNITCD,COUNTYCD, PLOT, MEASYEAR FROM gt""")
    plot_min_measurements_df = pdsql.sqldf("""SELECT STATECD,UNITCD,COUNTYCD, PLOT, min(MEASYEAR) min_measyear FROM gt
                                            GROUP BY STATECD,UNITCD,COUNTYCD, PLOT""")
    sql = """SELECT gt.*,(gt.MEASYEAR - mgt.min_measyear) SIM_YEAR FROM gt JOIN (SELECT STATECD, UNITCD, COUNTYCD, PLOT, min(MEASYEAR) min_measyear FROM gt) mgt
                                          on gt.STATECD = mgt.STATECD AND gt.UNITCD = mgt.UNITCD AND gt.COUNTYCD = mgt.COUNTYCD
                                            AND gt.PLOT = mgt.PLOT
                """
    gt_df = pdsql.sqldf(sql)
    return (gt_df,plots_df, plot_measurements_df, plot_min_measurements_df)

def run_comparison(plots_csv, prefix):
    df = read_csvs(prefix)
    gt_df,plots_df, plot_measurements_df, plot_min_measurements_df = load_gt(plots_csv)
    return compare_data(prefix, gt_df,plots_df, plot_measurements_df, plot_min_measurements_df, df)

if __name__ == '__main__':
    import sys
    prefix = sys.argv[1]
    from pprint import pprint
    pprint(run_comparison('./data_fl5_plot_genus_sp_ba_age_agb_20.csv', prefix))




    #year, plot, ecoregion, species, cohortAge, cohortBiomass

