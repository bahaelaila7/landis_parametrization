import sys
import pandas as pd
import pandasql as pdsql
sys.path.append('.')
from generate_experiment_files import *

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
def compare_data(plots,total_df):
    df = pd.read_csv('./data_fl5_plot_genus_sp_ba_age_agb_20.csv')
    df = df.sort_values(by=['invyr','measyear','measmon','measday','statecd','unitcd','countycd','plot','species_symbol','ageclass','drybio_ag'])
    s = set()
    plots = []

    i = 0
    for _,row in df.iterrows():
        plot_id = (row['statecd'],row['unitcd'],row['countycd'], row['plot'])
        if plot_id not in s:
            plots.append((i, *plot_id))
            s.add(plot_id)
            i+=1
    plots_df = pd.DataFrame(plots , columns=['MapCode','STATECD','UNITCD','COUNTYCD','PLOT'])
    #print(total_df)
    #print(len(total_df))
    total_p_df = pdsql.sqldf( """ SELECT p.STATECD, p.UNITCD, p.COUNTYCD, p.PLOT,d.* FROM total_df d
                                JOIN plots_df p ON d.MapCode = p.MapCode""")

    gt = pd.read_csv('./data_fl5_plot_genus_sp_ba_age_agb_20.csv')
    plot_measurements_df = pdsql.sqldf("""SELECT distinct STATECD,UNITCD,COUNTYCD, PLOT, MEASYEAR FROM gt""")
    plot_min_measurements_df = pdsql.sqldf("""SELECT STATECD,UNITCD,COUNTYCD, PLOT, min(MEASYEAR) min_measyear FROM gt
                                            GROUP BY STATECD,UNITCD,COUNTYCD, PLOT""")
    sql = """SELECT gt.*,(gt.MEASYEAR - mgt.min_measyear) SIM_YEAR FROM gt JOIN (SELECT STATECD, UNITCD, COUNTYCD, PLOT, min(MEASYEAR) min_measyear FROM gt) mgt
                                          on gt.STATECD = mgt.STATECD AND gt.UNITCD = mgt.UNITCD AND gt.COUNTYCD = mgt.COUNTYCD
                                            AND gt.PLOT = mgt.PLOT
                """
    dff= pdsql.sqldf(sql)
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

                        FULL OUTER JOIN dff d ON
                        sim.STATECD = d.STATECD AND
                        sim.UNITCD = d.UNITCD AND
                        sim.COUNTYCD = d.COUNTYCD AND
                        sim.PLOT = d.PLOT AND
                        sim.SpeciesName = d.SPECIES_SYMBOL AND 
                        sim.ageclass = d.ageclass AND
                        sim.SIM_YEAR  = d.SIM_YEAR 

                        """)
    import statsmodels.api as sm
    from scipy.stats import pearsonr
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
        print(g['drybio_ag'])
        print(len(g['drybio_ag']))

        # Direct (restricted): drybio_ag = CohortBiomass + err
        se_restricted = (g['drybio_ag'] - g['CohortBiomass']).pow(2)
        sse_restricted = se_restricted.sum()
        variance = (g['drybio_ag'] - g['drybio_ag'].mean()).pow(2).sum()

        # OLS (unrestricted): drybio_ag = intercept + slope* CohortBiomass + err
        X = sm.add_constant(g['CohortBiomass'])
        print(len(X))
        model = sm.OLS(g['drybio_ag'], X).fit()
        se_unrestricted = model.resid.pow(2)
        sse_unrestricted = se_unrestricted.sum()


        df_restricted = len(X)
        df_unrestricted = len(X) - 2 # slope and intercept
        df_diff =  df_restricted - df_unrestricted  # 2

        sse_per_degree_of_freedom = sse_unrestricted/df_unrestricted

        #print(sse_restricted, sse_unrestricted, sse_restricted - sse_unrestricted)
        deviation_in_sse_per_degree_of_freedom_diff = (sse_restricted - sse_unrestricted) / df_diff
        F = deviation_in_sse_per_degree_of_freedom_diff / sse_per_degree_of_freedom
        pval = 1- fdtr(df_diff, df_unrestricted, F)
        #print(deviation_in_sse_per_degree_of_freedom_diff)
        #print(sse_per_degree_of_freedom)
        #print(F)

        return pd.Series({'F':F,'pval':pval, 'sse_restricted':sse_restricted,'var':variance, 'r2_restricted': 1.0 - sse_restricted/variance, 'sse_unrestricted':sse_unrestricted, 'r2_unrestricted':model.rsquared, 'pval_unrestricted':model.pvalues['CohortBiomass']})


    def stratified_apply(g,func,stratify_by_cols = []):
        if stratify_by_cols:
            return g.groupby(stratify_by_cols ).apply(func)
        else:
            return func(g)


    output = stratified_apply(h,r2pval)
    print(output)
    print(len(output))
    return



    def corr_stats(g):
        r, p = pearsonr(g['drybio_ag'], g['CohortBiomass'])
        return pd.Series({"r2":r**2, "pval":p})
    print(h)
    print(h.groupby('species_symbol').apply(linreg))
    return


    mse_df = pdsql.sqldf("""SELECT species_symbol, SIM_YEAR, sum(pow(coalesce(drybio_ag,0) - coalesce(cohortbiomass,0),2)) agb_mse
                            FROM h
                          GROUP BY species_symbol, SIM_YEAR""")
    year_mse = pdsql.sqldf("""SELECT SIM_YEAR, sum(agb_mse) agb_mse FROM mse_df GROUP BY SIM_YEAR""")
    total_mse = pdsql.sqldf("""SELECT sum(agb_mse) FROM mse_df""")

    print(mse_df)
    print(year_mse)
    print(total_mse)
    #total_df['ageclass'] = total_df['CohortAge']
    #print(total_df.join(gt,lsuffix='_sim', how= 'outer'))
    #plots = load_plots_data('./data_fl5_plot_genus_sp_ba_age_agb_20.csv')
    #print(len(plots))
    #mse = 0
    #for num, plot in enumerate(plots):
    #    measyear0 = plot.measurements[0].measyear
    #    print(num, [ meas.measyear - measyear0 for meas in plot.measurements])


if __name__ == '__main__':
    import sys
    prefix = sys.argv[1]

    df = read_csvs(prefix)
    #print(df)
    compare_data(None,df)


    #year, plot, ecoregion, species, cohortAge, cohortBiomass

