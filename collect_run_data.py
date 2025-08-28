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
    return total_df
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
                                Coalesce(SIM.SpeciesName,d.SPECIES_SYMBOL) SPECIES_SYMBOL,
                                Coalesce(SIM.CohortAge,d.ageclass) ageclass,
                                Coalesce(SIM.SIM_YEAR, d.SIM_YEAR) SIM_YEAR,
                                min_df.min_measyear,
                                d.MEASYEAR,
                                SIM.CohortBiomass,
                                d.drybio_ag
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
                        sim.CohortAge = d.ageclass AND
                        sim.SIM_YEAR  = d.SIM_YEAR 

                        """)

    mse = pdsql.sqldf(""" SELECT sum(pow(coalesce(drybio_ag,0) - coalesce(cohortbiomass,0),2)) FROM h""")
    print(mse)
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

