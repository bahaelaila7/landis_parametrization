from collections import namedtuple, OrderedDict
from typing import List
import subprocess
from dataclasses import dataclass
import sys
import numpy as np
import pandas as pd
from PIL import Image
RNG = np.random.default_rng(seed=1)

SpeciesData =namedtuple('SpeciesData',['species_symbol', 'longevity', 'sex_mat', 'seed_disp_eff', 'seed_disp_max', 'reprod_prob', 'sprout_age_min', 'sprout_age_max', 'fire_regen']) 

def generate_species_data(species_symbol):
    longevity =RNG.integers(2,61)*5
    seed_disp_max = 10*RNG.integers(100)
    return SpeciesData(species_symbol = species_symbol, 
                       longevity = longevity,
                       sex_mat = RNG.integers(min(41,longevity)), 
                       seed_disp_eff = RNG.integers(1,int(0.5*(seed_disp_max/5)))*5, 
                       seed_disp_max = seed_disp_max, 
                       reprod_prob = 1,
                       sprout_age_min = 0, 
                       sprout_age_max= min(80,longevity), 
                       fire_regen = 'none')
def generate_species_params(species_symbol):
    return SpeciesParams(SpeciesCode= species_symbol,
                         LeafLongevity= RNG.uniform(1,10),
                         WoodDecayRate= RNG.uniform()*0.4,
                         MortalityCurve=RNG.uniform(5,25),
                         GrowthCurve = RNG.uniform(),
                         LeafLignin = RNG.uniform()*0.4,
                         ShadeTolerance = 1+RNG.integers(5), 
                         FireTolerance = 1+RNG.integers(5),
                         )
def generate_spp_params(species_symbol, ecoregion, year = 0):
    return SppParams (Year = year,
                      EcoregionName= ecoregion,
                      SpeciesCode = species_symbol, 
                      ProbEstablish= 1,
                      ProbMortality = 0,
                      ANPPmax = 400,
                      BiomassMax=26000)



@dataclass
class SpeciesDataClass:
    species_symbol: str
    longevity: int
    sex_mat: int
    seed_disp_eff: int
    seed_disp_max: int
    reprod_prob: float
    sprout_age_min: int
    sprout_age_max: int
    fire_regen: str

@dataclass
class FIAPlotMeasurementDataClass:
    species_symbol: str
    ageclass: int
    drybio_ag: float

@dataclass
class FIAPlotMeasurementClass:
    invyr: int
    measyear: int
    measmon: int
    measday: int
    private_owner: bool
    data: List[FIAPlotMeasurementDataClass]

@dataclass
class FIAPlotClass:
    statecd: int
    unitcd: int
    countycd: int
    plot: int
    lon: str
    lat: str
    epa_l4: str
    measurements: List[FIAPlotMeasurementClass]

FIAPlot = namedtuple('FIAPlot', ['STATECD', 'UNITCD', 'COUNTYCD', 'PLOT', 'LON','LAT', 'EPA_L4','MEASUREMENTS'])
FIAMeasurement = namedtuple('FIAMeasurement', ['INVYR','MEASYEAR','MEASMON','MEASDAY','data'])
FIAMeasurementData = namedtuple('FIAMeasurementData',['species_symbol','age','drybio_ag'])


def load_plots_data(csv_file):
    df = pd.read_csv(csv_file)
    df = df.sort_values(by=['invyr','measyear','measmon','measday','statecd','unitcd','countycd','plot','species_symbol','ageclass','drybio_ag'])
    i = 0
    d = dict()
    plot_list = []
    for _,row in df.iterrows():
        if np.isnan(row['drybio_ag']) or row['drybio_ag'] == 0:
            continue
        plot_id = (row['statecd'],row['unitcd'],row['countycd'], row['plot'])
        if plot_id not in d:
            d[plot_id] = i
            idx = i 
            i+=1
            plot_list.append(FIAPlotClass(statecd=row['statecd'],unitcd = row['unitcd'],countycd=row['countycd'],plot=['plot'], lon = row['lon'], lat=['lat'], epa_l4= None, measurements = list()))
        else:
            idx = d[plot_id]
        plot_obj = plot_list[idx]
        if len(plot_obj.measurements) == 0 or row['invyr'] != plot_obj.measurements[-1].invyr:
            m = FIAPlotMeasurementClass(invyr = row['invyr'], measyear=row['measyear'], measmon = row['measmon'], measday=row['measday'], private_owner = True, data = list())
            plot_obj.measurements.append(m)
        plot_obj.measurements[-1].data.append(FIAPlotMeasurementDataClass(species_symbol = row['species_symbol'], ageclass = row['ageclass'], drybio_ag = row['drybio_ag']))

    return plot_list

            


def generate_initial_communities_and_ecoregions(prefix, plots):
    mapcodes = np.arange(len(plots)).reshape(1, -1).astype(np.uint8)
    img = Image.fromarray(mapcodes)
    img.save(f'{prefix}/initial-communities.tif', format='TIFF')
    ecoregions_dict = OrderedDict()
    eco_count = 0
    eco_name_offset = 1000
    ecoregions_dict['Inactive'] = (0,eco_name_offset,False)
    ecoregions = np.full(len(plots),0).astype(np.uint8)
    with open(f"{prefix}/biomass-succession_InitialCommunities.csv", 'w') as f:
        f.write('MapCode,SpeciesName,CohortAge,CohortBiomass\n')
        for plot_id, plot in enumerate(plots,start = 1):
            if plot.epa_l4 not in ecoregions_dict:
                eco_count += 1
                ecoregions_dict[plot.epa_l4] = (eco_count, eco_count+eco_name_offset, True)
                eco_idx = eco_count
            else:
                eco_idx, _, _ = ecoregions_dict[plot.epa_l4]
            ecoregions[plot_id-1] = eco_idx
            for m in plot.measurements[0].data:
                f.write(f'{plot_id},{m.species_symbol},{m.ageclass},{int(m.drybio_ag)}\n')
    with open(f'{prefix}/ecoregions.txt', 'w') as f:
        f.write('''LandisData	"Ecoregions"		

>> 	Map Code		
>>Active	Name	Description
>>--------------------------------------------------------			
''')
        for eco_desc, (eco_idx, eco_name, eco_active) in ecoregions_dict.items():
            f.write(f"{'yes' if eco_active else 'no'}  {eco_idx}  {eco_name}  {eco_desc}\n")
    print(ecoregions_dict)
    img = Image.fromarray(ecoregions.reshape(1,-1))
    img.save(f'{prefix}/ecoregions.tif', format='TIFF')
    return ecoregions_dict






SPP_HEADER = 'Year,EcoregionName,SpeciesCode,ProbEstablish,ProbMortality,ANPPmax,BiomassMax'
def generate_spp_file(prefix, spps):
    with open(f'{prefix}/SppEcoregionData.csv','w') as f:
        f.write(f'''{SPP_HEADER}
{'\n'.join(
','.join(str(x) for x in [r.Year, r.EcoregionName, r.SpeciesCode, r.ProbEstablish, r.ProbMortality, r.ANPPmax, r.BiomassMax])
for r in spps
)
}
''')

def generate_core_species_file(prefix, rs):
    s= f'''
    LandisData  Species

>>                      Sexual  Seed Disperal Dist  Vegetative   Sprout Age  Post-Fire
>> Name      Longevity  Maturity  Effective  Maximum  Reprod Prob  Min   Max   Regen
>> ----      ---------  --------  ---------  -------  -----------  ----------  --------
{'\n'.join(
'\t'.join(str(x) for x in [r.species_symbol, r.longevity, r.sex_mat, r.seed_disp_eff, r.seed_disp_max, r.reprod_prob, r.sprout_age_min, r.sprout_age_max, r.fire_regen])
for r in rs
)
}
'''
    with open(f'{prefix}/Core_species_data.txt''','w') as f:
        f.write(s)


def get_core_species_params():
    rs_txt='''QUFA	185	25	20	50	1	1	10	resprout
    QUST	300	25	20	50	1	1	60	resprout
    QUAL	400	20	30	100	1	1	80	resprout
    LIST2	160	25	60	183	1	1	50	resprout
    FAGR	300	40	20	50	1	1	100	resprout
    JUVI	200	10	100	1000	0	0	0	none
    QUNI	140	20	20	50	1	1	40	resprout
    QUVE	175	20	20	50	1	1	50	resprout
    PRSE2	200	10	100	1000	1	1	30	resprout
    ULAL	160	15	50	100	1	1	40	resprout
    LITU	300	20	30	150	0	0	0	none
    MORU2	100	10	50	300	1	1	50	resprout
    FRAM2	200	20	50	100	1	1	60	resprout
    CAAL27	500	25	20	50	1	1	150	resprout
    DIVI5	100	6	50	300	1	1	50	resprout
    CAGL8	250	25	20	50	1	1	100	resprout
    OSVI	105	20	30	60	1	1	40	resprout
    CECA4	90	5	10	30	1	1	20	resprout
    NYSY	650	20	100	1000	1	1	150	resprout
    QUPA5	150	25	20	50	1	1	60	resprout
    PIEC2	175	20	60	100	1	1	30	resprout
    FRPE	150	20	50	100	1	1	50	resprout
    MAPO	135	25	20	100	1	1	75	resprout
    CELA	150	15	50	200	1	1	60	resprout
    SASAD	85	8	30	100	1	1	30	resprout
    SILAL3	100	10	50	200	1	1	50	resprout
    QUMA2	300	35	20	50	1	1	80	resprout
    OXAR	150	5	20	50	1	1	50	resprout
    ACRU	155	4	50	150	1	1	60	resprout
    QURU	250	25	30	80	1	1	50	resprout
    ILOP	125	10	50	500	1	1	80	resprout
    PIEL	200	10	60	90	0	0	0	none
    MAGR4	200	10	60	90	0	0	0	none
    TAAS	200	10	60	90	0	0	0	none
    QUIN	200	10	60	90	0	0	0	none
    NYBI	200	10	60	90	0	0	0	none
    COFL2	100	5	20	100	1	1	30	resprout
    PITA	200	10	60	90	0	0	0	none
    MEAZ	140	3	50	300	1	1	20	resprout
    QUCO2	190	20	20	50	1	1	40	resprout
    QUPR2	350	20	20	50	1	1	80	resprout
    ACBA3	200	30	40	100	1	1	50	resprout
    SAAL5	100	10	50	200	1	1	100	resprout
    PIVI2	130	5	50	80	0	0	0	none
    MAMA2	120	10	30	100	0	1	20	resprout
    MAVI2	90	10	50	200	1	1	50	resprout
    CAOV2	300	40	20	50	1	1	100	resprout
    MAAC	150	30	50	200	0	1	50	none
    CACO15	200	30	20	50	1	1	60	resprout
    ULRU	250	15	50	100	1	1	60	resprout
    QUMI	200	25	20	50	1	1	80	resprout
    ULAM	200	15	60	150	1	1	50	resprout
    QUPH	175	20	20	50	1	1	40	resprout
    QUMA3	125	20	20	50	1	1	40	resprout
    JUNI	225	12	30	80	1	1	40	resprout
    CEOC	175	15	50	200	1	1	60	resprout
    CACA18	120	15	30	60	1	1	40	resprout
    PISE	105	10	60	90	0	0	0	serotiny
    GOLA	75	10	30	80	0	0	0	none
    QULA3	145	15	20	50	1	1	30	resprout
    CATE9	150	25	20	50	1	1	60	resprout
    CRATA	75	5	20	100	1	1	20	resprout
    ROPS	100	6	50	200	1	1	90	resprout
    AIAL	90	5	100	300	1	1	50	resprout
    ACSA2	190	11	60	150	1	1	40	resprout
    ACNE2	115	5	50	150	1	1	30	resprout
    PIPA2	300	30	20	50	1	1	5	resprout
    CAOV3	200	25	20	50	1	1	60	resprout
    QUSH	200	25	20	50	1	1	60	resprout
    QULA2	130	20	20	50	0	0	0	none
    QULY	300	25	20	50	1	1	80	resprout
    TIAM	145	15	30	80	1	1	50	resprout
    PRAM	90	3	20	50	1	1	20	resprout
    QUMU	300	20	20	50	1	1	80	resprout
    TIAMC	140	15	30	80	1	1	50	resprout
    PIST	300	20	60	200	0	0	0	none
    TSCA	400	25	40	200	0	0	0	none
    BENI	80	15	100	200	1	1	30	resprout
    SANI	50	5	200	1000	1	1	50	resprout
    CAIL2	300	20	30	80	1	1	60	resprout
    PRUNU	50	5	50	300	1	1	50	resprout
    MAAN3	70	5	30	100	1	1	50	resprout
    CALA21	300	40	20	50	1	1	100	resprout
    PLOC	250	25	100	300	1	1	60	resprout
    ALJU	100	3	50	200	1	1	20	resprout'''
    rs = []
    for x in rs_txt.split('\n'):
        rs.append(SpeciesData(*x.split()))

    return rs
SppParams = namedtuple('SppParams',SPP_HEADER.split(','))
def get_spp_params():
    spp_txt = '''0,1001,QUFA,0.58,0.05,494,26000
0,1001,QUST,0.7,0.02,639,26000
0,1001,QUAL,0.62,0.07,529,26000
0,1001,LIST2,0.58,0.06,299,26000
0,1001,FAGR,0.5,0.03,523,26000
0,1001,JUVI,0.65,0.06,161,26000
0,1001,QUNI,0.51,0.02,441,26000
0,1001,QUVE,0.82,0.06,144,26000
0,1001,PRSE2,0.87,0.09,188,26000
0,1001,ULAL,0.47,0.04,133,26000
0,1001,LITU,0.75,0.07,745,26000
0,1001,MORU2,0.57,0.02,332,26000
0,1001,FRAM2,0.6,0.07,136,26000
0,1001,CAAL27,0.85,0.04,559,26000
0,1001,DIVI5,0.25,0.03,390,26000
0,1001,CAGL8,0.26,0.06,297,26000
0,1001,OSVI,0.21,0.01,354,26000
0,1001,CECA4,0.78,0.08,180,26000
0,1001,NYSY,0.74,0.01,546,26000
0,1001,QUPA5,0.81,0.07,236,26000
0,1001,PIEC2,0.89,0.03,289,26000
0,1001,FRPE,0.76,0.08,229,26000
0,1001,MAPO,0.52,0.1,309,26000
0,1001,CELA,0.75,0.03,468,26000
0,1001,SASAD,0.28,0.06,391,26000
0,1001,SILAL3,0.65,0.06,476,26000
0,1001,QUMA2,0.3,0.06,447,26000
0,1001,OXAR,0.86,0.03,268,26000
0,1001,ACRU,0.57,0.1,276,26000
0,1001,QURU,0.49,0.05,637,26000
0,1001,ILOP,0.39,0.09,423,26000
0,1001,PIEL,0.74,0.07,608,26000
0,1001,MAGR4,0.74,0.07,608,26000
0,1001,TAAS,0.74,0.07,608,26000
0,1001,QUIN,0.74,0.07,608,26000
0,1001,NYBI,0.74,0.07,608,26000
0,1001,COFL2,0.52,0.04,214,26000
0,1001,PITA,0.6,0.08,129,26000
0,1001,MEAZ,0.21,0.05,389,26000
0,1001,QUCO2,0.63,0.09,246,26000
0,1001,QUPR2,0.63,0.06,373,26000
0,1001,ACBA3,0.63,0.09,321,26000
0,1001,SAAL5,0.86,0.07,440,26000
0,1001,PIVI2,0.68,0.08,609,26000
0,1001,MAMA2,0.45,0.06,614,26000
0,1001,MAVI2,0.51,0.1,169,26000
0,1001,CAOV2,0.69,0.07,457,26000
0,1001,MAAC,0.24,0.05,656,26000
0,1001,CACO15,0.67,0.06,422,26000
0,1001,ULRU,0.67,0.01,723,26000
0,1001,QUMI,0.35,0.04,191,26000
0,1001,ULAM,0.29,0.07,441,26000
0,1001,QUPH,0.42,0.04,523,26000
0,1001,QUMA3,0.45,0.07,651,26000
0,1001,JUNI,0.6,0.05,762,26000
0,1001,CEOC,0.51,0.02,514,26000
0,1001,CACA18,0.89,0.04,757,26000
0,1001,PISE,0.27,0.06,374,26000
0,1001,GOLA,0.35,0.06,143,26000
0,1001,QULA3,0.31,0.06,183,26000
0,1001,CATE9,0.66,0.07,533,26000
0,1001,CRATA,0.38,0.07,193,26000
0,1001,ROPS,0.53,0.05,274,26000
0,1001,AIAL,0.37,0.09,505,26000
0,1001,ACSA2,0.31,0.04,301,26000
0,1001,ACNE2,0.28,0.05,598,26000
0,1001,PIPA2,0.66,0.09,580,26000
0,1001,CAOV3,0.3,0.08,429,26000
0,1001,QUSH,0.34,0.07,640,26000
0,1001,QULA2,0.46,0.02,309,26000
0,1001,QULY,0.77,0.09,717,26000
0,1001,TIAM,0.27,0.07,675,26000
0,1001,PRAM,0.79,0.1,116,26000
0,1001,QUMU,0.27,0.02,206,26000
0,1001,TIAMC,0.88,0.09,264,26000
0,1001,PIST,0.53,0.02,706,26000
0,1001,TSCA,0.88,0.07,636,26000
0,1001,BENI,0.62,0.02,728,26000
0,1001,SANI,0.72,0.09,291,26000
0,1001,CAIL2,0.23,0.08,295,26000
0,1001,PRUNU,0.4,0.06,407,26000
0,1001,MAAN3,0.28,0.05,236,26000
0,1001,CALA21,0.41,0.02,193,26000
0,1001,PLOC,0.28,0.07,187,26000
0,1001,ALJU,0.42,0.05,260,26000'''
    spps = []
    for x in spp_txt.split('\n'):
        spps.append(SppParams(*x.split(',')))
    return spps

def generate_biomass_succession_file(prefix, ecoregions):
    s=f'''LandisData  "Biomass Succession"


>>------------------
>> REQUIRED INPUTS
>>------------------

Timestep  			5

SeedingAlgorithm  		WardSeedDispersal

InitialCommunities      	./biomass-succession_InitialCommunities.csv
InitialCommunitiesMap   	initial-communities.tif
ClimateConfigFile		./biomass-succession_ClimateGenerator.txt

CalibrateMode 		no


>>----------------------------
>> LIFE HISTORY PARAMETERS
>>----------------------------

MinRelativeBiomass
>> Shade	Percent Max Biomass
>> Class	by Ecoregions
>> ----------	--------------------	
            {
            '\n'.join(f"""
                {eco_name}        
	1	25%        
	2	45%     
	3	56%     
	4	70%     
	5	90% 	
        """ for _,(_,eco_name,is_active) in ecoregions.items() if is_active)
        }


SufficientLight
>> Spp Shade	Probability
>> Class	by Actual Shade
>> ----------	--------------------	
>>		0	1	2	3	4	5
	1	1.0	0.5	0.25	0.0	0.0	0.0
	2	1.0	1.0	0.5	0.25	0.0	0.0
	3	1.0	1.0	1.0	0.5	0.25	0.0
	4	1.0	1.0	1.0	1.0	0.5	0.25
	5	0.1	0.5	1.0	1.0	1.0	1.0


SpeciesDataFile		SpeciesData.csv

EcoregionParameters
>>	AET (mm)
{'\n'.join(f"{eco_name}	600" for _,(_,eco_name,is_active) in ecoregions.items() if is_active)}

SpeciesEcoregionDataFile   SppEcoregionData.csv 

FireReductionParameters
>>	Severity	WoodLitter	Litter	
>>	Fire		Reduct		Reduct	
	1		0.0		0.5	
	2		0.0		0.75	
	3		0.0		1.0	

HarvestReductionParameters
>>	Name		WoodLitter	Litter	Cohort		Cohort
>>			Reduct		Reduct	WoodRemoval	LeafRemoval
	MaxAgeClearcut	0.5		0.15	0.8		0.0
	PatchCutting	1.0		1.0	1.0		0.0
'''




    with open(f'{prefix}/biomass-succession.txt', 'w') as f:
        f.write(s)

def generate_scenario_file(prefix):
    s = '''LandisData  "Scenario"


>> ---------------------------------------------
>> DEFINING A SCENARIO FOR A SINGLE LANDIS-II RUN
>>----------------------------------------------

>>	1. Provide the Required Inputs
>>	2. Select ONE Succession Extension
>>	3. Select ONE OR MORE Disturbance Extensions (but only ONE harvest extension)
>>	4. Select ONE OR MORE (compatible) Output Extensions

>>	A selection is made active by uncommenting a line (ie, remove the >> symbols) 



>>-------------------
>> REQUIRED INPUTS
>>-------------------

Duration  	50

Species   	Core_species_data.txt

Ecoregions      ./ecoregions.txt
EcoregionsMap   ./ecoregions.tif

CellLength  	30 << meters, 100 x 100 m = 1 ha



>> -----------------------
>> SUCCESSION EXTENSIONS
>> -----------------------

>> 	Succession Extension     Initialization File
>> 	--------------------     -------------------
	"Biomass Succession"	biomass-succession.txt


>> --------------------------
>> DISTURBANCE EXTENSIONS
>> -------------------------

>> 	Disturbance Extension	Initialization File
>>	--------------------	-------------------


>>   DisturbancesRandomOrder  yes  	<< optional
                         		<< Commented (default) is "no"

>> ------------------------
>> OUTPUT EXTENSONS
>> ----------------------

>> 	Output Extension		Initialization File
>> 	----------------		-------------------
	"Output Biomass"		output_Biomass.txt
	"Output Biomass Community"	output_Biomass_community.txt
	"Output Cohort Statistics"	output_CohortStats.txt



 RandomNumberSeed  1337  << optional parameter; uncomment for reproducibilty tests
                          << Commented (default) is a RandomNumberSeed generated using the current time
'''
    with open(f'{prefix}/scenario.txt', 'w') as f:
        f.write(s)


SPECIES_PARAMS_HEADER = "SpeciesCode,LeafLongevity,WoodDecayRate,MortalityCurve,GrowthCurve,LeafLignin,ShadeTolerance,FireTolerance"
SpeciesParams = namedtuple('SpeciesParams', SPECIES_PARAMS_HEADER.split(','))

def generate_species_params_file(prefix, sps):
    with open(f'{prefix}/SpeciesData.csv','w') as f:
        f.write(f'''{SPECIES_PARAMS_HEADER}
{'\n'.join(
','.join(str(x) for x in [r.SpeciesCode,r.LeafLongevity,r.WoodDecayRate,r.MortalityCurve,r.GrowthCurve,r.LeafLignin,r.ShadeTolerance,r.FireTolerance])
for r in sps
)
}
''')
def get_species_params():
    sp_txt = '''QUFA,2.1,0.64,17.7,0.5,0.29,3,3
QUST,3.9,0.43,15.7,0.21,0.25,3,4
QUAL,3.2,0.24,6.8,0.94,0.21,1,4
LIST2,2.8,0.42,21.7,0.54,0.22,4,1
FAGR,1.5,0.43,11.4,0.97,0.18,4,1
JUVI,1.5,0.71,8.7,0.97,0.15,5,1
QUNI,1.2,0.65,5.8,0.88,0.17,1,1
QUVE,3.6,0.82,16.8,0.44,0.25,3,3
PRSE2,2.8,0.53,18.6,0.51,0.1,4,2
ULAL,3.1,0.28,5.3,0.88,0.12,1,3
LITU,1.1,0.7,15.2,0.45,0.11,4,1
MORU2,3.9,0.73,9.5,0.34,0.11,4,1
FRAM2,3.5,0.59,17.9,0.65,0.27,3,3
CAAL27,1.6,0.74,8.5,0.95,0.24,2,1
DIVI5,1.5,0.55,18.8,0.76,0.19,5,4
CAGL8,1.6,0.57,12.7,0.66,0.12,5,3
OSVI,1.9,0.5,23.7,0.28,0.2,3,1
CECA4,2.6,0.22,7.8,0.69,0.19,4,1
NYSY,2.3,0.28,11.8,0.99,0.13,1,2
QUPA5,1.9,0.22,7.3,0.31,0.19,4,1
PIEC2,2.8,0.65,23.5,0.61,0.18,3,5
FRPE,1.4,0.42,22.5,0.9,0.22,5,2
MAPO,1.9,0.56,10.2,0.79,0.23,4,4
CELA,2.1,0.84,18.2,0.76,0.11,5,1
SASAD,2.4,0.37,21.3,0.76,0.17,1,2
SILAL3,3.4,0.49,16.1,0.49,0.23,5,3
QUMA2,1.6,0.73,15.6,0.43,0.2,5,5
OXAR,2.5,0.36,9.8,0.85,0.27,2,2
ACRU,2.8,0.25,6.9,0.85,0.23,2,2
QURU,1.1,0.4,22.9,0.89,0.13,2,3
ILOP,2.8,0.31,23,0.93,0.11,5,1
PIEL,1.5,0.85,17.7,0.61,0.23,3,4
MAGR4,1.5,0.85,17.7,0.61,0.23,3,4
TAAS,1.5,0.85,17.7,0.61,0.23,3,4
QUIN,1.5,0.85,17.7,0.61,0.23,3,4
NYBI,1.5,0.85,17.7,0.61,0.23,3,4
COFL2,1.2,0.77,11.8,0.6,0.11,5,1
PITA,3.8,0.64,12,0.84,0.22,3,3
MEAZ,3.9,0.81,19.5,0.72,0.29,3,1
QUCO2,3.4,0.76,22.9,0.76,0.22,2,4
QUPR2,1.9,0.33,22.7,0.84,0.18,4,4
ACBA3,1.3,0.82,20.6,0.91,0.23,1,1
SAAL5,3.1,0.58,17.8,0.47,0.19,2,2
PIVI2,2.3,0.77,6.7,0.5,0.21,2,2
MAMA2,1.4,0.83,8.2,0.28,0.29,4,1
MAVI2,2.5,0.42,23,0.66,0.18,1,1
CAOV2,1.1,0.28,17.1,0.23,0.29,5,3
MAAC,3.7,0.36,5.2,0.57,0.28,5,1
CACO15,1.8,0.5,7,0.63,0.14,2,1
ULRU,3,0.77,18.3,0.43,0.11,1,2
QUMI,1.9,0.8,5.1,0.67,0.12,2,2
ULAM,2.6,0.2,8.2,0.22,0.1,3,1
QUPH,2.6,0.56,16,0.23,0.12,2,1
QUMA3,1.6,0.49,18.8,0.86,0.24,2,5
JUNI,3.9,0.36,18,0.49,0.11,5,2
CEOC,3.3,0.28,9.5,0.3,0.16,5,1
CACA18,3.8,0.44,19.2,0.62,0.27,5,4
PISE,3.7,0.86,9.7,0.82,0.1,3,1
GOLA,2.8,0.43,11.5,0.37,0.26,5,2
QULA3,3.8,0.56,19.9,0.7,0.16,1,2
CATE9,1.3,0.69,18,0.27,0.12,4,5
CRATA,1.6,0.45,22,0.24,0.24,1,3
ROPS,1.1,0.88,18.2,0.63,0.23,1,2
AIAL,2,0.87,16.4,0.63,0.28,5,4
ACSA2,2.2,0.38,6.9,0.71,0.25,4,5
ACNE2,1.8,0.55,12.4,0.78,0.26,4,2
PIPA2,3.5,0.41,10.3,0.98,0.16,4,2
CAOV3,2.1,0.4,9.9,0.61,0.14,3,2
QUSH,1.8,0.23,24.5,0.46,0.25,5,3
QULA2,2.6,0.63,12.9,0.84,0.26,4,2
QULY,1.4,0.55,22.8,0.42,0.3,3,3
TIAM,3.4,0.24,17.6,0.55,0.18,2,2
PRAM,1.2,0.4,20.9,0.26,0.17,2,3
QUMU,4,0.84,15.1,0.22,0.26,3,2
TIAMC,3.3,0.37,16.5,0.97,0.17,3,1
PIST,1.6,0.3,14.9,0.87,0.29,5,1
TSCA,1,0.54,8.9,0.76,0.27,5,1
BENI,3.4,0.89,19.4,0.53,0.19,2,1
SANI,3.1,0.37,10.6,0.34,0.25,4,2
CAIL2,3.2,0.67,5.5,0.33,0.25,2,2
PRUNU,3.3,0.73,17.9,0.4,0.12,4,2
MAAN3,1.2,0.37,8.5,0.64,0.28,4,1
CALA21,2.1,0.71,23.8,0.77,0.2,5,1
PLOC,1.3,0.46,24.1,0.73,0.27,1,1
ALJU,3.6,0.64,23.3,0.42,0.16,1,1'''
    sps = []
    for x in sp_txt.split('\n'):
        sps.append(SpeciesParams(*x.split(',')))
    return sps

def generate_prism_data(prefix, ecoregions):
    f = f'{prefix}/PRISM_data.csv'
    df = pd.read_csv(f)
    print(df)
    d = df['101']
    df = df.drop('101', axis=1)
    for _,(_,eco_name, is_active) in ecoregions.items():
        if is_active:
            df[str(eco_name)] = d
    df.to_csv(f, index=False)

def get_fl5_species():
    return ["PIEL",
"QUVI",
"QULA3",
"NYBI",
"PITA",
"PIPA2",
"TAAS",
"QUNI",
"GOLA",
"MAVI2",
"ACRU",
"LIST2",
"PRSE2",
"CAGL8",
"NYSY",
"QUHE2",
"QUST",
"PEBO",
"CAAL27",
"MAGR4"]

def replace_species(rs, species):
    import copy
    species = {x:"" for x in species}
    out = []
    non = []
    species_field = 'species_symbol' if hasattr(rs[0], 'species_symbol') else 'SpeciesCode'
    for r in rs:
        rspecies = getattr(r, species_field) 
        if rspecies in species:
            del species[rspecies]
            out.append(r)
        else:
            non.append(r)
    for s,r in zip(species, non):
        out.append(copy.replace(r, **{species_field:s}))
    return out


def generate_output_biomass_file(prefix, species):
    txt = f'''LandisData  "Output Biomass"

Timestep  5

MakeTable yes  << Optional parameter

Species {'\n\t'.join(species)}

MapNames  outputs/biomass/biomass-{{species}}-{{timestep}}.tif

DeadPools both
	  
MapNames  outputs/biomass/biomass-{{pool}}-{{timestep}}.tif'''
    with open(f'{prefix}/output_Biomass.txt','w') as f:
        f.write(txt)





if __name__ == '__main__':

    import sys

    prefix = f'./{sys.argv[1]}'
    print(subprocess.run(['cp', '-r', './template',prefix]))
    plots = load_plots_data('./data_fl5_plot_genus_sp_ba_age_agb_20.csv')
    print(len(plots))
    ecoregions = generate_initial_communities_and_ecoregions(prefix,plots)
    ##rs = get_core_species_params()
    ##rs = replace_species(rs,get_fl5_species())
    fl5_species = get_fl5_species()
    rs = [generate_species_data(sp) for sp in fl5_species]
    from pprint import pprint
    pprint(rs)
    generate_core_species_file(prefix, rs)
    ##sps = get_species_params()
    ##sps = replace_species(sps,get_fl5_species())
    sps = [generate_species_params(sp) for sp in fl5_species]
    pprint(sps)
    generate_species_params_file(prefix, sps)
    #spps = get_spp_params()
    #spps = replace_species(spps,get_fl5_species())
    print(ecoregions)
    spps = [generate_spp_params(sp,eco_name) for sp in fl5_species for _,(_,eco_name,is_active) in ecoregions.items() if is_active]
    pprint(spps)
    generate_spp_file(prefix, spps)
    generate_biomass_succession_file(prefix, ecoregions)
    generate_scenario_file(prefix)
    generate_prism_data(prefix, ecoregions)
    generate_output_biomass_file(prefix, get_fl5_species())


