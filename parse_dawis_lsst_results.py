import pandas as pd
import numpy as np

res = pd.read_excel('/home/ellien/LSST_ICL/analysis/out4/run1/dawis_lsst_results.xlsx')
res.drop([144])

T1 = ['Wavelet separation', 'Size separation', 'Surface brightness threshold']
T2 = ['ICL Luminosity', 'ICL+BCG Luminosity', 'ICL Fraction', 'ICL+BCG Fraction' ]
T3 = ['xy', 'yz', 'xz']
col = [ ('Simulation', '', ''), ('Cluster ID','','') ]
for t1 in T1:
    for t2 in T2:
        for t3 in T3:
            col.append((t1,t2,t3))

col_list = pd.MultiIndex.from_tuples(col)

r = pd.DataFrame( columns = col_list)
for i in range(0, 152, 3):
    cluster_res = res.iloc[i:i+3]
    name = cluster_res['name'].iloc[0]
    name = name.split('.')[0]
    name = name.split('_')[1]
    sim = cluster_res['dir'].iloc[0]
    print(name, sim)
    xy_res = cluster_res.iloc[0]
    yz_res = cluster_res.iloc[2]
    zx_res = cluster_res.iloc[1]
    nr = pd.DataFrame( [[ sim, name, \
                          xy_res['FICL_wavsep'], yz_res['FICL_wavsep'], zx_res['FICL_wavsep'], xy_res['FICL_bcgwavsep'], yz_res['FICL_bcgwavsep'], zx_res['FICL_bcgwavsep'], \
                          xy_res['fICL_wavsep'], yz_res['fICL_wavsep'], zx_res['fICL_wavsep'], xy_res['fICL_bcgwavsep'], yz_res['fICL_bcgwavsep'], zx_res['fICL_bcgwavsep'], \
                          xy_res['FICL_sizesep'], yz_res['FICL_sizesep'], zx_res['FICL_sizesep'], xy_res['FICL_bcgsizesep'], yz_res['FICL_bcgsizesep'], zx_res['FICL_bcgsizesep'], \
                          xy_res['fICL_sizesep'], yz_res['fICL_sizesep'], zx_res['fICL_sizesep'], xy_res['fICL_bcgsizesep'], yz_res['fICL_bcgsizesep'], zx_res['fICL_bcgsizesep'], \
                          xy_res['FICL_sbt'], yz_res['FICL_sbt'], zx_res['FICL_sbt'], 0, 0, 0, \
                          xy_res['fICL_sbt'], yz_res['fICL_sbt'], zx_res['fICL_sbt'], 0, 0, 0  ]],  columns = col_list )
    r = pd.concat( [ r, nr ], ignore_index = True )

r.to_excel('/home/ellien/LSST_ICL/analysis/out4/run1/parsed_dawis_lsst_results.xlsx')
