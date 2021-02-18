import pandas as pd
from decimal import Decimal

df = pd.read_csv(
    "/Users/alessandraberretta/JetFit/2012_ergcm2s/GRB120712A.csv")
del df["TimeBnds"]
del df["TimeBndsNeg"]
del df["FluxErrsNeg"]
df["Times"] = ['%.6E' % Decimal(x) for x in df['Times']]
df["Fluxes"] = ['%.6E' % Decimal(y) for y in df['Fluxes']]
df["FluxErrs"] = ['%.6E' % Decimal(z) for z in df['FluxErrs']]
df.insert(loc=1, column='Freqs', value='%.6E' % Decimal(1E+17))
df.to_csv("GRB_121024A_2.298_def.csv", index=False)
