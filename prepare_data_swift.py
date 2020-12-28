import pandas as pd
from decimal import Decimal

df = pd.read_csv("GRB190114A.csv")
del df["TimeBnds"]
del df["TimeBndsNeg"]
del df["FluxErrsNeg"]
df["Fluxes"] = 1000 * df["Fluxes"]
df["FluxErrs"] = 1000 * df["FluxErrs"]
df["Times"] = ['%.6E' % Decimal(x) for x in df['Times']]
df["Fluxes"] = ['%.6E' % Decimal(y) for y in df['Fluxes']]
df["FluxErrs"] = ['%.6E' % Decimal(z) for z in df['FluxErrs']]
df.insert(loc=1, column='Freqs', value='%.6E' % Decimal(2E+18))
df.to_csv("GRB190114A_new.csv", index=False)
