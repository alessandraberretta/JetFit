import pandas as pd
from decimal import Decimal


def manipulate_txt(file, ergcm2s=False):

    with open(file, 'r') as file_input:
        read = file_input.readlines()

    last_row = -1
    found = False

    for idx, i in enumerate(read):
        if ergcm2s:
            if "PC" in i and not "PCUL" in i:
                first_row = idx+1
                found = True
                continue
            if "NO" in i:
                if found:
                    last_row = idx
        else:
            if "xrtpcflux" in i:
                first_row = idx+1
                found = True
            if "NO" in i:
                if found:
                    last_row = idx

    rows = read[first_row:last_row]

    if len(rows) <= 10:
        print("Def csv file not produced... len < 10")
        return ""

    with open(file[:file.rfind('.')]+".csv", 'w') as file_output:
        file_output.write(
            'Times,TimeBnds,TimeBndsNeg,Fluxes,FluxErrs,FluxErrsNeg\n')
        for elm in rows:
            file_output.write(elm.replace('\t', ','))

    df = pd.read_csv(file_output.name)

    del df["TimeBnds"]
    del df["TimeBndsNeg"]
    del df["FluxErrsNeg"]

    df["Fluxes"] = df["Fluxes"]
    df["FluxErrs"] = df["FluxErrs"]
    df["Times"] = ['%.6E' % Decimal(x) for x in df['Times']]
    df["Fluxes"] = ['%.6E' % Decimal(y) for y in df['Fluxes']]
    df["FluxErrs"] = ['%.6E' % Decimal(z) for z in df['FluxErrs']]
    df.insert(loc=1, column='Freqs', value='%.6E' % Decimal(2E+18))

    df.to_csv(file_output.name[:file_output.name.rfind(
        '.')] + "_def.csv", index=False)

    return file_output.name[:file_output.name.rfind(
        '.')] + "_def.csv"
