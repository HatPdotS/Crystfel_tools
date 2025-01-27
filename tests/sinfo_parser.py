import pandas as pd
import subprocess

def parse_sinfo_to_dataframe() -> pd.DataFrame:
    data = str(subprocess.check_output(["sinfo"]))
    data = data.replace("\\n", "\n")
    lines = [line.strip("'").strip() for line in data.split("\n") if line.strip()][1:]  
    lines = [line for line in lines if line]
    parsed_data = []
    for line in lines:
        parsed_data.append(line.split())
    df = pd.DataFrame(parsed_data, columns=["Partition", "Avail", "TimeLimit", "Nodes", "State", "Nodelist"])
    df["TimeLimit"] = pd.to_timedelta(df["TimeLimit"].str.replace("-", " days "))
    df["Nodes"] = df["Nodes"].astype(int)
    df["State"] = df["State"].astype("category")
    default_partition = [i for i in df.Partition if '*' in i][1]
    df.Partition = df.Partition.str.replace("*", "")
    df = df.loc[df.Partition != "admin"]
    print('I found the following partitions:')
    df_dense = df.groupby("Partition").first()
    df_dense.Nodes = df.groupby("Partition").Nodes.sum()
    print(df_dense)
    print('The default partition is:' , default_partition.strip('*'))   
    df_dense = df_dense.loc[df_dense.Nodes > 5]
    df_dense = df_dense.sort_values(by="TimeLimit", ascending=False)
    return list(df_dense.index)


print(parse_sinfo_to_dataframe())