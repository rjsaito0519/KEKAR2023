import pandas as pd
import os
import numpy as np
import sys

# --- dfを書き込む関数  ----------------------------------------------------------
def write_data(file_path, df):
    if os.path.exists(file_path):
        data = pd.read_csv(file_path, index_col=0, header=0)

        for tmp_index in df.index:
            if tmp_index in data.index:
                data.loc[tmp_index, :] = df.loc[ tmp_index, : ]
            else:
                data = pd.concat([data, df.loc[tmp_index].to_frame().T])
        sorted_data = data.sort_index()
        sorted_data.to_csv(file_path)
    else:
        df.to_csv(file_path)

# +------------+
# |  pedestal  |
# +------------+
# -------------------------------------------------------------
def write_pedestal_indiv():

    file_path = "./csv_data/pedestal_indiv.csv"
    df = pd.read_csv("./tmp_data.csv")
    df.index = ["ch1", "ch2", "ch3", "ch4", "ch5", "ch6", "ch7", "ch8"]
    print("\n--- load dataframe ----------------------------------")
    print(df)
    print("-----------------------------------------------------")
    index1 = "{}ch{}_{:0=4}".format(df.at["ch1", "counter"], int(df.at["ch1", "ch"]), int(df.at["ch1", "run"]))
    index2 = "{}ch{}_{:0=4}".format(df.at["ch2", "counter"], int(df.at["ch2", "ch"]), int(df.at["ch2", "run"]))
    index3 = "{}ch{}_{:0=4}".format(df.at["ch3", "counter"], int(df.at["ch3", "ch"]), int(df.at["ch3", "run"]))
    index4 = "{}ch{}_{:0=4}".format(df.at["ch4", "counter"], int(df.at["ch4", "ch"]), int(df.at["ch4", "run"]))
    index5 = "{}ch{}_{:0=4}".format(df.at["ch5", "counter"], int(df.at["ch5", "ch"]), int(df.at["ch5", "run"]))
    index6 = "{}ch{}_{:0=4}".format(df.at["ch6", "counter"], int(df.at["ch6", "ch"]), int(df.at["ch6", "run"]))
    index7 = "{}ch{}_{:0=4}".format(df.at["ch7", "counter"], int(df.at["ch7", "ch"]), int(df.at["ch7", "run"]))
    index8 = "{}ch{}_{:0=4}".format(df.at["ch8", "counter"], int(df.at["ch8", "ch"]), int(df.at["ch8", "run"]))
    df.index = [index1, index2, index3, index4, index5, index6, index7, index8]

    write_data(file_path, df)

# -------------------------------------------------------------
def write_pedestal():

    file_path = "./csv_data/pedestal.csv"
    df = pd.read_csv("./tmp_data.csv")
    df.index = ["SAC", "BAC", "KVC"]
    print("\n--- load dataframe ----------------------------------")
    print(df)
    print("-----------------------------------------------------")
    index1 = "{}_{:0=4}".format(df.at["SAC", "counter"], int(df.at["SAC", "run"]))
    index2 = "{}_{:0=4}".format(df.at["BAC", "counter"], int(df.at["BAC", "run"]))
    index3 = "{}_{:0=4}".format(df.at["KVC", "counter"], int(df.at["KVC", "run"]))
    df.index = [index1, index2, index3]

    write_data(file_path, df)

# +--------------+
# |  efficiency  |
# +--------------+
# -------------------------------------------------------------
def write_efficiency():

    log1 = np.genfromtxt(
        'csv_data/eff_data1.csv',
        delimiter=',',
        skip_header=0
    )
    log2 = np.genfromtxt(
        'csv_data/eff_data2.csv',
        delimiter=',',
        skip_header=0
    )

    file_path1 = "./csv_data/efficiency1.csv"
    file_path2 = "./csv_data/efficiency2.csv"
    df = pd.read_csv("./tmp_data.csv")
    df.index = ["1"]

    if int(df.at["1", "run"]) in log1[:, 0]:
        df["x"] = log1[log1[:, 0] == int(df.at["1", "run"])][0][1]
        df["y"] = log1[log1[:, 0] == int(df.at["1", "run"])][0][2]
        print("\n--- load dataframe ----------------------------------")
        print(df)
        print("-----------------------------------------------------")
        index1 = "run{:0=4}".format(int(df.at["1", "run"]))
        df.index = [index1]

        write_data(file_path1, df)
    
    if int(df.at["1", "run"]) in log2[:, 0]:
        df["x"] = log2[log2[:, 0] == int(df.at["1", "run"])][0][1]
        df["y"] = log2[log2[:, 0] == int(df.at["1", "run"])][0][2]
        print("\n--- load dataframe ----------------------------------")
        print(df)
        print("-----------------------------------------------------")
        index1 = "run{:0=4}".format(int(df.at["1", "run"]))
        df.index = [index1]

        write_data(file_path2, df)

# +-------+
# |  SAC  |
# +-------+
# -------------------------------------------------------------
def write_SAC_one_photon_gain():

    file_path = "./csv_data/SAC_one_photon_gain_data.csv"

    df = pd.read_csv("./tmp_data.csv")
    df.index = ["ch1", "ch2", "ch3", "ch4", "ch5", "ch6"]
    index1 = "{:0=4}_{}".format(int(df.at["ch1", "run"]), int(df.at["ch1", "ch"]))
    index2 = "{:0=4}_{}".format(int(df.at["ch2", "run"]), int(df.at["ch2", "ch"]))
    index3 = "{:0=4}_{}".format(int(df.at["ch3", "run"]), int(df.at["ch3", "ch"]))
    index4 = "{:0=4}_{}".format(int(df.at["ch4", "run"]), int(df.at["ch4", "ch"]))
    index5 = "{:0=4}_{}".format(int(df.at["ch5", "run"]), int(df.at["ch5", "ch"]))
    index6 = "{:0=4}_{}".format(int(df.at["ch6", "run"]), int(df.at["ch6", "ch"]))
    df.index = [index1, index2, index3, index4, index5, index6]
    
    write_data(file_path, df)

def write_one_photon_gain():

    file_path = "./csv_data/SAC_one_photon_gain_data.csv"
    df = pd.read_csv("./tmp_data.csv")
    df.index = ["ch1", "ch2", "ch3", "ch4", "ch5", "ch6", "ch7", "ch8"]
    print("\n--- load dataframe ----------------------------------")
    print(df)
    print("-----------------------------------------------------")
    index1 = "{:.1f}_{}".format(df.at["ch1", "LED_vol"], int(df.at["ch1", "ch"]))
    index2 = "{:.1f}_{}".format(df.at["ch2", "LED_vol"], int(df.at["ch2", "ch"]))
    index3 = "{:.1f}_{}".format(df.at["ch3", "LED_vol"], int(df.at["ch3", "ch"]))
    index4 = "{:.1f}_{}".format(df.at["ch4", "LED_vol"], int(df.at["ch4", "ch"]))
    index5 = "{:.1f}_{}".format(df.at["ch5", "LED_vol"], int(df.at["ch5", "ch"]))
    index6 = "{:.1f}_{}".format(df.at["ch6", "LED_vol"], int(df.at["ch6", "ch"]))
    index7 = "{:.1f}_{}".format(df.at["ch7", "LED_vol"], int(df.at["ch7", "ch"]))
    index8 = "{:.1f}_{}".format(df.at["ch8", "LED_vol"], int(df.at["ch8", "ch"]))
    df.index = [index1, index2, index3, index4, index5, index6, index7, index8]

    write_data(file_path, df)

# def generate_tex_table(output_filename):
#     # テーブルのデータを定義
#     table_data = np.genfromtxt(
#                     "./tmp_data.csv",
#                     skip_header=0,
#                     delimiter=","
#                 )

#     table = "\\begin{table}\n"
#     table += "\\begin{tabular}{c|ccc}\n"
#     table += "\t & conditon1 & conditon2 & conditon3 \\\\\n"
#     table += "\t\\hline\n"
#     for row in table_data[1:]:
#         table += "\tPMT{:.0f} (LED: {:.1f} V) ".format(row[0], table_data[1][-1])
#         for i in range(3):
#             table += "& ${:.1f} \pm {:.1f}$ ".format( row[2*i+1], row[2*i+2] )
#         table += "\\\\\n"
#     table += "\t\\hline\n"
#     table += "\\end{tabular}\n"
#     table += "\\end{table}"

#     with open(output_filename, "w") as f:
#         f.write(table)

def write_SAC_summary():

    # params = np.genfromtxt(
    #     'csv_data/SAC.csv',
    #     delimiter=',',
    #     skip_header=1
    # )

    file_path = "./csv_data/SAC_summary.csv"
    df = pd.read_csv("./tmp_data.csv")
    df.index = ["sum"]
    # param = params[ params[:, 0] == int(df.at["sum", "run"]) ][0]
    # df["x_pos"] = param[1]
    # df["y_pos"] = param[2]
    # df["HV_cond"] = param[3]
    # df["Vth"] = param[4]
    print("\n--- load dataframe ----------------------------------")
    print(df)
    print("-----------------------------------------------------")
    index = "run{:0=4}".format(int(df.at["sum", "run"]))
    df.index = [index]
    
    write_data(file_path, df)

# +-------+
# |  BAC  |
# +-------+
# -------------------------------------------------------------
def write_BAC_one_photon_gain():

    params = np.genfromtxt(
        'csv_data/LED_BAC.csv',
        delimiter=',',
        skip_header=1
    )

    file_path = "./csv_data/BAC_one_photon_gain.csv"
    df = pd.read_csv("./tmp_data.csv")
    df.index = ["ch1", "ch2", "ch3", "ch4"]
    param = params[ params[:, 0] == int(df.at["ch1", "run"]) ][0]
    df["MPPC_ch"] = param[1]
    df["LED_vol"] = param[2]
    df["HV"] = param[3]
    df["use"] = [int(i) for i in str(int(param[4]))[1:]]
    print("\n--- load dataframe ----------------------------------")
    print(df)
    print("-----------------------------------------------------")
    index1 = "run{:0=4}_{}".format(int(df.at["ch1", "run"]), int(df.at["ch1", "ch"]))
    index2 = "run{:0=4}_{}".format(int(df.at["ch2", "run"]), int(df.at["ch2", "ch"]))
    index3 = "run{:0=4}_{}".format(int(df.at["ch3", "run"]), int(df.at["ch3", "ch"]))
    index4 = "run{:0=4}_{}".format(int(df.at["ch4", "run"]), int(df.at["ch4", "ch"]))
    df.index = [index1, index2, index3, index4]
    
    write_data(file_path, df)

if __name__ == '__main__':
    print()
    args = sys.argv
    if args[1] == "SAC":
        write_SAC_one_photon_gain()
    elif args[1] == "SAC_LED":
        write_one_photon_gain()
    elif args[1] == "pedestal":
        write_pedestal()
    elif args[1] == "pedestal_indiv":
        write_pedestal_indiv()
    elif args[1] == "SAC_summary":
        write_SAC_summary()
    elif args[1] == "BAC_opg":
        write_BAC_one_photon_gain()
    elif args[1] == "efficiency":
        write_efficiency()

    print("write data")