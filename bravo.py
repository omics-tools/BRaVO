#!/usr/bin/env python
# -*- coding:utf-8 -*-

'''
BRaVO: Bacterial Relative Variation Outliers
Created on Mar. 2019
Verion 0.0.1b
Lisense MIT Lisense
For more details, please check the documentation at https://github.com/omics-tools/BRaVO
'''

import os.path
import pkg_resources
from collections import Counter
import numpy as np
import pandas as pd
from scipy import stats
from skbio.stats.composition import multiplicative_replacement
from statsmodels.stats.multitest import fdrcorrection_twostage
from skbio.stats.composition import clr
import argparse
from argparse import RawTextHelpFormatter
from itertools import chain

def main():

    desc="""
    ==============================================================

                            BRaVO
             Bacterial Relative Variation Outliers
                      > Version: 0.0.1b

    =============================================================="""

    parser = argparse.ArgumentParser(usage='Usage: python {} -t count_table.txt -g group_label.txt'.format(__file__),description=desc,formatter_class=RawTextHelpFormatter)
    parser.add_argument('-t',action='store',dest='count_table', help='A taxonomic abundance table (count table)',default="false",required=True)
    parser.add_argument('-g',action='store',dest='group_label', help='Group labels',default="false",required=True)
    parser.add_argument('-a',action='store',dest='alpha_level', help='Alpha level at which to control false discoveries (default:0.05)',default=0.05)
    parser.add_argument('-out_dir',action='store',dest='output',help='Output directory (default: The directory in the input file)',default='')
    parser.add_argument('-v',action='version',version='Version : 0.0.1b')
    args = parser.parse_args()

    class colorset:
        RED = '\033[31m'
        BLUE = '\033[34m'
        UNDERLINE = '\033[4m'
        END = '\033[0m'

    #Check for requirement packages
    requirement_pkg=["numpy","scipy","pandas","scikit-bio","statsmodels"]
    installed_pkg=[dist.project_name for dist in pkg_resources.working_set]
    missing_pkg = [i for i in requirement_pkg if i not in installed_pkg]
    if len(missing_pkg)>0:
        print colorset.UNDERLINE+colorset.RED+"\nENVIROMENTAL ERROR"+colorset.END
        print "\nWe could not found {0}.".format(missing_pkg)
        print"\nBefore running this script, you need to install it on your Python enviroment."
        print "\n(e.g., $ pip install {0})".format(" ".join(missing_pkg))
        sys.exit()

    #Check for input data
    group_label=pd.read_csv(os.path.abspath(args.group_label))
    label_class = Counter(group_label["Group"]).keys()
    count_table = pd.read_csv(os.path.abspath(args.count_table),index_col=0)
    for col in ["Group","Sample"]:
        if col not in group_label.columns:
            print colorset.UNDERLINE+colorset.RED+"\nINPUT ERROR"+colorset.END
            print "\nYour group label is an invalid format.\nWe could found {0} in your group label."
            sys.exit()
    missing_samples=[i for i in group_label['Sample'] if i not in count_table.columns]
    if len(missing_samples)>0:
        print colorset.UNDERLINE+colorset.RED+"\nINPUT ERROR"+colorset.END
        print "\nYour taxonomic abundance table is an invalid format.\nWe could not found sample [{0}] in the table.".format(", ".join(missing_samples))
        sys.exit()
    if len(count_table.columns) != len(group_label):
        print colorset.UNDERLINE+colorset.RED+"\nINPUT ERROR"+colorset.END
        print "The number of samples between the group label and the count table is different."
        sys.exit()
    print colorset.UNDERLINE+colorset.BLUE+"\nINPUT DATA"+colorset.END
    for i in sorted(Counter(group_label["Group"]).items(), key=lambda x: x[0]):
        print "\nGroup: {0} ({1} samples)".format(i[0],i[1])
    print "\n"

    #Output setting
    if len(args.output) == 0:
        root, ext = os.path.splitext(os.path.abspath(args.count_table))
    else:
        out_dir = os.path.abspath(args.output)
        root = out_dir +"/" + os.path.splitext(os.path.basename(args.count_table))[0]

    #Detect outliers
    def outliers(df_dat):
        df_T_array = df_dat.T.values
        df_mr=pd.DataFrame(multiplicative_replacement(df_T_array).T,index=df_dat.index,columns=df_dat.columns)
        df_mr_array = df_mr.T.values
        df_mr_clr= pd.DataFrame(clr(df_mr_array).T,index=df_mr.index, columns=df_mr.columns)
        rv_list=[]
        for i in df_mr_clr.index:
            rv_list.append(float(np.std(df_mr_clr.T[i]))/abs(np.mean(df_mr_clr.T[i])))
        p_dic={}
        rv_niqr = stats.iqr(rv_list)*0.7413
        for k,v in enumerate(df_mr_clr.index):
            robust_z = float(rv_list[k]-np.median(rv_list))/rv_niqr
            p_dic[v] = stats.norm.sf(abs(robust_z))*2
        p_dic_new = {"name":[],"p":[]}
        for i in df_mr.index:
            p_dic_new["name"].append(i)
            p_dic_new["p"].append(p_dic[i])
        corrected_pvals=fdrcorrection_twostage(p_dic_new["p"],method='bky',alpha=args.alpha_level)[1]
        outliers_list=fdrcorrection_twostage(p_dic_new["p"],method='bky',alpha=args.alpha_level)[0]
        return pd.DataFrame({'Feature':df_mr.index,'Outlier':outliers_list, 'Adjusted-Pvalue':corrected_pvals,'RV':rv_list}).sort_values('Adjusted-Pvalue')[['Feature','RV','Adjusted-Pvalue','Outlier']]

    if len(label_class) == 1:
        print "Detecting outliers ....\n"
        print "The FDR threshold (alpha) = {0}\n".format(args.alpha_level)
        print "Output File: {0}\n".format(root+"_bravo_outliers.csv")
        res = outliers(count_table)
        if len(res[res["Outlier"] == True]['Feature']) == 0:
            print "Any outliers were not found.\n"
            sys.exit()
        res[res["Outlier"] == True]['Feature','RV','Adjusted-Pvalue'].read_csv(root+"_bravo_outliers.csv",index=False)

    elif len(label_class) >= 2:
        print "Detecting outliers ....\n"
        print "The FDR threshold (alpha) = {0}\n".format(args.alpha_level)
        print "Output File: {0}\n".format(root+"_bravo_outliers.csv")
        group_sample = {k:[] for k in label_class}
        for i in group_label.iterrows():
            group_sample[i[1]["Group"]].append(i[1]["Sample"])
        group_var=list(chain.from_iterable([list(outliers(count_table[group_sample[i]])[outliers(count_table[group_sample[i]])["Outlier"] == True]["Feature"]) for i in label_class]))
        var_dup = [x for x in set(group_var) if group_var.count(x) > 1]
        res = outliers(count_table)
        if len(res[res["Outlier"] == True]['Feature']) == 0:
            print "Any outliers were not found."
            sys.exit()
        res = res[res['Outlier'] == True][['Feature','RV','Adjusted-Pvalue']]
        res = res[res['Feature'].isin(var_dup) == False].to_csv(root+"_bravo_outliers.csv",index=False)
        print "Finish!"

if __name__ == "__main__":
    print """
==============================================================

                            BRaVO
             Bacterial Relative Variation Outliers
                      > Version: 0.0.1b

=============================================================="""
    main()
