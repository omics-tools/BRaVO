#!/usr/bin/env python
# -*- coding:utf-8 -*-

import numpy as np
import pandas as pd
from scipy import stats
from skbio.stats.composition import multiplicative_replacement
from statsmodels.stats.multitest import local_fdr
from statsmodels.stats.multitest import NullDistribution
from skbio.stats.composition import clr
import matplotlib.pyplot as plt
import numpy.random as rd
import argparse

parser = argparse.ArgumentParser(description='BRaVO\n v.0.0.1')
parser.add_argument('-i',action='store',dest='input_file', help='A taxonomic abundance data (count table).',default="false")
parser.add_argument('-g',action='store',dest='group', help='A taxonomic abundance data (count table).',default="false")
parser.add_argument('-out_dir',action='store',dest='output',help='Output directory (default: The directory in the input file)',default='')
parser.add_argument('-v',action='version',version='Version : 0.0.1b')
args = parser.parse_args()

print "\n"
print "============================================================================"
print "                                                                            "
print "                                  BRaVO                                     "
print "  　　　　　　　    Bacterial Relative Variation Outliers                     "
print "                          > Version: 0.0.1 (beta)                           "
print "                          > Date: 5, Dec., 2018                             "
print "                          > Written by K.I.                                 "
print "                                                                            "
print "============================================================================"

def outliers(df_dat):
  #入力データは列にサンプル名，行にOTUsとする
  #後ほど，Qiime等の出力に対応予定

  #1. ゼロ値置換（乗法置換法：δ=1/N**2; NはOTU数)

  #転置後，arrayへの変換
  df_T_array = df_dat.T.values
  df_mr=pd.DataFrame(multiplicative_replacement(df_T_array).T,index=df_dat.index,columns=df_dat.columns)

  #2. 有心対数比変換（構成比を実数として扱う）
  df_mr_array = df_mr.T.values
  df_mr_clr= pd.DataFrame(clr(df_mr_array).T,index=df_mr.index, columns=df_mr.columns)

  #3. 変動係数計算
  cv_list=[]
  for i in df_mr_clr.index:
    cv_list.append(np.std(df_mr_clr.T[i])/np.mean(df_mr_clr.T[i]))

  cv_dic={i:np.std(df_mr_clr.T[i])/np.mean(df_mr_clr.T[i]) for i in df_mr_clr.index}

  #4. Robust Z
  z_dic = {k:(v-np.median(cv_dic.values()))/stats.iqr(cv_dic.values()) for k,v in cv_dic.items()}

  #順番を維持したままのリストに変換
  z_dic_new = {"name":[],"z":[]}
  for i in df_mr.index:
    z_dic_new["name"].append(i)
    z_dic_new["z"].append(z_dic[i])

  #5. Local FDR
  z_array=np.array(z_dic_new["z"])
  null=NullDistribution(z_array)
  local_fdr_list = local_fdr(z_array,nbins=len(z_array)).tolist()

  return pd.DataFrame({k:v for k,v in zip(z_dic_new["name"],local_fdr_list)}.items(),columns=["Feature","FDR"]).sort_values('FDR')
