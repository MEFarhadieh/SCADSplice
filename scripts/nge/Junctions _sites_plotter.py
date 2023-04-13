### scripts to visualise junctions on chromosome regions

import argparse
from collections import defaultdict
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.patches import Patch
import os
import pandas as pd
import seaborn as sns
from tqdm import tqdm
import warnings
warnings.filterwarnings("ignore")

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def plot_exon(locs,ax,h=0.25,offset = 0,bin_size=400,alpha=1,color="k",ecolor="k"):
  ax.add_patch(Rectangle((locs[0], -h + offset), locs[1] - locs[0], 2*h,edgecolor=ecolor,color=color,alpha=alpha,linewidth=0))

def load_gtf(gtf_file,filt_chr):
  gtf = pd.read_csv(gtf_file, names = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"], sep="\t")
  if 'NC_000023.11' in gtf["seqname"].unique():
    can_chrom = [x for x in gtf["seqname"].unique() if x.startswith("NC_")]
    name_dict = {x : "chr" + str(int(x.split("_")[1].split(".")[0])) for x in can_chrom}
    name_dict['NC_000023.11'] = "chrX"
    name_dict['NC_000024.10'] = "chrY"
    name_dict['NC_012920.1'] = "chrM"
    gtf["seqname"] = gtf["seqname"].map(name_dict)
    gtf = gtf[~gtf["seqname"].isna()]
    filt_chr = False
  try:
    gtf["gene_id"] = gtf["attribute"].str.split("gene_name").str[1].str.split(";").str[0].str.split('"').str[1]
  except:
    gtf["gene_id"] = gtf["attribute"].str.split("gene_id").str[1].str.split(";").str[0].str.split('"').str[1]
  gtf["transcript_id"] = gtf["attribute"].str.split("transcript_id").str[1].str.split(";").str[0].str.split('"').str[1]
  if filt_chr:
    chromosomes = gtf["seqname"].unique()
    chromosomes = [x for x in chromosomes if "_" not in x and not x.startswith("KN")]
    gtf = gtf[gtf["seqname"].isin(chromosomes)]
  gtf["chr_gene"] = gtf["seqname"] + gtf["gene_id"]
  return gtf 

def annotation_plot(gtf, domains, gene, end,outpath, suff):
  rev_dict = {"A" : "B","B" : "A"} 
  don_df = pd.read_csv("{}{}_{}_{}_coords.tsv".format(outpath,gene,end),sep="\t")
  don_df = don_df.rename(columns={"rank_acc" : "rank", "rank_don" : "rank"}).astype(int)
  if don_df["juncPosR1A"].nunique() == 1:
    let = "A"
  else:
    let = "B"
  shared_ends = list(don_df["juncPosR1" + rev_dict[let]])
  zoom = True      
  gene_gtf = gtf[gtf["gene_id"] == gene]
  gene_gtf = gene_gtf[gene_gtf["feature"].isin(["exon"])]
  legend = True
  if don_df["juncPosR1" + rev_dict[let]].nunique() > 1:
    colors = 5*[u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf',u'#1f77b4', u'#ff7f0e', u'#2ca02c', u'#d62728', u'#9467bd', u'#8c564b', u'#e377c2', u'#7f7f7f', u'#bcbd22', u'#17becf']
    someX, someY = 0.5, 0.5
    plt.figure(figsize=(12, 6))
    h = 1
    offset = 1
    currentAxis = plt.gca()
    count = 1
    y_labels = []
    arc_height = 10
    y_ticks = []
    arcs = False
    chromosome = gene_gtf["seqname"].iloc[0]
    gene_min_all = gene_gtf["start"].min()
    gene_max_all = gene_gtf["end"].max()
    legend_elements = []
    gene_domains = domains[(domains[1] == chromosome) & (domains[2] < gene_max_all) & (domains[3] > gene_min_all)]
    if gene_gtf["strand"].iloc[0] == "+":
      asc = True
    else:
      asc = False
    for ind, row in don_df.iterrows():
        plt.text(row["juncPosR1" + rev_dict[let]],gene_gtf["transcript_id"].nunique() + 1,row["rank"],horizontalalignment="center")
        plt.plot([row["juncPosR1" + rev_dict[let]],row["juncPosR1" + rev_dict[let]]],[ -0.5,gene_gtf["transcript_id"].nunique() + 0.5],color="red")
    plt.plot([row["juncPosR1" + let],row["juncPosR1" + let]],[ -0.5,gene_gtf["transcript_id"].nunique() + 0.5],color="blue")
    count = 0
    for transcript, trans_df in gene_gtf.groupby("transcript_id"):
      y_labels.append(transcript)
      y_ticks.append(count * offset)
      gene_min = trans_df["start"].min()
      gene_max = trans_df["end"].max()
      plt.plot([gene_min,gene_max],[offset * count,offset * count],color="k")

      for exons in set([tuple(x) for x in trans_df[["start","end"]].to_numpy()]):
        plot_exon(exons,currentAxis,offset = count * offset)
      i = 0
      for d in set([tuple(x) for x in gene_domains[[2,3,4]].to_numpy()]):
        plot_exon(d[:2],currentAxis,offset = count * offset,alpha = 0.4,color=colors[i],ecolor=None,h=0.5)
        legend_elements.append(Patch(facecolor=colors[i], edgecolor=None,label=d[2],alpha=0.4))
        i += 1
      count += 1
    if arcs:
      plt.ylim([-1,gene_gtf["transcript_id"].nunique() + arc_height + 1])
    else:
      plt.ylim([-1,gene_gtf["transcript_id"].nunique() + 2])
    currentAxis.ticklabel_format(useOffset=False,style="plain")
    if legend:
      currentAxis.legend(handles=legend_elements[:gene_domains.shape[0]],bbox_to_anchor=(1., 1.0))
    plt.yticks(y_ticks,y_labels)
    if zoom:
      print("end",end,"shared_ends",shared_ends)
      buff = max([abs(x - end) for x in shared_ends])/12
      plt.xlim(min([end] + shared_ends) - buff,max([end] + shared_ends) + buff)
    else:
      plt.xlim([gene_min_all,gene_max_all])
    plt.title("{} {} {}".format(gene,chromosome,end))
    plt.show()
    plt.savefig("{}{}{}_{}_don_ann{}.png".format(outpath,gene,end,suff),bbox_inches = "tight")
    plt.close()
    
def group_colors(groups):
    group_color_dict = {comp : col for comp, col in zip(groups,sns.color_palette("deep",len(groups)))}
    palette = sns.color_palette("deep",3)
    group_color_dict["control"] = palette[0]
    group_color_dict["ADM"] = palette[1]
    group_color_dict["ADH"] = palette[2]
    return group_color_dict
    
def type_colors():
    
    type_color_dict = {'Astro': '#e7969c',
             'EN': '#d6616b',
             'Endo': '#cedb9c',
             'IN': '#7b4173',
             'Micro': '#31a354',
             'Oligo': '#3182bd',
             'OPC': '#8c6d31',
             }    
    return cell_color_dict
    
def get_args():
  parser = argparse.ArgumentParser(description="make box plots for each donor and acceptor")
  parser.add_argument("--letter",help="which letter; donor: A, acceptor: B")
  parser.add_argument("--params",help="file of parameters to make plots for")
  parser.add_argument("--subfolder",help="subfolder to save in")
  parser.add_argument("--cell_lim",type=int,default=20,help="number of cells to limit ontologies to")
  args = parser.parse_args()
  return args

def dot_plot(don_df, let, let_dict, palette, onts, outpath, gene, don, tiss, rev_dict, domains,gtf, alpha, suff):
  don_df["ontology_rank"] = don_df["ontology"] + don_df["rank_" + let_dict[let]].astype(str)
  don_df["rank_count"] = don_df["ontology_rank"].map(don_df.groupby("ontology_rank")["numReads"].sum())
  ont_dict = {o : i for o, i in zip(onts,range(1,don_df["ontology"].nunique() + 1))}
  don_df["ont_num"] = don_df["ontology"].map(ont_dict)
  pdf = don_df.drop_duplicates("ontology_rank")
  pdf["rank_sum"] = pdf["ontology"].map(pdf.groupby("ontology")["rank_count"].sum())
  pdf["frac_rank"] = pdf["rank_count"] / pdf["rank_sum"]
  print("ranks",pdf["rank_" + let_dict[let]].value_counts())
  pdf["rank_" + let_dict[let]] = pdf["rank_" + let_dict[let]].astype(int)
  ann_dict = pd.Series(pdf.splice_ann.values,index=pdf["rank_" + let_dict[let]]).to_dict() 
  coords = pdf.drop_duplicates("rank_" + let_dict[let])[["rank_" + let_dict[let],"juncPosR1A","juncPosR1B"]]
  coords = coords[coords["rank_" + let_dict[let]].isin(range(1,int(coords["rank_" + let_dict[let]].max() + 1)))].sort_values("rank_" + let_dict[let])
  end = int(coords["juncPosR1" + let].iloc[0])
  coords.to_csv("{}{}_{}_{}_coords.tsv".format(outpath,gene,end),sep="\t",index=False)
  try:
    annotation_plot(gtf, domains, gene, end,outpath,suff)
  except Exception as e:
    print("tried annotation plot",e)
  print("ann_dict",ann_dict)
  g = sns.relplot(x="rank_" + let_dict[let], y="ont_num", size="frac_rank",
              sizes=(10, 400), alpha=alpha, palette=palette,hue="group",
              height=max(4,pdf["ontology"].nunique()*0.3), data=pdf)
  plt.yticks(range(1,don_df["ontology"].nunique() + 1),onts)
  plt.title("{}\n{} {} {} {}".format(,gene,tiss, don, let_dict[rev_dict[let]]))
  plt.show()
  plt.savefig("{}{}_{}_{}_{}_{}_dot{}.png".format(outpath, gene, don, tiss, let_dict[rev_dict[let]],suff),bbox_inches="tight")
  plt.close()
  return 0

def plot_df(df, let, cell_lim, outpath, gene, let_dict, palette, rev_dict, don,tiss, don_df, comp_sort,suff):
    expand_dict = {"don" : "donor","acc" : "acceptor"}
    df["ont_rank"] = df["ontology_" + let_dict[rev_dict[let]]].map(df.groupby("ontology_" + let_dict[rev_dict[let]])["avg_rank"].median())
    df["ont_75"] = df["ontology_" + let_dict[rev_dict[let]]].map(df.groupby("ontology_" + let_dict[rev_dict[let]])["avg_rank"].quantile(0.75))
    df["ont_25"] = df["ontology_" + let_dict[rev_dict[let]]].map(df.groupby("ontology_" + let_dict[rev_dict[let]])["avg_rank"].quantile(0.25))
    df["ont_max"] = df["ontology_" + let_dict[rev_dict[let]]].map(df.groupby("ontology_" + let_dict[rev_dict[let]])["avg_rank"].max())
    df["ont_min"] = df["ontology_" + let_dict[rev_dict[let]]].map(df.groupby("ontology_" + let_dict[rev_dict[let]])["avg_rank"].min())
    if comp_sort:
      df = df.sort_values(by=["group","free_annotation"])
    else:
      df = df.sort_values(["ont_rank","ont_75","ont_25","ont_max","ont_min"])
    ann_dict = pd.Series(df.splice_ann.values,index=df["rank_" + let_dict[let]]).to_dict() 
    df["sort_num"] = df["ontology"].map({x : y for x, y in zip(focus,range(len(focus)))})
    df = df.sort_values("sort_num")
    num_cells = list(df.drop_duplicates("ontology")["num_cells"])
    df["sum_reads"] = df["ontology"].map(df.groupby("ontology")["numReads"].sum())
    num_reads = list(df.drop_duplicates("ontology")["sum_reads"])
    medians = list(df.drop_duplicates("ontology")["ont_rank"])
    comp_dict = {"Endothelial" : u'#1f77b4', "Epithelial" : u'#ff7f0e', "Immune" : u'#2ca02c',"Stromal" :  u'#d62728'}
    if df.shape[0] > 0:
      fig,(ax2) = plt.subplots(1,figsize=(7,5))
      g = sns.boxplot(x="avg_rank", y="ontology",hue="group",dodge=False,
                       data=df,ax=ax2,palette = comp_dict)
      g.axes.set_ylabel('') 
      g.axes.set_xlabel('average {} rank per cell'.format(expand_dict[let_dict[let]])) 
      for i,artist in enumerate(ax2.artists):
          col = artist.get_facecolor()
          artist.set_edgecolor(col)
          artist.set_facecolor('None')
          for j in range(i*6,i*6+6):
              line = ax2.lines[j]
              line.set_color(col)
              line.set_mfc(col)
              line.set_mec(col)
      for legpatch in ax2.get_legend().get_patches():
          col = legpatch.get_facecolor()
          legpatch.set_edgecolor(col)
          legpatch.set_linewidth(10)
          legpatch.set_facecolor('None')
      for i in range(len(num_reads)):
        plt.text(df["avg_rank"].max() + (df["avg_rank"].max() - df["avg_rank"].min())/12,i, int(num_reads[i]))
        plt.scatter([medians[i]],[i],color = "k",s=20,zorder=100)
      g.legend_.remove()
      plt.title("{}\n{} {} {} {}\nmean: {:0.2f} median: {:0.2f}", don_df["avg_rank"].mean(), don_df["avg_rank"].median()))      
      plt.savefig("{}{}_{}_{}_{}_{}{}_{}.pdf".format(outpath, gene, don, tiss, let_dict[rev_dict[let]],suff,cell_lim),bbox_inches="tight",transparent=True)
      plt.savefig("{}{}_{}_{}_{}_{}{}_{}.png".format(outpath, gene, don, tiss, let_dict[rev_dict[let]],suff,cell_lim),bbox_inches="tight")
      print("saved","{}{}_{}_{}_{}_{}{}_{}.pdf".format(outpath, gene, don, tiss, let_dict[rev_dict[let]],suff,cell_lim))
      plt.show()
      plt.close()
      return df

def box(df, let, cell_lim, outpath, gene, let_dict, palette, rev_dict, domains, gtf, alpha, comp_sort, suff):
  for don, don_df in df.groupby("pos{}_group".format(let)):
    temp = plot_df(don_df.drop_duplicates("pos{}_cell".format(let)), let, cell_lim, outpath, gene, let_dict, palette, rev_dict, don,"all", don_df, comp_sort,suff)
    if not temp is None:
      if temp["ontology"].nunique() > 0:
        onts = list(temp.sort_values(["ont_rank","ont_75","ont_25","ont_max","ont_min"]).drop_duplicates("ontology")["ontology"].unique())
        onts.reverse()
        response = dot_plot(don_df, let, let_dict, palette,onts, outpath, gene, don, "all", rev_dict,domains,gtf, alpha, suff)
    for tiss, tiss_df in don_df.groupby("type"):
      temp = plot_df(tiss_df.drop_duplicates("pos{}_cell".format(let)), let, cell_lim, outpath, gene, let_dict, palette, rev_dict, don,tiss, don_df, comp_sort,suff)
      if not temp is None:
        if temp["ontology"].nunique() > 0:
          onts = list(temp.sort_values(["ont_rank","ont_75","ont_25","ont_max","ont_min"]).drop_duplicates("ontology")["ontology"].unique())
          onts.reverse()
          response = dot_plot(tiss_df, let, let_dict, palette,onts,outpath, gene, don, tiss, rev_dict, domains,gtf, alpha, suff)
          if type(response) != int:
            print("ERROR")
            return response

df = pd.read_parquet("../SpliZ_values/{}_sym_SVD_normdonor_S_0.1_z_0.0_b_5.pq",columns=["type","group","sub_type","geneR1A_uniq","posA_group","posB_group","cell", "juncPosR1A", "juncPosR1B", "numReads", "splice_ann"])
outpath = "../plots/"

domains = pd.read_csv("~/SCADSplice/references/ucscGenePfam.txt",sep="\t",header=None)
gtf_file = "~/SCADSplice/references/SICILIAN_human_hg38_Refs/gtf_file/grch38_known_genes.gtf"
gtf = load_gtf(gtf_file,False)

cell_lim = 0
df["ontology"] = df["type"] + df["group"] + df["free_annotation"]
df["ontology_gene"] = df["ontology"] + df["geneR1A_uniq"]
df["ontology_don"] = df["ontology"] + df["posA_group"]
df["ontology_acc"] = df["ontology"] + df["posB_group"]
df["num_cells"] = df["ontology_gene"].map(df.groupby("ontology_gene")["cell"].nunique())
df = df[df["num_cells"] > cell_lim]
df["num_cell_ont"] = df["ontology_gene"].map(df.groupby("ontology_gene")["cell"].nunique())
df = df[df["num_cell_ont"] > cell_lim]

params = "../files/CLU.tsv"
letter = "B"
focus = ['control Astro',
 'ADM Astro',
 'ADH Astro',
 'control EN',
 'ADM EN',
 'ADH EN',
 'conrol Endo',
 'ADM Endo',
 'ADH Endo',
 'conrol IN',
 'ADM IN',
 'ADH IN',
 'conrol Micro',
 'ADM Micro',
 'ADH Micro',
 'conrol Oligo',
 'ADM Oligo',
 'ADH Oligo',
 'control OPC',
 'ADM OPC',
 'ADH OPC']
focus.reverse()

params = pd.read_csv(params,sep="\t")
pos_dict = {"CLU" : [27609086],"PTGDS": [136979946]}
for gene, gene_df in params.groupby("gene"):
  pos_dict[gene] = list(gene_df["end"].unique())
print("params",params)

sub = True
const_color = False
comp_sort = False

let_dict = {"A" : "acc", "B" : "don"}
rev_dict = {"A" : "B", "B" : "A"}

exact = {"CLU" : ["juncPosR1B",[27610475,27614655],["control Astro","ADM Astro","ADH Astro"]}

groups = sorted([x for x in list(df["group"].unique()) if x != None])
alpha = 1
if const_color:
  palette = defaultdict(lambda : "black")
else:
  palette = group_colors(groups)

genes = list(pos_dict.keys())

for gene in tqdm(genes):
    print("gene",gene)
    genepath = outpath + gene + "/"
    if not os.path.exists(genepath):
      os.makedirs(genepath)
    gene_df = df[df["geneR1A_uniq"] == gene]
    
    print("size before sub",gene_df.shape)
    if sub:
      gene_df = gene_df[gene_df["juncPosR1{}".format(letter)].isin(pos_dict[gene])]
      print("size gene_df",gene_df.shape)
    suff = ""    
    if gene in exact.keys():
      gene_df = gene_df[gene_df[exact[gene][0]].isin(exact[gene][1])]
      gene_df = gene_df[gene_df["free_annotation"].isin(exact[gene][2])]
      suff = "_exact"
    for let in [letter]:
      if gene_df.shape[0] > 0: 
        gene_df = gene_df[gene_df["ontology"].isin(focus)]
        print("num ontologies",gene_df["ontology"].nunique())
        gene_df["sort_num"] = gene_df["ontology"].map({x : y for x, y in zip(focus,range(len(focus)))})
        gene_df = gene_df.sort_values("sort_num")
        gene_df["num_" + let_dict[let]] = gene_df["pos{}_group".format(let)].map(gene_df.groupby("pos{}_group".format(let))["pos{}_group".format(rev_dict[let])].nunique())
        gene_df = gene_df[gene_df["num_" + let_dict[let]] > 1]
        gene_df["pos{}_cell".format(let)] = gene_df["pos{}_group".format(let)] + gene_df["cell"]
        gene_df["rank_{}".format(let_dict[let])] = gene_df.groupby("pos{}_group".format(let))["juncPosR1{}".format(rev_dict[let])].rank(method="dense")
        gene_df["scaled_rank"] = gene_df["rank_" + let_dict[let]] * gene_df["numReads"]
        gene_df["num"] = gene_df["pos{}_cell".format(let)].map(gene_df.groupby("pos{}_cell".format(let))["scaled_rank"].sum())
        gene_df["denom"] = gene_df["pos{}_cell".format(let)].map(gene_df.groupby("pos{}_cell".format(let))["numReads"].sum())
        gene_df["avg_rank"] = gene_df["num"]/gene_df["denom"]
        print("gene size 2",gene_df.shape)
        box(gene_df, let, cell_lim, genepath, gene, let_dict, palette, rev_dict, domains, gtf,alpha, comp_sort, suff)
                  
